package clump;

import java.util.ArrayList;
import java.util.Arrays;

import dna.AminoAcid;
import align2.Shared;
import align2.Tools;

import stream.Read;

/**
 * A list of reads sharing a kmer.
 * @author Brian Bushnell
 * @date Nov 7, 2015
 *
 */
public class Clump extends ArrayList<Read> {
			
	public Clump(long kmer_){
		this(kmer_, 8);
	}

	public Clump(long kmer_, int size){
		super(size);
		kmer=kmer_;
	}

	public boolean add(Read r){
		long[] obj=(long[]) r.obj;
		assert(obj[0]==kmer);
		return super.add(r);
	}
	
	/** This will create a count consensus of the bases at each position in the cluster. */
	private int[][] count(final boolean qualityScores){
		int maxLeft=-1, maxRight=-1;
		for(Read r : this){
			long[] obj=(long[]) r.obj;
			int pos=(int)obj[1];
			maxLeft=Tools.max(maxLeft, pos);
			maxRight=Tools.max(maxRight, r.length()-pos);
		}
		final int width=maxLeft+maxRight;
//		assert(size()==1) : "\nleft="+maxLeft+", right="+maxRight+", width="+width+", "+k+"\n"+get(0).toFastq()+"\n"+get(size()-1).toFastq();
		
//		System.err.println("\n\n");
		final int[][] counts=new int[4][width];
		for(Read r : this){
			long[] obj=(long[]) r.obj;
			int pos=(int)obj[1];
			byte[] bases=r.bases, quals=r.quality;
//			System.err.println("pos="+pos+", maxLeft="+maxLeft);
			for(int cloc=0, rloc=maxLeft-pos; cloc<bases.length; cloc++, rloc++){
//				System.err.println("cloc="+cloc+"/"+bases.length+", rloc="+rloc+"/"+width);
				int x=AminoAcid.baseToNumber[bases[cloc]];
				if(x>-1){
					final int q;
					if(qualityScores){q=(quals==null ? 20 : quals[cloc]);}
					else{q=1;}
					counts[x][rloc]+=q;
				}
			}
		}
		return counts;
	}
	
	/*--------------------------------------------------------------*/
	
	public ArrayList<Read> condense_old(){
		Read r=makeSimpleConsensus();
		ArrayList<Read> list=new ArrayList<Read>();
		list.add(r);
		return list;
	}
	
	public Read makeSimpleConsensus(){
		if(size()==1){
			Read r=get(0);
			r.id=r.numericID+" size=1";
			return r;
		}
		final int[][] bcounts=baseCounts();
		final int[][] qcounts=qualityCounts();
		
		final byte[] bases=new byte[width], quals=new byte[width];
		for(int i=0; i<width; i++){
			final int x=getConsensusAtPosition(qcounts, i);
			final int y=getSecondHighest(qcounts, i);
			if(x<0){
//				System.err.println("q="+0+", x="+x+"; A="+counts[0][i]+", C="+counts[1][i]+", G="+counts[2][i]+", T="+counts[3][i]);
				bases[i]='N';
				quals[i]=0;
			}else{
				final long bmajor=bcounts[x][i];
				final long bminor=2*bcounts[x][i]-bcounts[0][i]-bcounts[1][i]-bcounts[2][i]-bcounts[3][i];
				final long bsecond=bcounts[y][i];
				
				final long qmajor=qcounts[x][i];
				final long qminor=2*qcounts[x][i]-qcounts[0][i]-qcounts[1][i]-qcounts[2][i]-qcounts[3][i];
				final long qsecond=qcounts[y][i];
				
				float bratio=bminor/(float)(bmajor+bminor);
				double phred=(bminor==0 ? Read.MAX_CALLED_QUALITY : -10*Math.log10(bratio));
				phred=Tools.min(phred, qmajor-qminor);
				assert(phred>=0 && phred<=127);
				byte q=(byte)Tools.mid(Read.MIN_CALLED_QUALITY, (long)Math.round(phred), Read.MAX_CALLED_QUALITY);
				bases[i]=AminoAcid.numberToBase[x];
				quals[i]=q;
			}
		}
		Read leftmost=this.get(0);
		Read r=new Read(bases, quals, 0, leftmost.numericID+" size="+size());
		//TODO: Attach the long pair, and make sure the kmer location is correct.
//		assert(false) : "\n"+r.toFastq()+"\nCheck kmer location.";
//		assert(size()==1) : "\n"+r.toFastq()+"\n"+get(0).toFastq()+"\n"+get(size()-1).toFastq()+"\n";
		return r;
	}
	
	/*--------------------------------------------------------------*/
	
	public ArrayList<Read> condense(){
		ArrayList<Read> list=makeConsensus();
		return list;
	}
	
	public ArrayList<Read> makeConsensus(){
		if(size()==1){
			Read r=get(0);
			r.id=r.numericID+" size=1";
			return this;
		}
		ArrayList<Read> list=new ArrayList<Read>(2);
		
		final Read consensus=consensusRead();
		
		//Find a position that differentiates two poulations of reads
		int pivot=findBestPivot();
		if(pivot<0){
			list.add(consensus);
			return list;
		}
		
		ArrayList<Clump> clumps=splitOnPivot(pivot);
		if(clumps==null || clumps.size()<2){
			list.add(consensus);
			return list;
		}
		
		for(Clump c : clumps){
			ArrayList<Read> temp=c.makeConsensus();
			list.addAll(temp);
		}
		return list;
	}
	
	/*--------------------------------------------------------------*/
	
	private ArrayList<Clump> splitOnPivot(final int pivot){
		//TODO
		return null;
//		ArrayList<Clump> list=new ArrayList<Clump>();
//		list.add(this);
//		return list;
	}
	
	private int findBestPivot(){
		final int[][] bcounts=baseCounts();
		final int[][] qcounts=qualityCounts();
		
		for(int i=0; i<width; i++){
			final int x=getConsensusAtPosition(bcounts, i);
			final int y=getSecondHighest(bcounts, i);
			if(x>=0){
				final long bmajor=bcounts[x][i];
				final long bsecond=bcounts[y][i];
				
				final long qmajor=qcounts[x][i];
				final long qsecond=qcounts[y][i];
				
//				float bratio=bsecond/(float)(bmajor+bsecond);
				float mult=5f;
				if(bsecond>1 && qsecond*mult>=qsecond+qmajor){//or bsecond*mult>=bmajor
					//do more tests
					return i;
				}
			}
		}
		return -1;
	}
	
	private float identity(byte[] consensus, Read r){
		long[] obj=(long[]) r.obj;
		int pos=(int)obj[1];
		byte[] bases=r.bases, quals=r.quality;
		int good=0, bad=0;
		for(int i=0, j=pos; i<bases.length; i++, j++){
			final byte a=bases[i], b=consensus[j];
			if(AminoAcid.isFullyDefined(a) && AminoAcid.isFullyDefined(b)){
				if(a==b){good++;}
				else{bad++;}
			}
		}
		return good==0 ? 0 : good/(float)(good+bad);
	}
	
	/*--------------------------------------------------------------*/
	
	private int getConsensusAtPosition(int[][] counts, int pos){
		int xMax=0;
		for(int x=1; x<4; x++){
//			System.err.println("x="+x+", max="+max+", Checking "+counts[x][pos]+" vs "+counts[x][max]);
			if(counts[x][pos]>counts[xMax][pos]){xMax=x;}
		}
//		assert(counts[max][pos]>=counts[0][pos]);
//		assert(counts[max][pos]>=counts[1][pos]);
//		assert(counts[max][pos]>=counts[2][pos]) : max+", "+counts[max][pos]+", ["+counts[0][pos]+", "+counts[1][pos]+", "+counts[2][pos]+", "+counts[3][pos]+"]";
//		assert(counts[max][pos]>=counts[3][pos]);
		return (counts[xMax][pos]>0 ? xMax : -1);
	}
	
	private int getSecondHighest(int[][] counts, int pos){
		int first=0;
		int second=1;
		if(counts[first][pos]<counts[second][pos]){
			first=1; second=0;
		}
		for(int x=2; x<4; x++){
//			System.err.println("x="+x+", max="+max+", Checking "+counts[x][pos]+" vs "+counts[x][max]);
			if(counts[x][pos]>counts[first][pos]){
				second=first;
				first=x;
			}else if(counts[x][pos]>counts[second][pos]){
				second=x;
			}
		}
//		assert(counts[max][pos]>=counts[0][pos]);
//		assert(counts[max][pos]>=counts[1][pos]);
//		assert(counts[max][pos]>=counts[2][pos]) : max+", "+counts[max][pos]+", ["+counts[0][pos]+", "+counts[1][pos]+", "+counts[2][pos]+", "+counts[3][pos]+"]";
//		assert(counts[max][pos]>=counts[3][pos]);
		
		return second; //may be actually 0.
		//return (counts[second][pos]>0 ? second : -1);
	}
	
	/*--------------------------------------------------------------*/
	
	public Read consensusRead(){
		if(consensusRead==null){
			consensusRead=makeSimpleConsensus();
		}
		return consensusRead;
	}
	
	public int width(){
		assert(width>=0);
		return width;
	}
	
	/*--------------------------------------------------------------*/
	
	private int[][] baseCounts(){
		if(baseCounts==null){
			baseCounts=count(false);
			int len=baseCounts[0].length;
			assert(width==-1 || width==len);
			width=len;
		}
		return baseCounts;
	}
	
	private int[][] qualityCounts(){
		if(qualityCounts==null){
			qualityCounts=count(false);
			int len=qualityCounts[0].length;
			assert(width==-1 || width==len);
			width=len;
		}
		return qualityCounts;
	}

	private void clearCounts(){
		baseCounts=qualityCounts=null;
	}
	
	private void clearConsensus(){
		consensusRead=null;
	}
	
	/*--------------------------------------------------------------*/
	
	public final long kmer;
	
	private int[][] baseCounts;
	private int[][] qualityCounts;
	private Read consensusRead;
	private int width=-1;
	
	private static final boolean countQuality=false;
	public static int k=31;
	private static final long serialVersionUID = 1L;
	
}
