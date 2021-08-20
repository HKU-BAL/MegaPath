package clump;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;

import jgi.Dedupe;
import align2.Shared;
import align2.Tools;

import kmer.KmerTableSet;
import kmer.Primes;

import stream.Read;

/**
 * @author Brian Bushnell
 * @date Nov 4, 2015
 *
 */
public class KmerComparator implements Comparator<Read> {
	
	public KmerComparator(int k_, long minDivisor_){
		k=k_;
		assert(k>0 && k<32);
		
		shift=2*k;
		shift2=shift-2;
		mask=~((-1L)<<shift);
		divisor=Primes.primeAtLeast(minDivisor_);
	}
	
	public void hashThreaded(ArrayList<Read> list){
		int threads=Shared.threads();
		ArrayList<HashThread> alt=new ArrayList<HashThread>(threads);
		for(int i=0; i<threads; i++){alt.add(new HashThread(i, threads, list));}
		for(HashThread ht : alt){ht.start();}
		
		/* Wait for threads to die */
		for(HashThread ht : alt){
			
			/* Wait for a thread to die */
			while(ht.getState()!=Thread.State.TERMINATED){
				try {
					ht.join();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
		}
	}
	
	public void hash(ArrayList<Read> list, KmerTableSet table, int minCount) {
		for(Read r : list){hash(r, table, minCount);}
	}
	
	private void hash(ArrayList<Read> list){
		for(Read r : list){hash(r);}
	}
	
	public long hash(Read r1, KmerTableSet table, int minCount){
		long[] kmers=new long[2];
		r1.obj=kmers;
		return fillLocalMax(r1, kmers, table, minCount);
	}
	
	private long hash(Read r1){
		long[] kmers=new long[2];
		r1.obj=kmers;
		return fillLocalMax(r1, kmers);
	}
	
	public void fuse(Read r1){
		Read r2=r1.mate;
		if(r2==null){return;}
		r1.mate=null;
		final int len1=r1.length(), len2=r2.length();
		int len=len1+len2+1;
		byte[] bases=new byte[len];
		for(int i=0; i<len1; i++){bases[i]=r1.bases[i];}
		bases[len1]='N';
		for(int i=0, j=len1+1; i<len2; i++){bases[j]=r2.bases[i];}
	}
	
	/* (non-Javadoc)
	 * @see java.util.Comparator#compare(java.lang.Object, java.lang.Object)
	 */
	@Override
	public int compare(Read a, Read b) {
		final long[] alist, blist;
		if(a.obj==null){
			alist=new long[2];
			a.obj=alist;
			fillLocalMax(a, alist);
		}else{alist=(long[])a.obj;}

		if(b.obj==null){
			blist=new long[2];
			b.obj=alist;
			fillLocalMax(b, blist);
		}else{blist=(long[])b.obj;}
		
		return compare(alist, blist);
	}

	
	/** Finds the global maximum */
	private long fillMax(Read r, long[] kmers){
		return fillMax(r, kmers, null, 0);
	}
	
	/** Finds the global maximum */
	public long fillMax(Read r, long[] kmers, KmerTableSet table, int minCount){
//		Arrays.fill(kmers, -1);
		kmers[0]=0;
		kmers[1]=k-1;
		final byte[] bases=r.bases;
		long kmer=0;
		long rkmer=0;
		int len=0;
		
		if(bases==null || bases.length<k){return -1;}
		
		long topMod=-1;
		boolean rcomp=false;
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			long x=Dedupe.baseToNumber[b];
			long x2=Dedupe.baseToComplementNumber[b];
			kmer=((kmer<<2)|x)&mask;
			rkmer=(rkmer>>>2)|(x2<<shift2);
			if(b=='N'){len=0;}else{len++;}
			if(len>=k){
				final long kmax=Tools.max(kmer, rkmer);
				final long mod=kmax%divisor;
				if(mod>topMod){
					if(minCount<2 || table.getCount(kmer, rkmer)>=minCount){
						topMod=mod;
						kmers[0]=kmax;
						kmers[1]=i;
						rcomp=(kmax!=kmer);
					}
				}
			}
		}
		rcomp&=rcompReads;
		
		if(topMod<0 && minCount>1){
			return fillMax(r, kmers, null, 0);
		}
		
//		r.id+=" "+kmers[1]+","+rcomp+","+(bases.length-kmers[1]+k-2);
		if(rcomp){
			r.reverseComplement();
			r.setSwapped(true);
			kmers[1]=bases.length-kmers[1]+k-2;
		}
		if(addName){r.id+=" "+kmers[1]+(rcomp ? ",t" : ",f")+","+kmers[0];}
		assert(kmers[0]>=0 && kmers[1]>=0) : Arrays.toString(kmers)+"\n"+r;
		return kmers[0];
	}
	
	/** Finds the highest local maximum */
	private long fillLocalMax(Read r, long[] kmers){
		return fillLocalMax(r, kmers, null, 0);
	}
	
	/** Finds the highest local maximum */
	public long fillLocalMax(Read r, long[] kmers, KmerTableSet table, int minCount){
		if(!LOCAL_MAX){return fillMax(r, kmers);}
		Arrays.fill(kmers, -1);//TODO: Note! 0, 0 can be detected and allowed to fall through.
		final byte[] bases=r.bases;
		long kmer=0;
		long rkmer=0;
		int len=0;
		
		if(bases==null || bases.length<k){return -1;}
		
		long topMod=-1;
		boolean rcomp=false;
		
		long mod1=-1, mod2=-1;
		long kmax1=-1;//, kmax2=-1;
		boolean rcomp1=false;//, rcomp2=false;
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			long x=Dedupe.baseToNumber[b];
			long x2=Dedupe.baseToComplementNumber[b];
			kmer=((kmer<<2)|x)&mask;
			rkmer=(rkmer>>>2)|(x2<<shift2);
			if(b=='N'){
				len=0;
				mod1=mod2=-1;
				kmax1=-1;//kmax2=-1;
			}else{len++;}
			if(len>=k){
				final long kmax0=Tools.max(kmer, rkmer);
				final long mod0=kmax0%divisor;
				final boolean rcomp0=(kmax0!=kmer);
				if(len>k+1 && mod1>topMod && mod1>mod2 && mod1>mod0){//Local maximum
					if(minCount<2 || table.getCount(kmer, rkmer)>=minCount){
						topMod=mod1;
						kmers[0]=kmax1;
						kmers[1]=i-1;
						rcomp=(rcomp1);
					}
				}
				mod2=mod1;
				mod1=mod0;
//				kmax2=kmax1;
				kmax1=kmax0;
//				rcomp2=rcomp1;
				rcomp1=rcomp0;
			}
		}
		
		if(topMod<0){//There was no local maximum
			if(minCount>1){return fillLocalMax(r, kmers, null, 0);}
			else{return fillMax(r, kmers, table, minCount);}
		}
		
		rcomp&=rcompReads;
//		r.id+=" "+kmers[1]+","+rcomp+","+(bases.length-kmers[1]+k-2);
		if(rcomp){
			r.reverseComplement();
			r.setSwapped(true);
			kmers[1]=bases.length-kmers[1]+k-2;
		}
		if(addName){r.id+=" "+kmers[1]+(rcomp ? ",t" : ",f")+","+kmers[0];}
		assert(kmers[0]>=0 && kmers[1]>=0) : Arrays.toString(kmers);
		return kmers[0];
	}
	
	private int compare(long[] alist, long[] blist){
		for(int i=0; i<alist.length; i++){
			final long a=alist[i], b=blist[i];
			if(a!=b){
				return a>b ? 1 : -1;
			}
		}
		return 0;
	}
	
	private class HashThread extends Thread{
		
		HashThread(int id_, int threads_, ArrayList<Read> list_){
			id=id_;
			threads=threads_;
			list=list_;
		}
		
		@Override
		public void run(){
			for(int i=id; i<list.size(); i+=threads){
				hash(list.get(i));
			}
		}
		
		final int id;
		final int threads;
		final ArrayList<Read> list;
		
	}
	
	public final int k;

	final int shift;
	final int shift2;
	final long mask;
	
	public final long divisor;
	public boolean addName=true;
	public boolean rcompReads=true;
	
	public static final boolean LOCAL_MAX=false; //Should improve compression, but decreases compression...?

}
