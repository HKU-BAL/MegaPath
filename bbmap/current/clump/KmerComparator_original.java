package clump;

import java.util.Arrays;
import java.util.Comparator;

import jgi.Dedupe;
import align2.Tools;

import kmer.Primes;

import stream.Read;

/**
 * @author Brian Bushnell
 * @date Nov 4, 2015
 *
 */
public class KmerComparator_original implements Comparator<Read> {
	
	public KmerComparator_original(int k_, int comparisons_, long minDivisor_){
		k=k_;
		comparisons=comparisons_;
		assert(k>0 && k<32);
		assert(comparisons>0 && comparisons<1000);
		
		shift=2*k;
		shift2=shift-2;
		mask=~((-1L)<<shift);
		
		divisors=new long[comparisons];
		divisors[0]=Primes.primeAtLeast(minDivisor_);
		for(int i=1; i<comparisons; i++){
			divisors[i]=Primes.primeAtLeast(divisors[i-1]+1);
		}
	}
	
	/* (non-Javadoc)
	 * @see java.util.Comparator#compare(java.lang.Object, java.lang.Object)
	 */
	@Override
	public int compare(Read a, Read b) {
		final long[] alist, blist;
		if(useCache){
			if(a.obj==null){
				alist=new long[comparisons];
				a.obj=alist;
				fill(a, alist);
			}else{alist=(long[])a.obj;}

			if(b.obj==null){
				blist=new long[comparisons];
				b.obj=alist;
				fill(b, blist);
			}else{blist=(long[])b.obj;}
		}else{
			long[][] matrix=local1.get();
			if(matrix==null){
				matrix=new long[2][comparisons];
				local1.set(matrix);
			}
			alist=matrix[0];
			blist=matrix[1];
			fill(a, alist);
			fill(b, blist);
		}
		
		return compare(alist, blist);
	}
	
	public void fill(Read r, long[] kmers){
		final byte[] bases=r.bases;
		long kmer=0;
		long rkmer=0;
		int len=0;
		
		if(bases==null || bases.length<k){return;}
		
		long[] mods=local2.get();
		if(mods==null){
			mods=new long[comparisons];
			local2.set(mods);
		}
		Arrays.fill(mods, -1);
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			long x=Dedupe.baseToNumber[b];
			long x2=Dedupe.baseToComplementNumber[b];
			kmer=((kmer<<2)|x)&mask;
			rkmer=(rkmer>>>2)|(x2<<shift2);
			if(b=='N'){len=0;}else{len++;}
			if(len>=k){
				final long kmax=Tools.max(kmer, rkmer);
				for(int j=0; j<comparisons; j++){
					final long div=divisors[j];
					final long mod=kmax%div;
					if(mod>mods[j]){
						mods[j]=mod;
						kmers[j]=kmax;
					}
				}
			}
		}
	}
	
	private int compare(long[] alist, long[] blist){
		for(int i=0; i<comparisons; i++){
			final long a=alist[i], b=blist[i];
			if(a!=b){
				return a>b ? 1 : -1;
			}
		}
		return 0;
	}
	
	public final int k;

	final int shift;
	final int shift2;
	final long mask;
	
	public final int comparisons;
	public final long[] divisors;
	public static boolean useCache=true;

	private ThreadLocal<long[][]> local1=new ThreadLocal<long[][]>();
	private ThreadLocal<long[]> local2=new ThreadLocal<long[]>();

}
