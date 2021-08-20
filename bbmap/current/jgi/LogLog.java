package jgi;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Random;
import java.util.concurrent.atomic.AtomicIntegerArray;

import kmer.Primes;

import dna.Parser;
import dna.Timer;

import fileIO.FileFormat;
import fileIO.ReadWrite;

import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadInputStream;
import stream.FastaReadInputStream;
import stream.Read;
import structures.ListNum;
import ukmer.Kmer;
import align2.ReadStats;
import align2.Shared;
import align2.Tools;

/**
 * @author Brian Bushnell
 * @date Sep 30, 2015
 *
 */
public class LogLog {
	
	public static void main(String[] args){
		LogLogWrapper llw=new LogLogWrapper(args);
		llw.process();
	}
	
	public final long cardinality(){
		long sum=0;
		for(int i=0; i<maxArray.length(); i++){
			sum+=maxArray.get(i);
		}
		double mean=sum/(double)buckets;
		long cardinality=(long)((((Math.pow(2, mean)-1)*buckets*SKIPMOD))/1.262);
		lastCardinality=cardinality;
		return cardinality;
	}
	
	public final long cardinalityH(){
		double sum=0;
		for(int i=0; i<maxArray.length(); i++){
			int x=Tools.max(1, maxArray.get(i));
			sum+=1.0/x;
		}
		double mean=buckets/sum;
		return (long)((Math.pow(2, mean)*buckets*SKIPMOD));
	}
	
	public LogLog(Parser p){
		this(p.loglogbuckets, p.loglogbits, p.loglogk, p.loglogseed);
	}
	
	public LogLog(int buckets_, int bits_, int k_, long seed){
//		hashes=hashes_;
		buckets=buckets_;
		bits=bits_;
		k=Kmer.getKbig(k_);
		maxArray=(atomic ? new AtomicIntegerArray(buckets) : null);
		maxArray2=(atomic ? null : new long[buckets]);
		steps=(63+bits)/bits;
		tables=new long[numTables][][];
		for(int i=0; i<numTables; i++){
			tables[i]=makeCodes(steps, bits, (seed<0 ? -1 : seed+i));
		}
		
//		assert(false) : "steps="+steps+", "+tables.length+", "+tables[0].length+", "+tables[0][0].length;
	}
	
	public long hash(final long value0, final long[][] table){
		long value=value0, code=value0;
		long mask=~((-1L)<<bits);
		
		for(int i=0; i<steps; i++){
			int x=(int)(value&mask);
			value>>=bits;
			code=Long.rotateLeft(code^table[i][x], 3);
		}
		return Long.rotateLeft(code, (int)(value0&31));
	}
	
	public void add(long number){
		hash(number);
	}
	
	public void hash(Read r){
		if(r!=null && r.length()>=k){hash(r.bases);}
		if(r.mateLength()>=k){hash(r.mate.bases);}
	}
	
	public void hash(byte[] bases){
		if(k<32){hashSmall(bases);}
		else{hashBig(bases);}
	}
	
	public void hashSmall(byte[] bases){
		final int shift=2*k;
		final int shift2=shift-2;
		final long mask=~((-1L)<<shift);
		int len=0;
		
		long kmer=0, rkmer=0;
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			long x=Dedupe.baseToNumber[b];
			long x2=Dedupe.baseToComplementNumber[b];
			kmer=((kmer<<2)|x)&mask;
			rkmer=(rkmer>>>2)|(x2<<shift2);
			if(b=='N'){len=0;}else{len++;}
			if(len>=k){
				add(Tools.max(kmer, rkmer));
			}
		}
	}
	
	public void hashBig(byte[] bases){
		
		Kmer kmer=getLocalKmer();
		int len=0;
		
		for(int i=0; i<bases.length; i++){
			byte b=bases[i];
			long x=Dedupe.baseToNumber[b];
			kmer.addRightNumeric(x);
			if(b=='N'){len=0;}else{len++;}
			if(len>=k){
				add(kmer.xor());
			}
		}
	}
	
	
	public void hash(final long number){
		if(number%SKIPMOD!=0){return;}
		long key=number;
		
		int i=(int)(number%5);
		key=Long.rotateRight(key, 1);
		key=hash(key, tables[i%numTables]);
		int leading=Long.numberOfLeadingZeros(key);
//		counts[leading]++;
		
		if(leading<3){return;}
		final int bucket=(int)((number&Integer.MAX_VALUE)%buckets);
		
//		if(maxArray!=null){
			int x=maxArray.get(bucket);
			while(leading>x){
				boolean b=maxArray.compareAndSet(bucket, x, leading);
				if(b){x=leading;}
				else{x=maxArray.get(bucket);}
			}
//		}else{
//			maxArray2[bucket]=Tools.max(leading, maxArray2[bucket]);
//		}
	}
	
	private static long[][] makeCodes(int length, int bits, long seed){
		Random randy;
		if(seed>=0){randy=new Random(seed);}
		else{randy=new Random();}
		int modes=1<<bits;
		long[][] r=new long[length][modes];
		for(int i=0; i<length; i++){
			for(int j=0; j<modes; j++){
				long x=randy.nextLong();
				while(Long.bitCount(x)>33){
					x&=(~(1L<<randy.nextInt(64)));
				}
				while(Long.bitCount(x)<31){
					x|=(1L<<randy.nextInt(64));
				}
				r[i][j]=x;
				
			}
		}
		return r;
	}
	
	public final int k;
	public final int numTables=4;
	public final int bits;
//	public final int hashes;
	public final int steps;
	private final long[][][] tables;
	public final AtomicIntegerArray maxArray;
	public final long[] maxArray2;
//	public final long[] counts=new long[64];
	public int buckets;
	private final ThreadLocal<Kmer> localKmer=new ThreadLocal<Kmer>();
	
	protected Kmer getLocalKmer(){
		Kmer kmer=localKmer.get();
		if(kmer==null){
			localKmer.set(new Kmer(k));
			kmer=localKmer.get();
		}
		kmer.clearFast();
		return kmer;
	}
	
	private static class LogLogWrapper{
		
		public LogLogWrapper(String[] args){

			Shared.READ_BUFFER_LENGTH=Tools.min(200, Shared.READ_BUFFER_LENGTH);
			Shared.capBuffers(4);
			ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
			ReadWrite.MAX_ZIP_THREADS=Shared.threads();
			

			Parser parser=new Parser();
			for(int i=0; i<args.length; i++){
				String arg=args[i];
				String[] split=arg.split("=");
				String a=split[0].toLowerCase();
				String b=split.length>1 ? split[1] : null;
				if(b==null || b.equalsIgnoreCase("null")){b=null;}
				while(a.startsWith("-")){a=a.substring(1);} //In case people use hyphens
				
				if(parser.parse(arg, a, b)){
					//do nothing
				}else if(a.equals("buckets") || a.equals("loglogbuckets")){
					long x=Tools.parseKMG(b);
					buckets=(int)Primes.primeAtLeast(Tools.min(1000000, x));
				}else if(a.equals("bits") || a.equals("loglogbits")){
					bits=Integer.parseInt(b);
				}else if(a.equals("k") || a.equals("loglogk")){
					k=Integer.parseInt(b);
				}else if(a.equals("seed") || a.equals("loglogseed")){
					seed=Long.parseLong(b);
				}else if(a.equals("verbose")){
					verbose=Tools.parseBoolean(b);
				}else if(a.equals("parse_flag_goes_here")){
					//Set a variable here
				}else if(in1==null && i==0 && !arg.contains("=") && (arg.toLowerCase().startsWith("stdin") || new File(arg).exists())){
					parser.in1=b;
				}else{
					outstream.println("Unknown parameter "+args[i]);
					assert(false) : "Unknown parameter "+args[i];
					//				throw new RuntimeException("Unknown parameter "+args[i]);
				}
			}

			{//Process parser fields
				Parser.processQuality();

				maxReads=parser.maxReads;

				overwrite=ReadStats.overwrite=parser.overwrite;
				append=ReadStats.append=parser.append;

				in1=(parser.in1==null ? null : parser.in1.split(","));
				in2=(parser.in2==null ? null : parser.in2.split(","));
				out=parser.out1;
			}
			
			assert(in1!=null && in1.length>0) : "No primary input file specified.";
			{
				ffin1=new FileFormat[in1.length];
				ffin2=new FileFormat[in1.length];
				
				for(int i=0; i<in1.length; i++){
					String a=in1[i];
					String b=(in2!=null && in2.length>i ? in2[i] : null);
					assert(a!=null) : "Null input filename.";
					if(b==null && a.indexOf('#')>-1 && !new File(a).exists()){
						b=a.replace("#", "2");
						a=a.replace("#", "1");
					}

					ffin1[i]=FileFormat.testInput(a, FileFormat.FASTQ, null, true, true);
					ffin2[i]=FileFormat.testInput(b, FileFormat.FASTQ, null, true, true);	
				}
			}

			assert(FastaReadInputStream.settingsOK());
		}
		
		
		void process(){
			Timer t=new Timer();
			LogLog log=new LogLog(buckets, bits, k, seed);
			
			
			for(int ffnum=0; ffnum<ffin1.length; ffnum++){
				ConcurrentReadInputStream cris=ConcurrentGenericReadInputStream.getReadInputStream(maxReads, false, ffin1[ffnum], ffin2[ffnum]);
				cris.start();

				LogLogThread[] threads=new LogLogThread[Shared.threads()];
				for(int i=0; i<threads.length; i++){
					threads[i]=new LogLogThread(log, cris);
				}
				for(LogLogThread llt : threads){
					llt.start();
				}
				for(LogLogThread llt : threads){
					while(llt.getState()!=Thread.State.TERMINATED){
						try {
							llt.join();
						} catch (InterruptedException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
					}
				}

				errorState|=ReadWrite.closeStreams(cris);
			}
			
			int[] copy=new int[log.maxArray.length()];
			for(int i=0; i<log.maxArray.length(); i++){
//				System.err.println(log.maxArray.get(i));
				copy[i]=log.maxArray.get(i);
			}
			
			t.stop();
			
			
			long cardinality=log.cardinality();
			
			if(out!=null){
				ReadWrite.writeString(cardinality+"\n", out);
			}
			
//			Arrays.sort(copy);
//			System.err.println("Median:        "+copy[Tools.median(copy)]);
			
//			System.err.println("Mean:          "+Tools.mean(copy));
//			System.err.println("Harmonic Mean: "+Tools.harmonicMean(copy));
			System.err.println("Cardinality:   "+log.cardinality());
//			System.err.println("CardinalityH:  "+log.cardinalityH());
			
//			for(long i : log.counts){System.err.println(i);}
			
			System.err.println("Time: \t"+t);
		}
		
		private class LogLogThread extends Thread{
			
			LogLogThread(LogLog log_, ConcurrentReadInputStream cris_){
				log=log_;
				cris=cris_;
			}
			
			@Override
			public void run(){
				ListNum<Read> ln=cris.nextList();
				ArrayList<Read> reads=(ln!=null ? ln.list : null);
				while(reads!=null && reads.size()>0){
					
					for(Read r : reads){
						log.hash(r);
					}
					
					cris.returnList(ln.id, ln.list.isEmpty());
					ln=cris.nextList();
					reads=(ln!=null ? ln.list : null);
				}
				cris.returnList(ln.id, ln.list.isEmpty());
			}
			
			private final LogLog log;
			private final ConcurrentReadInputStream cris;
			
		}
		
		/*--------------------------------------------------------------*/
		/*----------------            Fields            ----------------*/
		/*--------------------------------------------------------------*/
		
		private int buckets=1999;
		private int bits=8;
		private int k=31;
		private long seed=-1;
		
		
		private String[] in1=null;
		private String[] in2=null;
		private String out=null;
		
		/*--------------------------------------------------------------*/
		
		protected long readsProcessed=0;
		protected long basesProcessed=0;
		
		private long maxReads=-1;
		
		boolean overwrite=false;
		boolean append=false;
		boolean errorState=false;
		
		/*--------------------------------------------------------------*/
		/*----------------         Final Fields         ----------------*/
		/*--------------------------------------------------------------*/
		
		private final FileFormat[] ffin1;
		private final FileFormat[] ffin2;
		
		/*--------------------------------------------------------------*/
		/*----------------        Common Fields         ----------------*/
		/*--------------------------------------------------------------*/
	}
	
	private static PrintStream outstream=System.err;
	public static boolean verbose=false;
	public static boolean atomic=true;
	private static final long SKIPMOD=3;
	public static long lastCardinality=-1;
	
}
