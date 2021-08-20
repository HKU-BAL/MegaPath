package clump;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;

import kmer.KmerTableSet;

import stream.ConcurrentReadInputStream;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.ConcurrentReadOutputStream;
import stream.Read;
import structures.ListNum;
import dna.Parser;
import dna.Timer;
import fileIO.ByteFile;
import fileIO.ReadWrite;
import fileIO.FileFormat;
import align2.ReadStats;
import align2.Shared;
import align2.Tools;

/**
 * @author Brian Bushnell
 * @date June 20, 2014
 *
 */
public class KmerSplit {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		final boolean pigz=ReadWrite.USE_PIGZ, unpigz=ReadWrite.USE_UNPIGZ;
		Timer t=new Timer();
		KmerSplit ks=new KmerSplit(args);
		ks.process(t);
		ReadWrite.USE_PIGZ=pigz;
		ReadWrite.USE_UNPIGZ=unpigz;
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public KmerSplit(String[] args){
		
		args=Parser.parseConfig(args);
		if(Parser.parseHelp(args, true)){
			printOptions();
			System.exit(0);
		}
		
		outstream.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		
		boolean setInterleaved=false; //Whether it was explicitly set.

		
		
		Shared.READ_BUFFER_LENGTH=Tools.min(200, Shared.READ_BUFFER_LENGTH);
		ReadWrite.USE_PIGZ=false;
		ReadWrite.USE_UNPIGZ=true;
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
			}else if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
			}else if(a.equals("parse_flag_goes_here")){
				//Set a variable here
			}else if(a.equals("k")){
				k=Integer.parseInt(b);
				assert(k>0 && k<32);
			}else if(a.equals("mincount") || a.equals("mincr")){
				minCount=Integer.parseInt(b);
			}else if(a.equals("groups") || a.equals("g") || a.equals("sets")){
				groups=Integer.parseInt(b);
			}else if(a.equals("divisor") || a.equals("div") || a.equals("mindivisor")){
				minDivisor=Tools.parseKMG(b);
			}else if(a.equals("rename") || a.equals("addname")){
				addName=Tools.parseBoolean(b);
			}else if(a.equals("rcomp") || a.equals("reversecomplement")){
				//ignore rcomp=Tools.parseBoolean(b);
			}else if(a.equals("condense")){
				//ignore condense=Tools.parseBoolean(b);
			}else if(a.equals("prefilter")){
				KmerReduce.prefilter=Tools.parseBoolean(b);
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

			setInterleaved=parser.setInterleaved;
			
			in1=parser.in1;
			in2=parser.in2;

			out1=parser.out1;
			
			extin=parser.extin;
			extout=parser.extout;
		}
		
		if(groups>2){ReadWrite.USE_PIGZ=false;}
		
		if(in1!=null && in2==null && in1.indexOf('#')>-1 && !new File(in1).exists()){
			in2=in1.replace("#", "2");
			in1=in1.replace("#", "1");
		}
		if(in2!=null){
			if(FASTQ.FORCE_INTERLEAVED){outstream.println("Reset INTERLEAVED to false because paired input files were specified.");}
			FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null){
			printOptions();
			throw new RuntimeException("Error - at least one input file is required.");
		}
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		if(!setInterleaved){
			assert(in1!=null) : "\nin1="+in1+"\nin2="+in2+"\nout1="+out1+"\n";
			if(in2!=null){ //If there are 2 input streams.
				FASTQ.FORCE_INTERLEAVED=FASTQ.TEST_INTERLEAVED=false;
				outstream.println("Set INTERLEAVED to "+FASTQ.FORCE_INTERLEAVED);
			}
		}

		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		
		if(out1!=null){
			assert(out1.contains("%"));
			outArray=new String[groups];
			for(int i=0; i<groups; i++){
				outArray[i]=out1.replaceFirst("%", ""+i);
			}
			if(!Tools.testOutputFiles(overwrite, append, false, outArray)){
				outstream.println((out1==null)+", "+out1);
				throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+"\n");
			}
			ffout=new FileFormat[groups];
			for(int i=0; i<groups; i++){
				ffout[i]=FileFormat.testOutput(outArray[i], FileFormat.FASTQ, extout, true, overwrite, append, false);
			}
		}else{
			outArray=null;
			throw new RuntimeException("out is a required parameter.");
		}

		ffin1=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
		ffin2=FileFormat.testInput(in2, FileFormat.FASTQ, extin, true, true);
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Count kmers */
	void preprocess(){
		if(minCount>1 && ClumpTools.table==null){
			table=KmerReduce.getValidKmersFromReads(in1, k, minCount);
			ClumpTools.table=table;
		}
	}

	/** Create read streams and process all data */
	void process(Timer t){
		
		preprocess();
		
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, ffin2, null, null);
			cris.start();
			if(verbose){outstream.println("Started cris");}
		}
		boolean paired=cris.paired();
		if(!ffin1.samOrBam()){outstream.println("Input is being processed as "+(paired ? "paired" : "unpaired"));}
		if(cris.paired() && (in1==null || !in1.contains(".sam"))){
			outstream.println("Writing interleaved.");
		}

		final ConcurrentReadOutputStream ros[]=new ConcurrentReadOutputStream[groups];
		for(int i=0; i<groups; i++){
			final int buff=8;

			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1)) : "Input file and output file have same name.";
			
			ros[i]=ConcurrentReadOutputStream.getStream(ffout[i], null, null, null, buff, null, false);
			ros[i].start();
		}
		
		readsProcessed=0;
		basesProcessed=0;
		
		//Process the read stream
		processInner(cris, ros);
		
		errorState|=ReadStats.writeAll();
		
		t.stop();
		
		double rpnano=readsProcessed/(double)(t.elapsed);
		double bpnano=basesProcessed/(double)(t.elapsed);

		String rpstring=(readsProcessed<100000 ? ""+readsProcessed : readsProcessed<100000000 ? (readsProcessed/1000)+"k" : (readsProcessed/1000000)+"m");
		String bpstring=(basesProcessed<100000 ? ""+basesProcessed : basesProcessed<100000000 ? (basesProcessed/1000)+"k" : (basesProcessed/1000000)+"m");

		while(rpstring.length()<8){rpstring=" "+rpstring;}
		while(bpstring.length()<8){bpstring=" "+bpstring;}
		
		outstream.println("Time:                         \t"+t);
		outstream.println("Reads Processed:    "+rpstring+" \t"+String.format("%.2fk reads/sec", rpnano*1000000));
		outstream.println("Bases Processed:    "+bpstring+" \t"+String.format("%.2fm bases/sec", bpnano*1000));
		
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/** Collect and sort the reads */
	void processInner(final ConcurrentReadInputStream cris, final ConcurrentReadOutputStream[] ros){
		if(verbose){outstream.println("Making comparator.");}
		KmerComparator kc=new KmerComparator(k, minDivisor);
		kc.addName=addName;
		kc.rcompReads=false;
		
		if(verbose){outstream.println("Splitting reads.");}
		splitReads(cris, ros, kc);
		
		if(verbose){outstream.println("Done!");}
	}
	
	public void splitReads(final ConcurrentReadInputStream cris, final ConcurrentReadOutputStream[] ros, final KmerComparator kc){
		if(verbose){outstream.println("Making hash threads.");}
		final int threads=Shared.threads();
		ArrayList<HashThread> alht=new ArrayList<HashThread>(threads);
		for(int i=0; i<threads; i++){alht.add(new HashThread(i, cris, ros, kc));}
		
		if(verbose){outstream.println("Starting threads.");}
		for(HashThread ht : alht){ht.start();}
		
		
		if(verbose){outstream.println("Waiting for threads.");}
		/* Wait for threads to die */
		for(HashThread ht : alht){
			
			/* Wait for a thread to die */
			while(ht.getState()!=Thread.State.TERMINATED){
				try {
					ht.join();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
			readsProcessed+=ht.readsProcessedT;
			basesProcessed+=ht.basesProcessedT;
		}
		
		if(verbose){outstream.println("Closing streams.");}
		errorState=ReadWrite.closeStreams(cris, ros)|errorState;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	/** This is called if the program runs with no parameters */
	private void printOptions(){
		throw new RuntimeException("TODO");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	private class HashThread extends Thread{

		HashThread(int id_, ConcurrentReadInputStream cris_, ConcurrentReadOutputStream[] ros_, KmerComparator kc_){
			id=id_;
			cris=cris_;
			ros=ros_;
			kc=kc_;
		}

		@Override
		public void run(){
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			ArrayList<Read>[] array=new ArrayList[groups];
			for(int i=0; i<groups; i++){
				array[i]=new ArrayList<Read>(buffer);
			}
			
			while(reads!=null && reads.size()>0){
				kc.hash(reads, table, minCount);
				for(Read r : reads){
					readsProcessedT++;
					basesProcessedT+=r.length();
					long kmer=((long[])r.obj)[0];
					long mod=kmer%kc.divisor;
					int mod2=(int)(mod%groups);
					array[mod2].add(r);
					if(array[mod2].size()>=buffer){
						ros[mod2].add(array[mod2], 0);
						array[mod2]=new ArrayList<Read>(buffer);
					}
				}
				cris.returnList(ln.id, ln.list.isEmpty());
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
			for(int i=0; i<groups; i++){
				if(!array[i].isEmpty()){
					ros[i].add(array[i], 0);
				}
			}
		}

		final int id;
		final ConcurrentReadInputStream cris;
		final ConcurrentReadOutputStream[] ros;
		final KmerComparator kc;
		static final int buffer=200;
		
		protected long readsProcessedT=0;
		protected long basesProcessedT=0;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private int k=31;
	private long minDivisor=80000000;
	private int groups=16;
	private int minCount=0;
	
	KmerTableSet table=null;
	
	/*--------------------------------------------------------------*/
	/*----------------          I/O Fields          ----------------*/
	/*--------------------------------------------------------------*/
	
	private String in1=null;
	private String in2=null;

	private String out1=null;
	private String[] outArray=null;
	
	private String extin=null;
	private String extout=null;
	
	/*--------------------------------------------------------------*/
	
	protected long readsProcessed=0;
	protected long basesProcessed=0;
	
	private long maxReads=-1;
	private boolean addName=false;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private final FileFormat ffin1;
	private final FileFormat ffin2;
	
	private final FileFormat[] ffout;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=false;
	private boolean append=false;
	
}
