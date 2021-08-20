package clump;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import kmer.KmerTableSet;

import stream.ConcurrentReadInputStream;
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
public class KmerSort {
	
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
		KmerSort ks=new KmerSort(args);
		ks.process(t);
		ReadWrite.USE_PIGZ=pigz;
		ReadWrite.USE_UNPIGZ=unpigz;
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public KmerSort(String[] args){
		
		args=Parser.parseConfig(args);
		if(Parser.parseHelp(args, true)){
			printOptions();
			System.exit(0);
		}
		
		outstream.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		
		boolean setInterleaved=false; //Whether it was explicitly set.
		
		Shared.READ_BUFFER_LENGTH=Tools.min(200, Shared.READ_BUFFER_LENGTH);
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
			}else if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
			}else if(a.equals("parse_flag_goes_here")){
				//Set a variable here
			}else if(a.equals("k")){
				k=Integer.parseInt(b);
				assert(k>0 && k<32);
			}else if(a.equals("mincount") || a.equals("mincr")){
				minCount=Integer.parseInt(b);
			}else if(a.equals("comparisons") || a.equals("c")){
				comparisons=Integer.parseInt(b);
			}else if(a.equals("divisor") || a.equals("div") || a.equals("mindivisor")){
				minDivisor=Tools.parseKMG(b);
			}else if(a.equals("rename") || a.equals("addname")){
				addName=Tools.parseBoolean(b);
//			}else if(a.equals("cache")){
//				KmerComparator.useCache=Tools.parseBoolean(b);//Obsolete
			}else if(a.equals("rcomp") || a.equals("reversecomplement")){
				rcomp=Tools.parseBoolean(b);
			}else if(a.equals("condense") || a.equals("consensus") || a.equals("concensus")){//Note the last one is intentionally misspelled
				condense=Tools.parseBoolean(b);
			}else if(a.equals("prefilter")){
				KmerReduce.prefilter=Tools.parseBoolean(b);
			}else if(a.equals("groups") || a.equals("g") || a.equals("sets")){
				groups=Integer.parseInt(b);
				splitInput=(groups>1);
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

			out1=parser.out1;
			
			extin=parser.extin;
			extout=parser.extout;
		}
		
		assert(FastaReadInputStream.settingsOK());
		
		if(in1==null){
			printOptions();
			throw new RuntimeException("Error - at least one input file is required.");
		}
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}

		if(out1!=null && out1.equalsIgnoreCase("null")){out1=null;}
		
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output files "+out1+"\n");
		}
		
		ffout1=FileFormat.testOutput(out1, FileFormat.FASTQ, extout, true, overwrite, append, false);
		
		if(groups>1 && in1.contains("%") && (splitInput || !new File(in1).exists())){
			ffin=new FileFormat[groups];
			for(int i=0; i<groups; i++){
				ffin[i]=FileFormat.testInput(in1.replaceFirst("%", ""+i), FileFormat.FASTQ, extin, true, true);
			}
		}else{
			assert(!in1.contains("%") && groups==1) : "The % symbol must only be present in the input filename if groups>1.";
			ffin=new FileFormat[1];
			ffin[0]=FileFormat.testInput(in1, FileFormat.FASTQ, extin, true, true);
			groups=1;
		}
//		if(groups>1){ReadWrite.USE_UNPIGZ=false;} //Not needed since they are not concurrent
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
		
		final ConcurrentReadInputStream[] cris=new ConcurrentReadInputStream[groups];
		for(int i=0; i<cris.length; i++){
			cris[i]=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin[i], null, null, null);
		}

		final ConcurrentReadOutputStream ros;
		if(out1!=null){
			final int buff=1;
			
			if(cris[0].paired()  && (in1==null || !in1.contains(".sam"))){
				outstream.println("Writing interleaved.");
			}			
			
			assert(!out1.equalsIgnoreCase(in1) && !out1.equalsIgnoreCase(in1)) : "Input file and output file have same name.";
			
			ros=ConcurrentReadOutputStream.getStream(ffout1, null, null, null, buff, null, false);
			ros.start();
		}else{ros=null;}
		
		readsProcessed=0;
		basesProcessed=0;
		
		//Process the read stream
		processInner(cris, ros);
		
		table=ClumpTools.table=null;
		
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
	void processInner(final ConcurrentReadInputStream[] crisArray, final ConcurrentReadOutputStream ros){
		if(verbose){outstream.println("Making comparator.");}
		KmerComparator kc=new KmerComparator(k, minDivisor);
		kc.addName=addName;
		kc.rcompReads=rcomp;
		
		int i=0;
		for(ConcurrentReadInputStream cris : crisArray){
			i++;
			if(verbose){outstream.println("Starting cris "+i+".");}
			cris.start();
			
			if(verbose){outstream.println("Fetching reads.");}
			ArrayList<Read> reads=fetchReads(cris, kc);

			if(verbose){outstream.println("Sorting.");}
			Collections.sort(reads, kc);

			if(condense){
				if(verbose){outstream.println("Condensing.");}
				reads=condenseReads(reads);
			}

			if(ros!=null){
				if(verbose){outstream.println("Writing.");}
				ros.add(reads, 0);
			}
		}
		
		if(ros!=null){
			if(verbose){outstream.println("Waiting for writing to complete.");}
			errorState=ReadWrite.closeStream(ros)|errorState;
		}
		
		if(verbose){outstream.println("Done!");}
	}
	
	public ArrayList<Read> fetchReads(final ConcurrentReadInputStream cris, final KmerComparator kc){
		if(verbose){outstream.println("Making hash threads.");}
		final int threads=Shared.threads();
		ArrayList<HashThread> alht=new ArrayList<HashThread>(threads);
		for(int i=0; i<threads; i++){alht.add(new HashThread(i, cris, kc));}
		
		if(verbose){outstream.println("Starting threads.");}
		for(HashThread ht : alht){ht.start();}
		
		
		if(verbose){outstream.println("Waiting for threads.");}
		long readsThisPass=0;
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
			readsThisPass+=ht.readsProcessedT;
			basesProcessed+=ht.basesProcessedT;
		}
		readsProcessed+=readsThisPass;
		
		if(verbose){outstream.println("Closing input stream.");}
		errorState=ReadWrite.closeStream(cris)|errorState;
		
		if(verbose){outstream.println("Combining thread output.");}
		assert(readsProcessed<=Integer.MAX_VALUE);
		ArrayList<Read> list=new ArrayList<Read>((int)readsThisPass);
		for(int i=0; i<threads; i++){
			HashThread ht=alht.set(i, null);
			list.addAll(ht.storage);
		}
		
		assert(list.size()==readsThisPass);
		return list;
	}
	
	public ArrayList<Read> condenseReads(ArrayList<Read> list){
		ClumpList cl=new ClumpList(list);
		list.clear();
		ArrayList<Read> out=cl.condense();
		cl.clear();
		return out;
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
		
		HashThread(int id_, ConcurrentReadInputStream cris_, KmerComparator kc_){
			id=id_;
			cris=cris_;
			kc=kc_;
			storage=new ArrayList<Read>();
		}
		
		@Override
		public void run(){
			
			ListNum<Read> ln=cris.nextList();
			ArrayList<Read> reads=(ln!=null ? ln.list : null);
			
			while(reads!=null && reads.size()>0){
				kc.hash(reads, table, minCount);
				for(Read r : reads){
					readsProcessedT++;
					basesProcessedT+=r.length();
				}
				storage.addAll(reads);
				cris.returnList(ln.id, ln.list.isEmpty());
				ln=cris.nextList();
				reads=(ln!=null ? ln.list : null);
			}
			if(ln!=null){
				cris.returnList(ln.id, ln.list==null || ln.list.isEmpty());
			}
			Collections.sort(storage, kc);//Optimization for TimSort
		}

		final int id;
		final ConcurrentReadInputStream cris;
		final KmerComparator kc;
		final ArrayList<Read> storage;
		
		protected long readsProcessedT=0;
		protected long basesProcessedT=0;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private int k=31;
	private int minCount=0;
	private int comparisons=3;
	private long minDivisor=80000000;
	
	private int groups=1;
	
	KmerTableSet table=null;
	
	/*--------------------------------------------------------------*/
	/*----------------          I/O Fields          ----------------*/
	/*--------------------------------------------------------------*/
	
	private String in1=null;

	private String out1=null;
	
	private String extin=null;
	private String extout=null;
	
	/*--------------------------------------------------------------*/
	
	protected long readsProcessed=0;
	protected long basesProcessed=0;
	
	private long maxReads=-1;
	private boolean addName=true;
	private boolean rcomp=true;
	private boolean condense=false;
	private boolean splitInput=false;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private final FileFormat ffin[];

	private final FileFormat ffout1;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	private PrintStream outstream=System.err;
	public static boolean verbose=false;
	public boolean errorState=false;
	private boolean overwrite=false;
	private boolean append=false;
	
}
