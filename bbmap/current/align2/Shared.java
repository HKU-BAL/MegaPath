package align2;

import java.lang.management.ManagementFactory;
import java.util.List;

import dna.Data;

public class Shared {
	
	private static int THREADS=setThreads(-1);
	
	public static int READ_BUFFER_LENGTH=200;
	private static int READ_BUFFER_NUM_BUFFERS=setBuffers();
	public static long READ_BUFFER_MAX_DATA=400000;
	
	/** Temporary, for testing; should be made non-global */
	public static boolean AMINO_IN=false;
	
	//TODO:  For some reason, it seems as though GAPBUFFER must equal exactly 1/2 of GAPLEN.  Not good; 1/4 would be far better.
	
	public static final int GAPBUFFER=64; //TODO:  Seems to break less than 64, for some reason
	public static final int GAPBUFFER2=2*GAPBUFFER;
	public static final int GAPLEN=128; //TODO: May break when over 128
	public static final int MINGAP=GAPBUFFER2+GAPLEN;
	public static final int GAPCOST=Tools.max(1, GAPLEN/64);
	public static final byte GAPC='-';
	
	public static String BBMAP_VERSION_STRING="36.28";
	
	public static boolean TRIM_READ_COMMENTS=false;
	
	public static boolean USE_JNI=false;//Data.GENEPOOL;
	public static boolean USE_MPI=false;
	public static boolean MPI_KEEP_ALL=true;
	/** Use ConcurrentReadInputStreamMPI instead of D */
	public static boolean USE_CRISMPI=true;
	public static int MPI_RANK=0;
	public static int MPI_NUM_RANKS=1;
	
	public static int FASTA_WRAP=70;
	public static byte FAKE_QUAL=30;
	
	public static String BBMAP_CLASS=null;
	public static String[] COMMAND_LINE=null;
	public static List<String> JVM_ARGS(){
		return ManagementFactory.getRuntimeMXBean().getInputArguments();
	}
	
	public static long getAvailableMemory(){
		long usableMemory;
		{
			long memory=Runtime.getRuntime().maxMemory();
			double xmsRatio=Shared.xmsRatio();
			usableMemory=(long)Tools.max(((memory-96000000-(20*400000))*(xmsRatio>0.97 ? 0.82 : 0.72)), memory*0.45);
		}
		return usableMemory;
	}
	
	/** Directory in which to write temp files */
	public static String TMPDIR=(System.getenv("TMPDIR")==null ? null : (System.getenv("TMPDIR")+"/").replaceAll("//", "/"));
//	static{assert(false) : "TMPDIR="+TMPDIR;}
	
	/** Anomaly probably resolved as of v.20.1 
	 * This variable should be TRUE for normal users and FALSE for me. */
	public static boolean anomaly=!(System.getProperty("user.dir")+"").contains("/bushnell/") && !Data.WINDOWS;
	
	public static final char[] getTLCB(int len){
		char[] buffer=TLCB.get();
		if(buffer==null || buffer.length<len){
			buffer=new char[len];
			if(len<1000000){TLCB.set(buffer);}
		}
		return buffer;
	}
	private static final ThreadLocal<char[]> TLCB=new ThreadLocal<char[]>();
	
	public static int setThreads(String x){
		int y=Data.LOGICAL_PROCESSORS;
		if(x!=null && !x.equalsIgnoreCase("auto")){
			y=Integer.parseInt(x);
		}
		return setThreads(y);
	}
	
	public static int setThreads(int x){
		if(x>0){
			THREADS=x;
		}else{
			THREADS=Tools.max(1, Data.LOGICAL_PROCESSORS);
		}
		setBuffers();
		return THREADS;
	}
	
	public static int threads(){
		assert(THREADS>0);
		return THREADS;
	}
	
	public static int capBuffers(int num){
		return setBuffers(Tools.min(num, READ_BUFFER_NUM_BUFFERS));
	}
	
	public static int setBuffers(){
		return setBuffersFromThreads(THREADS);
	}
	
	public static int setBuffersFromThreads(int threads){
		return setBuffers(Tools.max(4, (threads*3)/2));
	}
	
	public static int setBuffers(int num){
		num=Tools.max(2, num);
		return READ_BUFFER_NUM_BUFFERS=num;
	}
	
	public static int numBuffers(){
		return READ_BUFFER_NUM_BUFFERS;
	}
	
	public static boolean LOW_MEMORY=false;
	
	/** Ratio of -Xms to -Xmx parameters */
	public static final double xmsRatio(){
		Runtime rt=Runtime.getRuntime();
		return rt.totalMemory()*1.0/rt.maxMemory();
	}
	
	/** Print statistics about current memory use and availability */
	public static final void printMemory(){
		try{
			if(GC_BEFORE_PRINT_MEMORY){
				System.gc();
				System.gc();
			}
			Runtime rt=Runtime.getRuntime();
			long mmemory=rt.maxMemory()/1000000;
			long tmemory=rt.totalMemory()/1000000;
			long fmemory=rt.freeMemory()/1000000;
			long umemory=tmemory-fmemory;
			System.err.println("Memory: "+"max="+mmemory+/*"m, total="+tmemory+*/"m, "+"free="+fmemory+"m, used="+umemory+"m");
		}catch(Throwable t){}
	}
	
	/** Do garbage collection prior to printing memory usage */
	private static final boolean GC_BEFORE_PRINT_MEMORY=false;
	
}
