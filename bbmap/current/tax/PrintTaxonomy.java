package tax;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import dna.Parser;
import dna.Timer;
import fileIO.ReadWrite;
import fileIO.FileFormat;
import fileIO.TextFile;
import fileIO.TextStreamWriter;
import stream.ConcurrentReadInputStream;
import stream.Read;
import structures.ListNum;
import structures.LongList;
import align2.ReadStats;
import align2.Shared;
import align2.Tools;

/**
 * Filters sequences according to their taxonomy,
 * as determined by the sequence name.  Sequences should
 * be labeled with a gi number or NCBI taxID.
 * 
 * @author Brian Bushnell
 * @date November 23, 2015
 *
 */
public class PrintTaxonomy {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		Timer t=new Timer();
		PrintTaxonomy as=new PrintTaxonomy(args);
		as.process(t);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public PrintTaxonomy(String[] args){
		
		//Process any config files
		args=Parser.parseConfig(args);
		
		//Detect whether the uses needs help
		if(Parser.parseHelp(args, true)){
			printOptions();
			System.exit(0);
		}
		
		//Print the program name and arguments
		outstream.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		
		//Set some shared static variables regarding PIGZ
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=Shared.threads();
		
		//Create a parser object
		Parser parser=new Parser();
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b==null || b.equalsIgnoreCase("null")){b=null;}
			while(a.startsWith("-")){a=a.substring(1);} //Strip leading hyphens
			
			
			if(a.equals("out")){
				out1=b;
			}else if(a.equals("counts")){
				countFile=b;
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
			}else if(a.equals("table") || a.equals("gi") || a.equals("gitable")){
				tableFile=b;
				if("auto".equalsIgnoreCase(b)){tableFile=TaxTree.DefaultTableFile;}
			}else if(a.equals("tree") || a.equals("taxtree")){
				treeFile=b;
				if("auto".equalsIgnoreCase(b)){treeFile=TaxTree.DefaultTreeFile;}
			}else if(a.equals("level") || a.equals("taxlevel")){
				if(b==null){taxLevel=-1;}
				else if(Character.isDigit(b.charAt(0))){
					taxLevel=Integer.parseInt(b);
				}else{
					taxLevel=TaxTree.stringToLevel(b.toLowerCase());
				}
			}else if(a.equals("printname")){
				printName=Tools.parseBoolean(b);
			}else if(b!=null && (a.equals("name") || a.equals("names") || a.equals("id") || a.equals("ids"))){
				for(String s : b.split(",")){
					names.add(s);
				}
			}else{
				names.add(arg);
			}
		}
		
		{//Process parser fields
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;

			in1=parser.in1;
			maxReads=parser.maxReads;
		}
		
		//Ensure output files can be written
		if(!Tools.testOutputFiles(overwrite, append, false, out1)){
			outstream.println((out1==null)+", "+out1);
			throw new RuntimeException("\n\noverwrite="+overwrite+"; Can't write to output file "+out1+"\n");
		}
		
		//Create output FileFormat objects
		ffout1=FileFormat.testOutput(out1, FileFormat.TEXT, null, true, overwrite, append, ordered);
		
		ffcount=FileFormat.testOutput(countFile, FileFormat.TEXT, null, true, overwrite, append, ordered);
		
		//Create input FileFormat objects
		ffin1=FileFormat.testInput(in1, FileFormat.TEXT, null, true, false);
		
		if(tableFile!=null){
			outstream.println("Loading gi table.");
			GiToNcbi.initialize(tableFile);
		}
		if(treeFile!=null){
			outstream.println("Loading tree.");
			tree=ReadWrite.read(TaxTree.class, treeFile, true);
			if(tree.nameMap==null){
				outstream.println("Hashing names.");
				tree.hashNames();
			}
			assert(tree.nameMap!=null);
		}else{
			tree=null;
			throw new RuntimeException("No tree specified.");
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	/** Create read streams and process all data */
	void process(Timer t){
		
		TextStreamWriter tsw=null;
		if(ffout1!=null){
			tsw=new TextStreamWriter(ffout1);
			tsw.start();
		}
		
		if(ffin1!=null){
			if(ffin1.fasta() || ffin1.fastq() || ffin1.samOrBam() || ffin1.scarf()){
				processReads(tsw);
			}else{
				processFile(new TextFile(ffin1), tsw);
			}
		}else{
			processNames(tsw);
		}
		
		if(tsw!=null){errorState|=tsw.poisonAndWait();}
		
		if(ffcount!=null){
			TextStreamWriter tswc=new TextStreamWriter(ffcount);
			tswc.start();
			for(TaxNode tn : tree.nodes){
				if(tn!=null && tn.countRaw>0){
					tswc.println(tn.countRaw+"\t"+tn.name);
				}
			}
			errorState|=tswc.poisonAndWait();
		}
		
		t.stop();
		
		//Throw an exception of there was an error in a thread
		if(errorState){
			throw new RuntimeException(getClass().getName()+" terminated in an error state; the output may be corrupt.");
		}
	}
	
	/** Iterate through the names */
	void processNames(final TextStreamWriter tsw){
		for(String name : names){
			if(taxLevel>-1){
				printTaxLevel(name, tsw);
			}else{
				printTaxonomy(name, tsw);
			}
		}
	}
	
	/** Iterate through the names */
	void processFile(final TextFile tf, final TextStreamWriter tsw){
		for(String name=tf.nextLine(); name!=null; name=tf.nextLine()){
			if(taxLevel>-1){
				printTaxLevel(name, tsw);
			}else{
				printTaxonomy(name, tsw);
			}
		}
	}
	
	/** Iterate through the names */
	void processReads(final TextStreamWriter tsw){
		final ConcurrentReadInputStream cris;
		{
			cris=ConcurrentReadInputStream.getReadInputStream(maxReads, true, ffin1, null);
			if(verbose){System.err.println("Started cris");}
			cris.start();
		}
		
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> reads=(ln!=null ? ln.list : null);
		
		while(reads!=null && reads.size()>0){

			for(Read r1 : reads){
				if(taxLevel>-1){
					printTaxLevel(r1.id, tsw);
				}else{
					printTaxonomy(r1.id, tsw);
				}
			}
			cris.returnList(ln.id, ln.list.isEmpty());
			ln=cris.nextList();
			reads=(ln!=null ? ln.list : null);
		}
		cris.returnList(ln.id, ln.list.isEmpty());
		ReadWrite.closeStreams(cris);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	void printTaxonomy(String name, final TextStreamWriter tsw){
		while(name.startsWith(">") || name.startsWith("@")){name=name.substring(1);}
		tsw.print("\n");
		if(printName){tsw.print(name+":\n");}
		TaxNode tn=tree.getNode(name);
		if(tn!=null){
			printTaxonomy(tn, tsw);
			return;
		}else{
			List<TaxNode> list=tree.getNodesByName(name);
			if(list!=null){
				String nl="";
				for(TaxNode tn2 : list){
					tsw.print(nl);
					printTaxonomy(tn2, tsw);
					nl="\n";
				}
				return;
			}
		}
		tsw.println("Could not find node" + (printName ? "." : " for '"+name+"'"));
		return;
	}
	
	void printTaxLevel(String name, final TextStreamWriter tsw){
		while(name.startsWith(">") || name.startsWith("@")){name=name.substring(1);}
		tsw.print("\n");
		if(printName){tsw.print(name+":\n");}
		TaxNode tn=tree.getNode(name);
		if(tn!=null){
			printTaxLevel(tn, tsw);
			return;
		}else{
			List<TaxNode> list=tree.getNodesByName(name);
			if(list!=null){
				for(TaxNode tn2 : list){
					printTaxLevel(tn2, tsw);
				}
				return;
			}
		}
		tsw.println("Could not find node" + (printName ? "." : " for '"+name+"'"));
		return;
	}
	
//	void printTaxCounts(String name, final TextStreamWriter tsw){
//		TaxNode tn=null;
//		tn=tree.getNode(name);
//		if(tn==null){tn=tree.getNodeByName(name);}
//		if(tn==null){tn=unknown;}
//		while(tn!=null && tn.id!=tn.pid && tn.level<taxLevel){tn=tree.getNode(tn.pid);}
//		if(tsw!=null)tsw.println(tn.name);
//		tn.incrementRaw(1);
//	}
	
	void printTaxonomy(TaxNode tn, final TextStreamWriter tsw){
		assert(tn!=null);
//		tsw.print("\n");
		do{
			if(tn.level<=taxLevel){tn.incrementRaw(1);}
			tsw.println(tn.levelString()+"\t"+tn.id+"\t"+tn.name);
			tn=tree.getNode(tn.pid);
		}while(tn!=null && tn.id!=tn.pid);
	}
	
	public static void printTaxonomy(TaxNode tn, final TextStreamWriter tsw, final TaxTree tree, final int maxLevel){
		assert(tn!=null);
//		tsw.print("\n");
		do{
			tsw.println(tn.levelString()+"\t"+tn.id+"\t"+tn.name);
			tn=tree.getNode(tn.pid);
		}while(tn!=null && tn.id!=tn.pid && tn.level<=maxLevel);
	}
	
	void printTaxLevel(TaxNode tn, final TextStreamWriter tsw){
		if(tn==null){tn=unknown;}
		while(tn!=null && tn.id!=tn.pid && tn.level<taxLevel){tn=tree.getNode(tn.pid);}
		if(tsw!=null)tsw.println(tn.name);
		tn.incrementRaw(1);
	}
	
//	void printTaxCounts(TaxNode tn, final TextStreamWriter tsw){
//		if(tn==null){tn=unknown;}
//		while(tn!=null && tn.id!=tn.pid && tn.level<taxLevel){tn=tree.getNode(tn.pid);}
//		if(tsw!=null)tsw.println(tn.name);
//		tn.incrementRaw(1);
//	}
	
	/** This is called if the program runs with no parameters */
	private void printOptions(){
		throw new RuntimeException("TODO");
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Optional input file path */
	private String in1=null;

	/** Primary output file path */
	private String out1="stdout.txt";
	
	private String countFile=null;

	private String tableFile=null;;
	private String treeFile=TaxTree.DefaultTreeFile;
	
	private final TaxTree tree;
	
	/** Level to print */
	private int taxLevel=-1;//TaxTree.stringToLevel("phylum");
	
	private ArrayList<String> names=new ArrayList<String>();
	
	private long maxReads=-1;
	
	boolean printName=true;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Optional input file */
	private final FileFormat ffin1;
	
	/** Primary output file */
	private final FileFormat ffout1;
	
	private final FileFormat ffcount;
	
	private final TaxNode unknown=new TaxNode(-99, -99, taxLevel, "UNKNOWN");
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print status messages to this output stream */
	private PrintStream outstream=System.err;
	/** Print verbose messages */
	public static boolean verbose=false;
	/** True if an error was encountered */
	public boolean errorState=false;
	/** Overwrite existing output files */
	private boolean overwrite=false;
	/** Append to existing output files */
	private boolean append=false;
	/** This flag has no effect on singlethreaded programs */
	private final boolean ordered=false;
	
}
