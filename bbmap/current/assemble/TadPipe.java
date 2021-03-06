package assemble;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import align2.Shared;
import align2.Tools;
import dna.Parser;
import fileIO.ReadWrite;
import jgi.BBDukF;
import jgi.BBMerge;

/**
 * Merges and error-corrects reads, then assembles.
 * @author Brian Bushnell
 * @date August 16, 2016
 *
 */
public class TadPipe {
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		TadPipe as=new TadPipe(args);
		as.process();
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public TadPipe(String[] args){
		
		//Process any config files
		args=Parser.parseConfig(args);
		
		//Detect whether the uses needs help
		if(Parser.parseHelp(args, true)){
			System.exit(1);
		}
		
		//Print the program name and arguments
		outstream.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		
		//Set some shared static variables regarding PIGZ
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=1+Shared.threads()/2;
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b==null || b.equalsIgnoreCase("null")){b=null;}
			while(a.startsWith("-")){a=a.substring(1);} //Strip leading hyphens
			
			if(Parser.isJavaFlag(arg)){
				//do nothing
			}else if(a.equals("in") || a.equals("in1")){
				in1=b;
			}else if(a.equals("in2")){
				in2=b;
			}else if(a.equals("out")){
				out=b;
			}else if(a.equals("phix") || a.equals("contam")){
				phix=b;
			}else if(a.equals("delete")){
				deleteTemp=Tools.parseBoolean(b);
			}else if(a.equals("gz")){
				gz=Tools.parseBoolean(b);
			}else if(a.equals("tmpdir") || a.equals("tempdir") || a.equals("temp")){
				tempdir=(b+"/").replaceAll("\\\\", "/").replaceAll("//", "/");
			}else if(a.startsWith("merge_")){
				mergeArgs.add(arg.substring(arg.indexOf('_')+1));
			}else if(a.startsWith("ecc_") || a.startsWith("correct_")){
				eccArgs.add(arg.substring(arg.indexOf('_')+1));
			}else if(a.startsWith("extend_") || a.startsWith("extend1_")){
				extendArgs.add(arg.substring(arg.indexOf('_')+1));
			}else if(a.startsWith("extend2_")){
				extend2Args.add(arg.substring(arg.indexOf('_')+1));
			}else if(a.startsWith("filter_")){
				filterArgs.add(arg.substring(arg.indexOf('_')+1));
			}else if(a.startsWith("trim_")){
				trimArgs.add(arg.substring(arg.indexOf('_')+1));
			}else if(a.startsWith("assemble_")){
				assembleArgs.add(arg.substring(arg.indexOf('_')+1));
			}else{
				assembleArgs.add(arg);
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/
	
	public void process(){
		
		File tmpfile=(tempdir==null ? null : new File(tempdir));
		if(!tmpfile.exists()){tmpfile.mkdirs();}
		
		String outFilter=null;
		
		String outAdapter=null;
		
		String outTrimmed=null;

		String outMerged=null;
		String outUnmerged=null;

		String outMergedEcc=null;
		String outUnmergedEcc=null;

		String outMergedExtended=null;
		String outUnmergedExtended=null;

		String outMergedExtended2=null;
		String outUnmergedExtended2=null;
		
		String outMultiK=null;
		
		String ext=(gz ? ".fq.gz" : ".fq");
		
		try {
			outFilter=File.createTempFile("filter_", ext, tmpfile).toString();
			outAdapter=File.createTempFile("adapters_", ".fa", tmpfile).toString();
			outTrimmed=File.createTempFile("trimmed_", ext, tmpfile).toString();
			outMerged=File.createTempFile("merged_", ext, tmpfile).toString();
			outUnmerged=File.createTempFile("unmerged_", ext, tmpfile).toString();
			outMergedEcc=File.createTempFile("m_ecc_", ext, tmpfile).toString();
			outUnmergedEcc=File.createTempFile("u_ecc_", ext, tmpfile).toString();
			outMergedExtended=File.createTempFile("m_extended_", ext, tmpfile).toString();
			outUnmergedExtended=File.createTempFile("u_extended_", ext, tmpfile).toString();
			if(extend2){
				outMergedExtended2=File.createTempFile("m_extended2_", ext, tmpfile).toString();
				outUnmergedExtended2=File.createTempFile("u_extended2_", ext, tmpfile).toString();
			}
			outMultiK=File.createTempFile("multik_%_", ".fa", tmpfile).toString();
			delete(outMultiK);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		{
			filterArgs.add("in="+in1);
			if(in2!=null){filterArgs.add("in2="+in2);}

			filterArgs.add("out="+outFilter);
			filterArgs.add("ftm="+5);
			filterArgs.add("ow");
			
			if(phix!=null){
				filterArgs.add("ref="+phix);
				filterArgs.add("k="+31);
			}
			
			Collections.reverse(filterArgs);
			
			BBDukF.main(filterArgs.toArray(new String[0]));
		}
		
		{
			adapterArgs.add("in="+outFilter);
			adapterArgs.add("outa="+outAdapter);
			
			adapterArgs.add("ow");
			
			Collections.reverse(adapterArgs);
			
			BBMerge.main(adapterArgs.toArray(new String[0]));
		}
		
		{
			trimArgs.add("in="+outFilter);
			trimArgs.add("ref="+outAdapter);
			trimArgs.add("out="+outTrimmed);
			
			ArrayList<byte[]> alb=Tools.toAdapterList(outAdapter, 23);
			boolean doTrim=false;
			if(alb!=null){
				for(byte[] array : alb){
					if(array.length>=11){
						doTrim=true;
					}
				}
			}
			
			if(doTrim){
				trimArgs.add("k=23");
				trimArgs.add("mink=11");
				trimArgs.add("hdist=1");
				trimArgs.add("ktrim=r");
			}

			trimArgs.add("qtrim=r");
			trimArgs.add("trimq=10");
			
			trimArgs.add("tbo");
			trimArgs.add("tpe");
			trimArgs.add("ow");
			trimArgs.add("minlen=62");
			
			Collections.reverse(trimArgs);
			
			BBDukF.main(trimArgs.toArray(new String[0]));
			if(deleteTemp){delete(outFilter);}
		}
		
		{
			mergeArgs.add("in="+outTrimmed);
			
			mergeArgs.add("out="+outMerged);
			mergeArgs.add("outu="+outUnmerged);

			mergeArgs.add("k=75");
			mergeArgs.add("extend2=120");
			mergeArgs.add("rem");
			mergeArgs.add("ecct");
//			mergeArgs.add("loose");
			mergeArgs.add("adapters="+outAdapter);
			mergeArgs.add("ow");
			
			Collections.reverse(mergeArgs);
			
			BBMerge.main(mergeArgs.toArray(new String[0]));
			if(deleteTemp){delete(outTrimmed, outAdapter);}
		}
		
		{
			eccArgs.add("in="+outMerged+","+outUnmerged);
			eccArgs.add("out="+outMergedEcc+","+outUnmergedEcc);

			eccArgs.add("k=50");
			eccArgs.add("ecc");
			eccArgs.add("tossjunk");
			eccArgs.add("deadzone=2");
			eccArgs.add("ow");
			
			Collections.reverse(eccArgs);
			
			Tadpole.main(eccArgs.toArray(new String[0]));
			if(deleteTemp){delete(outMerged, outUnmerged);}
		}
		
		{
			extendArgs.add("in="+outMergedEcc+","+outUnmergedEcc);
			extendArgs.add("out="+outMergedExtended+","+outUnmergedExtended);

			extendArgs.add("k=81");
			extendArgs.add("mode=extend");
			extendArgs.add("el=100");
			extendArgs.add("er=100");
			extendArgs.add("ow");
			extendArgs.add("tossjunk");
			extendArgs.add("deadzone=0");
			extendArgs.add("ecc");
			
			Collections.reverse(extendArgs);
			
			Tadpole.main(extendArgs.toArray(new String[0]));
			if(deleteTemp){delete(outMergedEcc, outUnmergedEcc);}
		}
		
		if(extend2){
			extend2Args.add("in="+outMergedExtended+","+outMergedExtended);
			extend2Args.add("out="+outMergedExtended2+","+outUnmergedExtended2);

			extend2Args.add("k=124");
			extend2Args.add("mode=extend");
			extend2Args.add("el=60");
			extend2Args.add("er=60");
			extend2Args.add("ow");
			extend2Args.add("tossjunk");
			extend2Args.add("deadzone=0");
			extend2Args.add("ecc");
			
			Collections.reverse(extend2Args);
			
			Tadpole.main(extend2Args.toArray(new String[0]));
			if(deleteTemp){delete(outMergedExtended, outMergedExtended);}
		}
		
		{
			if(extend2){
				assembleArgs.add("in="+outMergedExtended2+","+outUnmergedExtended2);
			}else{
				assembleArgs.add("in="+outMergedExtended+","+outUnmergedExtended);
			}
			assembleArgs.add("out="+outMultiK);
			assembleArgs.add("outfinal="+out);

			assembleArgs.add("k=210,250,290");
//			assembleArgs.add("quitearly");
			assembleArgs.add("expand");
			assembleArgs.add("bisect");
			if(deleteTemp){assembleArgs.add("delete");}
			
			Collections.reverse(assembleArgs);
			
			int best=TadpoleWrapper.process(assembleArgs.toArray(new String[0]));
			if(deleteTemp){delete(outMergedExtended, outUnmergedExtended);}
		}
		
//		outstream.println();
//		AssemblyStats2.main(new String[] {"in="+out});
		
//		{
//			assembleArgs.add("in="+outMergedExtended+","+outUnmergedExtended);
//			assembleArgs.add("out="+out);
//
//			assembleArgs.add("k=250");
//			
//			Collections.reverse(assembleArgs);
//			
//			Tadpole.main(assembleArgs.toArray(new String[0]));
//			if(deleteTemp){delete(outMergedExtended, outUnmergedExtended);}
//		}
		
	}
	
	/**
	 * Delete all non-null filenames.
	 * @param prefix Append this prefix to filenames before attempting to delete them
	 * @param names Filenames to delete
	 */
	private void delete(String...names){
		for(String s : names){
			if(s!=null){
				if(verbose){System.err.println("Trying to delete "+s);}
				File f=new File(s);
				if(f.exists()){
					f.delete();
				}
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private ArrayList<String> filterArgs=new ArrayList<String>();
	
	private ArrayList<String> adapterArgs=new ArrayList<String>();
	
	private ArrayList<String> trimArgs=new ArrayList<String>();
	
	private ArrayList<String> assembleArgs=new ArrayList<String>();
	
	private ArrayList<String> mergeArgs=new ArrayList<String>();
	
	private ArrayList<String> eccArgs=new ArrayList<String>();
	
	private ArrayList<String> extendArgs=new ArrayList<String>();
	
	private ArrayList<String> extend2Args=new ArrayList<String>();
	
	private String in1, in2;
	
	private String out="contigs.fa";
	
	private String tempdir=Shared.TMPDIR;

	private boolean deleteTemp=true;
	private boolean gz=false;
	
	private boolean extend2=false;
	
	private String phix="/global/dna/shared/rqc/ref_databases/qaqc/databases/phix174_ill.ref.fa";
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print status messages to this output stream */
	private PrintStream outstream=System.err;
	/** Print verbose messages */
	public static boolean verbose=false;
	/** True if an error was encountered */
	public boolean errorState=false;
	
}
