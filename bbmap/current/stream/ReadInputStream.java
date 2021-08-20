package stream;

import java.util.ArrayList;

public abstract class ReadInputStream {
	

	public abstract Read next();
	
//	public final ArrayList<Read> fetchAll(){
//		ArrayList<Read> out=new ArrayList<Read>();
//		for(ArrayList<Read> list=nextList(); list!=null && list.size()>0; list=nextList()){
//			out.addAll(list);
//		}
//		close();
//		return out;
//	}
	
	public abstract ArrayList<Read> nextList();
	
	public abstract boolean hasMore();

	public abstract void restart();
	
	/** Returns true if there was an error, false otherwise */
	public abstract boolean close();

	public abstract boolean paired();

	protected static final ArrayList<Read> toList(Read[] array){
		if(array==null || array.length==0){return null;}
		ArrayList<Read> list=new ArrayList<Read>(array.length);
		for(int i=0; i<array.length; i++){list.add(array[i]);}
		return list;
	}
	
	/** Return true if this stream has detected an error */
	public boolean errorState(){return errorState;}
	/** TODO */
	protected boolean errorState=false;
	
	public final boolean preferLists(){return true;}

	public abstract void start();
	
}
