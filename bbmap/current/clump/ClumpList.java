package clump;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicInteger;

import align2.Shared;

import stream.Read;

/**
 * A list of clumps, meaning a list of lists of reads.
 * Allows adding reads by streaming and generating new clumps as needed.
 * The input reads must be correctly ordered.
 * @author Brian Bushnell
 * @date Nov 9, 2015
 *
 */
public class ClumpList extends ArrayList<Clump> {

	public ClumpList(){}
	
	public ClumpList(ArrayList<Read> list){
		addReads(list);
	}
	
	public void addReads(ArrayList<Read> list){
		assert(list.getClass()!=Clump.class) : list.getClass();
		for(final Read r : list){
			final long[] obj=(long[])r.obj;
			final long kmer=obj[0];
			if(kmer!=currentKmer){
				currentKmer=kmer;
				currentClump=new Clump(kmer);
				add(currentClump);
			}
			currentClump.add(r);
		}
	}
	
	public ArrayList<Read> condense(){
		final int threads=Shared.threads();
		return condense(threads);
	}
	
	public ArrayList<Read> condense(final int threads){
		final ArrayList<CondenseThread> alct=new ArrayList<CondenseThread>(threads);
		for(int i=0; i<threads; i++){alct.add(new CondenseThread());}
		
		if(verbose){outstream.println("Starting condense threads.");}
		for(CondenseThread ct : alct){ct.start();}
		
		if(verbose){outstream.println("Waiting for threads.");}
		long readsThisPass=0;
		/* Wait for threads to die */
		for(CondenseThread ct : alct){
			
			/* Wait for a thread to die */
			while(ct.getState()!=Thread.State.TERMINATED){
				try {
					ct.join();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
			readsThisPass+=ct.storage.size();
		}
		
		if(verbose){outstream.println("Gathering reads.");}
		ArrayList<Read> list=new ArrayList<Read>((int)readsThisPass);
		for(int i=0; i<threads; i++){
			CondenseThread ct=alct.set(i, null);
			list.addAll(ct.storage);
		}
		
		assert(list.size()==readsThisPass);
		return list;
	}
	
	@Override
	public void clear(){
		super.clear();
		currentClump=null;
		currentKmer=Long.MIN_VALUE;
		ptr.set(0);
	}
	
	private class CondenseThread extends Thread{
		
		@Override
		public void run(){
			final int size=size();
			for(int i=ptr.getAndIncrement(); i<size; i=ptr.getAndIncrement()){
				Clump c=get(i);
				ArrayList<Read> list=c.makeConsensus();
				storage.addAll(list);
				c.clear();
				set(i, null);
			}
		}
		
		private ArrayList<Read> storage=new ArrayList<Read>();
		
	}

	private Clump currentClump=null;
	private long currentKmer=Long.MIN_VALUE;
	private final AtomicInteger ptr=new AtomicInteger(0);
	
	private static final long serialVersionUID = 1L;
	private static boolean verbose=false;
	private static final PrintStream outstream=System.err;
	
}
