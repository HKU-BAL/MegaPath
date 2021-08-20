package align2;

import java.util.Comparator;

import stream.Read;

/**
 * @author Brian Bushnell
 * @date Oct 27, 2014
 *
 */

public final class ReadComparatorID implements Comparator<Read>{
	
	@Override
	public int compare(Read r1, Read r2) {
		if(r1.numericID<r2.numericID){return -1;}
		else if(r1.numericID>r2.numericID){return 1;}
		
		if(!r1.id.equals(r2.id)){return r1.id.compareTo(r2.id);}
		return 0;
	}

	public static final ReadComparatorID comparator=new ReadComparatorID();
	
}
