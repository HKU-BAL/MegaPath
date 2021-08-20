package align2;

import java.util.Comparator;

import stream.Read;

/**
 * @author Brian Bushnell
 * @date Oct 27, 2014
 *
 */

public final class ReadComparatorName implements Comparator<Read>{
	
	@Override
	public int compare(Read r1, Read r2) {
		
		if(r1.id==null && r2.id==null){return r1.pairnum()-r2.pairnum();}
		if(r1.id==null){return -1;}
		if(r2.id==null){return 1;}
		int x=r1.id.compareTo(r2.id);
		if(x==0){return r1.pairnum()-r2.pairnum();}
		return x;
	}

	public static final ReadComparatorName comparator=new ReadComparatorName();
	
}
