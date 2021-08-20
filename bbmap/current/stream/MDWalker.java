package stream;

/**
 * @author Brian Bushnell
 * @date May 5, 2016
 * 
 *
 */
public class MDWalker {
	
	MDWalker(String tag){
		mdTag=tag;
		pos=(mdTag.startsWith("MD:Z:") ? 5 : 0);
		
		mpos=0;
		bpos=0;
		rpos=0;
		sym=0;
		current=0;
		mode=0;
	}
	
	boolean nextSub(){
		sym=0;
		while(pos<mdTag.length()){
			char c=mdTag.charAt(pos);
			pos++;
			
			if(Character.isDigit(c)){
				current=(current*10)+(c-'0');
				mode=NORMAL;
			}else{
				if(current>0){
					bpos+=current;
					rpos+=current;
					mpos+=current; 
					assert(mode==NORMAL) : mode+", "+current;
					current=0;	
				}
				if(c=='^'){mode=DEL;}
				else if(mode==DEL){
					rpos++;
					mpos++;
					sym=c;
				}else if(mode==NORMAL || mode==SUB){
					mode=SUB;
					bpos++;
					rpos++;
					mpos++;
					sym=c;
					return true;
				}
			}
		}
		return false;
	}
	
	public int matchPosition(){
		return mpos-1;
	}
	
	public int basePosition(){
		return bpos-1;
	}
	
	public int refPosition(){
		return rpos-1;
	}
	
	public char symbol(){
		assert(sym!=0);
		return sym;
	}

	/** Position in match string (excluding clipping and insertions) */
	private int mpos;
	/** Position in read bases (excluding clipping and insertions) */
	private int bpos;
	/** Position in reference bases (excluding clipping) */
	private int rpos;
	private char sym;
	
	private String mdTag;
	private int pos;
	private int current;
	private int mode;
	
//	private int dels=0, subs=0, normals=0;
	private static final int NORMAL=0, SUB=1, DEL=2;
	
}
