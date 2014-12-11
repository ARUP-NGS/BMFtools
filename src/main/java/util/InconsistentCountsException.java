package main.java.util;

/**
 * These are thrown when the expected number of fastq records do not
 * match between sets which should. 
 * E.g., Read 1 and Read 2 fastq files.
 * @author daniel
 *
 */
public class InconsistentCountsException extends Exception {
	
	public Object source = null;
	/**
	 * 
	 */
	private static final long serialVersionUID = -7795131398099898206L;

	public InconsistentCountsException(String message, Object source) {
		super(message);
		this.source = source;
	}
	
	public InconsistentCountsException(String message) {
		super(message);
		this.source = null; //Clear the source field, as none was provided at last use.
	}
}
