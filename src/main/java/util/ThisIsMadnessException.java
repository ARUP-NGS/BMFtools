package main.java.util;

/**
 * This is the base exception class for utility failures.
 * E.g., parsing of files, execution of commands, etc.
 * @author daniel
 *
 */
public class ThisIsMadnessException extends Exception {
	
	public Object source = null;
	/**
	 * 
	 */
	private static final long serialVersionUID = 1195131273499482206L;

	public ThisIsMadnessException(String message, Object source) {
		super(message);
		this.source = source;
	}
	
	public ThisIsMadnessException(String message) {
		super(message);
		this.source = null; //Clear the source field, as none was provided at last use.
	}
}
