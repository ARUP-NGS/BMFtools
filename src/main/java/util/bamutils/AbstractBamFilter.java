package main.java.util.bamutils;

import main.java.util.ThisIsMadnessException;

abstract class AbstractBamFilter {
	/*
	 */
	public abstract String FilterBam(String inputBAM, String passBAM,
			String failBAM) throws ThisIsMadnessException;

	/*
	 */
	public abstract String FilterBam(String inputBAM, String passBAM,
			String failBAM, boolean deleteInput) throws ThisIsMadnessException;

	/*
	 */
	public abstract String FilterBam(String inputBAM, String passBAM)
			throws ThisIsMadnessException;

	/*
	 */
	public abstract String FilterBam(String inputBAM, String passBAM,
			boolean deleteInput) throws ThisIsMadnessException;

	/*
	 * Abstract class for processing a BAM file, typically by reading in a BAM
	 * file, performing operations, and writing the finished product to another
	 * BAM file. Chooses an output filename based on the input filename.
	 */
	public abstract String FilterBam(String inputBAM)
			throws ThisIsMadnessException;

	/*
	 * Abstract class for processing a BAM file, typically by reading in a BAM
	 * file, performing operations, and writing the finished product to another
	 * BAM file. Chooses an output filename based on the input filename. Boolean
	 * deleteInput determines whether BAMFilter will delete the input BAM
	 * after completion.
	 */
	public abstract String FilterBam(String inputBAM, boolean deleteInput)
			throws ThisIsMadnessException;
}