package main.java.util.bamutils;

import main.java.util.ThisIsMadnessException;

abstract class BamIterator {
	/*
	 * Abstract class for processing a BAM file, typically by reading in a BAM
	 * file, performing operations, and writing the finished product to another
	 * BAM file. Explicitly requires outputBAM in call.
	 */
	public abstract String ProcessBam(String inputBAM, String outputBAM) throws ThisIsMadnessException;

	/*
	 * Abstract class for processing a BAM file, typically by reading in a BAM
	 * file, performing operations, and writing the finished product to another
	 * BAM file. Chooses an output filename based on the input filename.
	 */
	public abstract String ProcessBam(String inputBAM) throws ThisIsMadnessException;

	/*
	 * Abstract class for processing a BAM file, typically by reading in a BAM
	 * file, performing operations, and writing the finished product to another
	 * BAM file. Explicitly requires outputBAM in call. Boolean deleteInput
	 * determines whether BAMProcessor will delete the input BAM after
	 * completion.
	 */
	public abstract String ProcessBam(String inputBAM, String outputBAM,
			boolean deleteInput) throws ThisIsMadnessException;

	/*
	 * Abstract class for processing a BAM file, typically by reading in a BAM
	 * file, performing operations, and writing the finished product to another
	 * BAM file. Chooses an output filename based on the input filename. Boolean
	 * deleteInput determines whether BAMProcessor will delete the input BAM
	 * after completion.
	 */
	public abstract String ProcessBam(String inputBAM, boolean deleteInput) throws ThisIsMadnessException;
}