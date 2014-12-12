package main.java.util;

abstract class BamProcessor {
	/*
	 * Abstract class for processing a BAM file, typically by reading in a BAM
	 * file, performing operations, and writing the finished product to another
	 * BAM file. Explicitly requires outputBAM in call.
	 */
	public abstract String ProcessBam(String inputBAM, String outputBAM);

	/*
	 * Abstract class for processing a BAM file, typically by reading in a BAM
	 * file, performing operations, and writing the finished product to another
	 * BAM file. Chooses an output filename based on the input filename.
	 */
	public abstract String ProcessBam(String inputBAM);

	/*
	 * Abstract class for processing a BAM file, typically by reading in a BAM
	 * file, performing operations, and writing the finished product to another
	 * BAM file. Explicitly requires outputBAM in call. Boolean deleteInput
	 * determines whether BAMProcessor will delete the input BAM after
	 * completion.
	 */
	public abstract String ProcessBam(String inputBAM, String outputBAM,
			boolean deleteInput);

	/*
	 * Abstract class for processing a BAM file, typically by reading in a BAM
	 * file, performing operations, and writing the finished product to another
	 * BAM file. Chooses an output filename based on the input filename. Boolean
	 * deleteInput determines whether BAMProcessor will delete the input BAM
	 * after completion.
	 */
	public abstract String ProcessBam(String inputBAM, boolean deleteInput);
}