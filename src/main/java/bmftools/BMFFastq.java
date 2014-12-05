package main.java.bmftools;

import java.util.ArrayList;
import htsjdk.samtools.fastq.*;


/*
 * A Java re-write of BMFTools from Pipeline's BCFastq.py
 */

public class BMFFastq {

	/*
	 * Function which merges families of reads together. i7_index
	 * 
	 * @oaram inFastqs An ArrayList of Strings, populated from argparse4j
	 * 
	 * @return output An ArrayList of Strings of final fastqs. 1 (single-end) or
	 * 2 (paired-end)
	 */
	public FastqRecord MergeFamilies(ArrayList<FastqRecord> inFastqs) {
		FastqRecord output = new FastqRecord(null, null, null, null);
		return output;
	}
}