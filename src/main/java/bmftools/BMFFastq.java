package main.java.bmftools;

import java.io.File;
import java.util.ArrayList;
import java.util.NoSuchElementException;

import htsjdk.samtools.fastq.*;
import main.java.util.InconsistentCountsException;
import org.apache.commons.lang.StringUtils;

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
		FastqRecord output = new FastqRecord(inFastqs.get(0).getReadHeader(),
				null, null, null);
		return output;
	}

	/*
	 * Function which tags the original reads with a pass/fail and the barcode
	 * sequence, writing them to provided output files.
	 */
	public String[] TagShadesFastq(String inFastq1, String inFastq2,
			String indexFastq, String outFastq1, String outFastq2)
			throws InconsistentCountsException {
		FastqReader inFq1 = new FastqReader(new File(inFastq1));
		FastqReader inFq2 = new FastqReader(new File(inFastq2));
		FastqReader indexFq = new FastqReader(new File(indexFastq));
		FastqWriterFactory outBuilder = new FastqWriterFactory();
		FastqWriter outFq1 = outBuilder.newWriter(new File(outFastq1));
		FastqWriter outFq2 = outBuilder.newWriter(new File(outFastq2));
		FastqRecord read1 = null;
		FastqRecord read2 = null;
		String barcodeString = null;
		double ratio = 5 / 6.;
		int barcodeLength = -1;
		int HomopolymerMax = -1;
		boolean PassesQC = true;
		while (inFq1.hasNext()) {
			read1 = inFq1.next();
			try {
				read2 = inFq2.next();
			} catch (NoSuchElementException e) {
				inFq1.close();
				inFq2.close();
				indexFq.close();
				outFq1.close();
				outFq2.close();
				e.printStackTrace();
				throw new InconsistentCountsException(
						"Read 1 file and Read 2 file do not have the same number of reads.",
						inFq2.toString());
			}
			try {
				barcodeString = indexFq.next().getReadString();
				if (barcodeLength == -1) {
					barcodeLength = barcodeString.length();
					HomopolymerMax = (int) ((double) barcodeLength * ratio);
				}
				if (barcodeString.contains("N")
						|| barcodeString.contains(StringUtils.repeat("A",
								HomopolymerMax))
						|| barcodeString.contains(StringUtils.repeat("C",
								HomopolymerMax))
						|| barcodeString.contains(StringUtils.repeat("G",
								HomopolymerMax))
						|| barcodeString.contains(StringUtils.repeat("T",
								HomopolymerMax))) {
					PassesQC = false;
				}
				else {
					PassesQC = true;
				}
			} catch (NoSuchElementException e) {
				inFq1.close();
				inFq2.close();
				indexFq.close();
				outFq1.close();
				outFq2.close();
				e.printStackTrace();
				throw new InconsistentCountsException(
						"Index reads files does not have the same number of reads as read 1 and read 2.",
						indexFq.toString());
			}
			
			if(PassesQC){
				outFq1.write(new FastqRecord(read1.getReadHeader() + " #G~IndexPass #G~" + barcodeString, read1.getReadString(), read1.getBaseQualityHeader(), read1.getBaseQualityString()));
				outFq2.write(new FastqRecord(read2.getReadHeader() + " #G~IndexPass #G~" + barcodeString, read2.getReadString(), read2.getBaseQualityHeader(), read2.getBaseQualityString()));
			} else {
				outFq1.write(new FastqRecord(read1.getReadHeader() + " #G~IndexFail #G~" + barcodeString, read1.getReadString(), read1.getBaseQualityHeader(), read1.getBaseQualityString()));
				outFq2.write(new FastqRecord(read2.getReadHeader() + " #G~IndexFail #G~" + barcodeString, read2.getReadString(), read2.getBaseQualityHeader(), read2.getBaseQualityString()));
			}

		}
		inFq1.close();
		inFq2.close();
		indexFq.close();
		outFq1.close();
		outFq2.close();
		return new String[] { outFastq1, outFastq2 };
	}
}