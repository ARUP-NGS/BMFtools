package main.java.util;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.TreeSet;
import java.io.FileWriter;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.seekablestream.SeekableStream;

public class BamStats {

	public static SeekableStream myIndexSeekableStream() {
		throw new UnsupportedOperationException();
	}

	public void InsertSizeHistogram(String inputBamfile, String HistogramFile)
			throws IOException {
		
		// Open sam/bam file.
		final SamReader reader = SamReaderFactory.makeDefault().open(
				new File(inputBamfile));
		ArrayList<Integer> InsertSizes = new ArrayList<Integer>();
		
		// Get insert size
		for (final SAMRecord samRecord : reader) {
			InsertSizes.add(Math.abs(samRecord.getInferredInsertSize()));
		}

		// Get set of unique sizes for histogram.
		TreeSet<Integer> InsertSizeSet = new TreeSet<Integer>(InsertSizes);

		// Create a dictionary
		LinkedHashMap<Integer, Integer> InsertSizeDict = new LinkedHashMap<Integer, Integer>();
		FileWriter out = new FileWriter(HistogramFile);
		
		// Write histogram to file.
		for (Integer size : InsertSizeSet) {
			InsertSizeDict.put(size, Collections.frequency(InsertSizes, size));
			out.write(size + "\t" + Collections.frequency(InsertSizes, size));
		}
		out.close();
		reader.close();

	}
	

	public void FamilySizeHistogram(String inputBamfile, String FamilySizeHistogramFile)
			throws IOException {
		
		// Open sam/bam file.
		final SamReader reader = SamReaderFactory.makeDefault().open(
				new File(inputBamfile));
		ArrayList<Integer> FamilySizes = new ArrayList<Integer>();
		
		// Get family sizes
		for (final SAMRecord samRecord : reader) {
			FamilySizes.add(samRecord.getIntegerAttribute("FM"));
		}

		// Get set of unique sizes for histogram, sorted. (TreeSets are sorted.)
		TreeSet<Integer> FamilySizeSet = new TreeSet<Integer>(FamilySizes);

		// Create a dictionary
		LinkedHashMap<Integer, Integer> FamilySizeDict = new LinkedHashMap<Integer, Integer>();
		FileWriter out = new FileWriter(FamilySizeHistogramFile);
		
		// Write histogram to file.
		for (Integer size : FamilySizeSet) {
			FamilySizeDict.put(size, Collections.frequency(FamilySizes, size));
			out.write(size + "\t" + Collections.frequency(FamilySizes, size) + "\n");
		}
		out.close();
		reader.close();
	}
	
	

	public void BamBarcodeIndex(String inputBamfile, String BarcodeIndexFile)
			throws IOException {
		
		// Open sam/bam file.
		final SamReader reader = SamReaderFactory.makeDefault().open(
				new File(inputBamfile));
		ArrayList<String> Barcodes = new ArrayList<String>();
		
		// Get family sizes
		for (final SAMRecord samRecord : reader) {
			Barcodes.add(samRecord.getStringAttribute("BS"));
		}

		// Get set of unique sizes for histogram, sorted. (TreeSets are sorted.)
		TreeSet<String> FamilySizeSet = new TreeSet<String>(Barcodes);

		// Create a dictionary
		LinkedHashMap<String, Integer> BarcodeDict = new LinkedHashMap<String, Integer>();
		FileWriter out = new FileWriter(BarcodeIndexFile);
		
		// Write histogram to file.
		for (String barcode : FamilySizeSet) {
			BarcodeDict.put(barcode, Collections.frequency(Barcodes, barcode));
			out.write(Collections.frequency(Barcodes, barcode) + "\t" + barcode + "\n");
		}
		out.close();
		reader.close();
	}
	
	
	/*
	 * Counts the number of occurrences of each barcode without checking for 
	 * passing filters.
	 * Not complete - add processing and output!
	 */
	public void RetagFamilySize(String inputBamfile, String BarcodeIndexFile, String outputBamfile)
			throws IOException {
		
		// Open sam/bam file.
		final SamReader reader = SamReaderFactory.makeDefault().open(
				new File(inputBamfile));
		ArrayList<String> Barcodes = new ArrayList<String>();
		
		// Get family sizes
		for (final SAMRecord samRecord : reader) {
			Barcodes.add(samRecord.getStringAttribute("BS"));
		}

		// Get set of unique sizes for histogram, sorted. (TreeSets are sorted.)
		TreeSet<String> BarcodeSet = new TreeSet<String>(Barcodes);

		// Create a dictionary
		LinkedHashMap<String, Integer> BarcodeDict = new LinkedHashMap<String, Integer>();
		
		// Write histogram to file.
		for (String barcode : BarcodeSet) {
			BarcodeDict.put(barcode, Collections.frequency(Barcodes, barcode));
		}
		reader.close();
	}

}