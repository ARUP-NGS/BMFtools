package main.java.util.bamutils;

import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.SAMFileWriter;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import main.java.util.ThisIsMadnessException;
import main.java.util.Utilities;

/*
 * 
 * @author daniel
 */

public class SecondarySupplementaryFilter extends BamFilter {

	@Override
	public String FilterBam(String inputBAM, String passBAM, String failBAM,
			boolean deleteInput) throws ThisIsMadnessException {
		
		// Open sam/bam file.
		final SamReader reader = SamReaderFactory.makeDefault().open(
				new File(inputBAM));
		SAMFileWriterFactory WriterFactory = new SAMFileWriterFactory();
		SAMFileWriter PassHandle = WriterFactory.makeBAMWriter(reader.getFileHeader(), !reader.getFileHeader().getSortOrder().toString().equals("unsorted"), new File(passBAM));
		SAMFileWriter FailHandle = WriterFactory.makeBAMWriter(reader.getFileHeader(), !reader.getFileHeader().getSortOrder().toString().equals("unsorted"), new File(failBAM));
			
		// Get family sizes
		for (final SAMRecord samRecord : reader) {
				if (!samRecord.getSupplementaryAlignmentFlag() && !samRecord.getNotPrimaryAlignmentFlag()) {
					PassHandle.addAlignment(samRecord);
			}
				else{
					FailHandle.addAlignment(samRecord);
				}
		}
		PassHandle.close();
		FailHandle.close();
		try {
			reader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		Utilities utils = new Utilities();
		if(deleteInput){
			ArrayList<String> commandSet = new ArrayList<String>();
			commandSet.add("rm " + inputBAM);
			utils.Call(commandSet, null);
		}
		return passBAM;
		
	}
	
}