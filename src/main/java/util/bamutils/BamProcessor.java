package main.java.util.bamutils;

import main.java.util.ThisIsMadnessException;
import main.java.util.bamutils.BamIterator;;

public class BamProcessor extends BamIterator{

	@Override
	public String ProcessBam(String inputBAM, String outputBAM) throws ThisIsMadnessException {
		return ProcessBam(inputBAM, outputBAM, false);
	}

	@Override
	public String ProcessBam(String inputBAM) throws ThisIsMadnessException {
		String[] InputSplit = inputBAM.split(".");
		String outputBAM = "";
		for(int i=0; i<InputSplit.length - 1; i++){
			outputBAM += InputSplit[i] + ".";
		}
		outputBAM += ".out.bam";
		return ProcessBam(inputBAM, outputBAM, false);
	}

	@Override
	public String ProcessBam(String inputBAM, String outputBAM,
			boolean deleteInput) throws ThisIsMadnessException {
		throw new ThisIsMadnessException("This is meant to be the base class for processing BAMs. Do not run this directly!");
	}

	@Override
	public String ProcessBam(String inputBAM, boolean deleteInput) throws ThisIsMadnessException {
		String[] InputSplit = inputBAM.split(".");
		String outputBAM = "";
		for(int i=0; i<InputSplit.length - 1; i++){
			outputBAM += InputSplit[i] + ".";
		}
		outputBAM += ".out.bam";
		return ProcessBam(inputBAM, outputBAM, deleteInput);
	}
}