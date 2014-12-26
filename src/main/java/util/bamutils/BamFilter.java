package main.java.util.bamutils;

import main.java.util.ThisIsMadnessException;

public class BamFilter extends AbstractBamFilter{

	@Override
	public String FilterBam(String inputBAM, String passBAM, String failBAM,
			boolean deleteInput) throws ThisIsMadnessException {
		throw new ThisIsMadnessException("This is meant to be the base BAM Filter class. This must be run through inheritance.");
	}
	
	@Override
	public String FilterBam(String inputBAM, String passBAM) throws ThisIsMadnessException {
		String[] InputSplit = inputBAM.split(".");
		String failBAM = "";
		for(int i=0; i<InputSplit.length - 1; i++){
			failBAM += InputSplit[i] + ".";
		}
		failBAM += ".fail.bam";
		return FilterBam(inputBAM, passBAM, failBAM, false);
	}

	@Override
	public String FilterBam(String inputBAM) throws ThisIsMadnessException {
		String[] InputSplit = inputBAM.split(".");
		String passBAM = "";
		String failBAM = "";
		for(int i=0; i<InputSplit.length - 1; i++){
			passBAM += InputSplit[i] + ".";
			failBAM += InputSplit[i] + ".";
		}
		passBAM += ".pass.bam";
		failBAM += ".fail.bam";
		return FilterBam(inputBAM, passBAM, failBAM, false);
	}

	@Override
	public String FilterBam(String inputBAM, String passBAM,
			boolean deleteInput) throws ThisIsMadnessException {
		String[] InputSplit = inputBAM.split(".");
		String failBAM = "";
		for(int i=0; i<InputSplit.length - 1; i++){
			failBAM += InputSplit[i] + ".";
		}
		failBAM += ".fail.bam";
		return FilterBam(inputBAM, passBAM, failBAM, deleteInput);
	}

	@Override
	public String FilterBam(String inputBAM, boolean deleteInput) throws ThisIsMadnessException {
		String[] InputSplit = inputBAM.split(".");
		String passBAM = "";
		String failBAM = "";
		for(int i=0; i<InputSplit.length - 1; i++){
			passBAM += InputSplit[i] + ".";
			failBAM += InputSplit[i] + ".";
		}
		passBAM += ".pass.bam";
		failBAM += ".fail.bam";
		return FilterBam(inputBAM, passBAM, failBAM, deleteInput);
	}

	@Override
	public String FilterBam(String inputBAM, String passBAM, String failBAM)
			throws ThisIsMadnessException {
		return FilterBam(inputBAM, passBAM, failBAM, false);
	}

}