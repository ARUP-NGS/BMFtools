package main.java.bmftools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.logging.Logger;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.internal.LinkedTreeMap;
import com.google.gson.reflect.TypeToken;

import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;
import main.java.util.DeepMagic;
import main.java.util.DeepMagic.FileParser;
import main.java.util.DeepMagic.ValueErrorThrower;
import main.java.bmftools.BMFFastq;
import htsjdk;


/**
 * BMFTools is a suite of tools under active development for 
 * ultra-low frequency variant detection and analysis and 
 * removal/accounting of PCR Duplicates
 * Fully operational: SNP detection
 * Also included is an experimental technique (XLIF) for translocation
 * TODO: CNV detection
 * @author daniel
 */

public class BMFTools {

	public static final String BEDFILE = "bedfile";
	
	public String bedfile = "default";
	
	public static void main(String[] arguments) throws IOException{
	Namespace args = null;
	Logger log = Logger.getLogger("Main Log");
	log.info("Parsing arguments.");
    ArgumentParser parser = ArgumentParsers.newArgumentParser("BMFTools")
            .description("Run BMFTools analysis");
    parser.addArgument("--conf").metavar("configFile").type(String.class).required(true).help("File path to config file (required)");
    parser.addArgument("--run", "-r").metavar("Analysis run XML").type(String.class).required(true).help("Path to xml for the run.");
    parser.addArgument("-i", "--input").required(true).nargs("+").type(String.class).help("Path(s) to input fastq.");
	try {
			args = parser.parseArgs(arguments);
		} catch (ArgumentParserException e) {
			parser.handleError(e);
		}
	/*
	 * Parsing configuration files.
	 */
	
	DeepMagic deepMagic = new DeepMagic();
	ValueErrorThrower VET = deepMagic.new ValueErrorThrower();
	FileParser confParser = deepMagic.new FileParser();
	//Currently using only JSON. Will like change in the future.
	//LinkedTreeMap<String, String> configuration = confParser.parseConfig(args.getString("conf"));
	@SuppressWarnings("unchecked")
	LinkedTreeMap<String, Object> runProtocol = deepMagic.ParseRunJson((args.getString("run")));
	Gson gson = new Gson();
	System.out.println(runProtocol.get("Config").getClass());
	System.out.println("Trying to load global settings. Current json string for Config: " + runProtocol.get("Config"));
	LinkedTreeMap<String, Object> GlobalSettings = (LinkedTreeMap<String, Object>) runProtocol.get("Config");
	System.out.println("Successfully loaded global settings.");
	boolean paired = Boolean.parseBoolean((String) GlobalSettings.get("paired"));
	String bedfile = ((String) GlobalSettings.get("bedfile"));
	if(bedfile == null)
		VET.ValueError("Bed file required.");
	
	LinkedTreeMap<String, Object> FastqSettings = (LinkedTreeMap<String, Object>) GlobalSettings.get("fastqConfig");
	String test = (String) FastqSettings.get("i7_index");
	System.out.println("i7_index value: " + test);
	boolean i7_index = Boolean.parseBoolean((String) FastqSettings.get("i7_index"));
	String homing = (String) FastqSettings.get("homing");
	if(!i7_index && homing == null){
		VET.ValueError("Neither i7 index nor a homing sequence were provided. Abort mission!");
	}
	System.out.println(Boolean.toString(i7_index));
	
	LinkedTreeMap<String, Object> BamSettings = (LinkedTreeMap<String, Object>)GlobalSettings.get("bamConfig");
	
	ArrayList<String> FastqAL = args.get("input");
	System.out.println(FastqAL);
	int expectedFastqCount = 0;
	if(paired)
		expectedFastqCount+=2;
	else
		expectedFastqCount++;
	if(i7_index)
		expectedFastqCount += 1;
	if(FastqAL.size() != expectedFastqCount) {
		VET.ValueError("Expected " + expectedFastqCount + " fastq files but received " + FastqAL.size() + "!");
	}
	BMFFastq bf = new BMFFastq();
	if(i7_index){
		MergedFamilies = bf.MergeFamilies(FastqAL);
	}
	}

	public void setBedfile(String bedfile){
		this.bedfile = bedfile;
	}
	
	public String getBedfile(){
		return this.bedfile;
	}
}