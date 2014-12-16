package main.java.bmftools;

import java.io.IOException;
import java.util.ArrayList;
import java.util.logging.Logger;

import main.java.util.ThisIsMadnessException;
import main.java.util.Utilities;
import main.java.util.Utilities.FileParser;
import main.java.util.Utilities.ValueErrorThrower;
import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;

import com.google.gson.Gson;
import com.google.gson.internal.LinkedTreeMap;

import main.java.bmftools.BMFAnalysis;

/**
 * BMFTools is a suite of tools under active development for ultra-low frequency
 * variant detection and analysis and removal/accounting of PCR Duplicates Fully
 * operational: SNP detection Also included is an experimental technique (XLIF)
 * for translocation TODO: CNV detection
 * 
 * @author daniel
 */

public class BMFToolsMain {

	public static final String BEDFILE = "bedfile";

	public static void main(String[] arguments) throws IOException,
			ThisIsMadnessException {
		Namespace args = null;
		Logger log = Logger.getLogger("Main Log");
		log.info("Parsing arguments.");
		ArgumentParser parser = ArgumentParsers.newArgumentParser("BMFTools")
				.description("Run BMFTools analysis");
		parser.addArgument("--conf").metavar("configFile").type(String.class)
				.required(true).help("File path to config file (required)");
		parser.addArgument("--run", "-r").metavar("Analysis run XML")
				.type(String.class).required(true)
				.help("Path to xml for the run.");
		parser.addArgument("-i", "--input").required(true).nargs("+")
				.type(String.class).help("Path(s) to input fastq.");
		parser.addArgument("-l", "--log").metavar("Logger name")
				.type(String.class).required(true)
				.help("String name for Logger");
		String LoggerName = "Default_Loggername";
		try {
			args = parser.parseArgs(arguments);
		} catch (ArgumentParserException e) {
			parser.handleError(e);
		}
		/*
		 * Parsing configuration files.
		 */

		Utilities deepMagic = new Utilities();
		ValueErrorThrower VET = deepMagic.new ValueErrorThrower();
		FileParser confParser = deepMagic.new FileParser();
		// Currently using only JSON. Will like change in the future.
		// LinkedTreeMap<String, String> configuration =
		// confParser.parseConfig(args.getString("conf"));
		@SuppressWarnings("unchecked")
		LinkedTreeMap<String, Object> runProtocol = deepMagic
				.ParseRunJson((args.getString("run")));
		Gson gson = new Gson();
		LinkedTreeMap<String, String> configDict = confParser.parseConfig(args
				.getString("conf"));
		System.out.println(runProtocol.get("Config").getClass());
		System.out
				.println("Trying to load global settings. Current json string for Config: "
						+ runProtocol.get("Config"));
		LinkedTreeMap<String, Object> GlobalSettings = (LinkedTreeMap<String, Object>) runProtocol
				.get("Config");
		System.out.println("Successfully loaded global settings.");
		ArrayList<String> FastqAL = args.get("input");
		BMFAnalysis RunBMF = new BMFAnalysis(args.getString("log"),runProtocol, FastqAL);

	}

}