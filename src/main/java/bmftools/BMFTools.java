/*
*/

package main.java.bmftools;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.impl.Arguments;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;

// TODO: Auto-generated Javadoc
//import org.apache.log4j.Logger;


/**
 * BMFTools: Main Class
 * A developing suite of tools for ultra-low frequency variant detection and analysis.
 * Currently supported: SNP calling.
 * Experimental: Translocation
 * TODO: CNV detection
 * @author daniel
 */
