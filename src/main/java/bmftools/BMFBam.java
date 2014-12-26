package main.java.bmftools;

import java.io.File;
import java.util.ArrayList;
import java.util.NoSuchElementException;

import main.java.util.InconsistentCountsException;

import org.apache.commons.lang.StringUtils;

import com.google.gson.internal.LinkedTreeMap;

/*
 * A Java re-write of BMFTools from Pipeline's BCBam.py
 */

public class BMFBam {

	public static final String ABRA_PATH = "abra.path";
	public String AbraPath = "/mounts/bin/abra-0.86-SNAPSHOT-jar-with-dependencies.jar";
	
	public void Abracadabra(File inbam, File outbam, LinkedTreeMap<String, Object> GlobalSettings){
		
		String AbraAttr = (String) GlobalSettings.get(ABRA_PATH);
		if(AbraAttr != null)
			this.AbraPath = AbraAttr;
	/*
	 *     command = ("java {} -jar {} --in {}".format(memStr, jar, inbam) +
               " --out {} --ref {} --targets".format(outbam, ref) +
               " {} --threads {} ".format(bed, threads) +
               "--working {}".format(working))
    pl("Command: {}.".format(command))	
	 */
		
		
	}
	
}