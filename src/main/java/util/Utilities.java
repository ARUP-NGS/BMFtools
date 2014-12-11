/*
 * Miscellaneous utilities for BMFTools
 * 
 * 
*/
package main.java.util;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Scanner;

import com.google.gson.JsonObject;
import com.google.gson.internal.LinkedTreeMap;
import com.google.gson.reflect.TypeToken;
import com.google.gson.stream.JsonReader;
import com.google.gson.Gson;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParseException;
import com.google.gson.JsonParser;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

public class Utilities {
		String OutputStr = "";
		String ErrorStr = "";
		LinkedTreeMap<String, Object> Output = new LinkedTreeMap<String, Object>();
		

	/*
	 * A generalized call to the system.
	 * @param	call	List of strings for a function call
	 * @param	metadata	Object for holding extra data, included for future extending classes.
	 * @return	Object	Default: LinkedTreeMap.
	 */
	private void Call(List<String> call, Object metadata){
		
		//Resetting output and error strings for each new Call
		this.OutputStr = "";
		this.ErrorStr = "";
	    try {
		      String line;
		      //OutputStream stdin = null;
		      //InputStream stderr = null;
		      //InputStream stdout = null;

		      // launch the command and grab stdin/stdout and stderr
		      ProcessBuilder pb = new ProcessBuilder(call);
		      Process process = pb.start();
		      BufferedReader stdout = new BufferedReader(new  InputStreamReader(process.getInputStream()));
		      BufferedReader stderr = new BufferedReader(new InputStreamReader(process.getErrorStream()));

		      // get process stdout
		      while ((line = stdout.readLine()) != null) {
		        this.OutputStr = this.OutputStr + line + "\n";
		      }
		      stdout.close();

		      // get process stderr
		      while ((line = stderr.readLine()) != null) {
		        this.ErrorStr = this.ErrorStr + line + "\n";
		      }
		      stderr.close();
		      process.destroy();
		    } catch (Exception err) {
		      err.printStackTrace();
		    }
	    Output.put("out", OutputStr);
	    Output.put("err", ErrorStr);
	}

	/*
	 * Tool for reading JSONs
	 */

	public class FileParser {

	/*
	 * @param	filename	Path to file for reading in.
	 */
	    public String[] readLines(String filename) throws IOException {
	        FileReader fileReader = new FileReader(filename);
	        BufferedReader bufferedReader = new BufferedReader(fileReader);
	        List<String> lines = new ArrayList<String>();
	        String line = null;
	        while ((line = bufferedReader.readLine()) != null) {
	            lines.add(line);
	        }
	        bufferedReader.close();
	        return lines.toArray(new String[lines.size()]);
	    }
	   
	    
	    /*
	    * Parses a configuration file. Lines can be commented out with a #.
	    * First token becomes key, second value.
	    * Additional tokens are ignored. 
	    */
	    public LinkedTreeMap<String, String> parseConfig(String filename) throws IOException {
	    	LinkedTreeMap<String, String> config = new LinkedTreeMap<String, String>();
	    	String[] lines = readLines(filename);
	    	for(String line: lines){
	    		if(!line.startsWith("#"))
	    			config.put(line.split("=")[0], line.split("=")[1]);
	    	}
	    	return config;
	    }
	}
	
	public LinkedTreeMap<String, Object> ParseRunJson(String jsonPath) throws JsonParseException, IOException {
		Scanner scanner = new Scanner(new File(jsonPath));
		String jsonStr = scanner.useDelimiter("\\Z").next();
		Gson gson = new Gson();
		JsonReader reader = new JsonReader(new StringReader(jsonStr));
		reader.setLenient(true);
		LinkedTreeMap<String, Object> map = new LinkedTreeMap<String, Object>();
		ValueErrorThrower Pitcher = new ValueErrorThrower();
		scanner.close();
		try {
			map = gson.fromJson(reader, new TypeToken<LinkedTreeMap<String, Object>>() {}.getType());
		} catch (JsonParseException e) {
			e.printStackTrace();
			Pitcher.ValueError("The JSON was not successfully parsed. Check your path and file for validity.");
		}
		System.out.println("Json string: " + jsonStr);
		System.out.println("Map: " + map);
		return map;
	}
	
	public class ValueErrorThrower{
		public void ValueError(String message){
			String fp = "............................................________" + "\n";
			fp += "....................................,.-'\"...................``~.," + "\n";
			fp += ".............................,.-\"...................................-.," + "\n";
			fp += ".........................,/...............................................," + "\n";
			fp += ".....................,?......................................................," + "\n";
			fp += ".................../...........................................................,}" + "\n";
			fp += "................./......................................................,:`^`..}" + "\n";
			fp += ".............../...................................................,:\"........./" + "\n";
			fp += "..............?.....__.........................................:`.........../" + "\n";
			fp += "............./__.(.....\"~-,_..............................,:`........../" + "\n";
			fp += ".........../(_....\"~,_........\"~,_....................,:`........_/" + "\n";
			fp += "..........{.._$;_......\"=,_.......\"-,_.......,.-~-,},.~\";/....}" + "\n";
			fp += "...........((.....*~_.......\"=-._......\";,,./`..../\"............../" + "\n";
			fp += "...,,,___.`~,......\"~.,....................`.....}............../"  + "\n";
			fp += "............(....`=-,,.......`........................(......;_,,-" + "\n";
			fp += "............/.`~,......`-...................................../" + "\n";
			fp += ".............`~.*-,.....................................|,./.....,__" + "\n";
			fp += ",,_..........}.>-._...................................|..............`=~-," + "\n";
			fp += ".....`=~-,__......`,................................." + "\n";
			fp += "...................`=~-,,.,..............................." + "\n";
			fp += "................................`:,,...........................`..............__" + "\n";
			fp += ".....................................`=-,...................,%`>--==``" + "\n";
			fp += "........................................_..........._,-%.......`" + "\n";
			fp += "..................................., " + "\n";
			System.out.println(fp);
			throw new IllegalArgumentException(message);
		}

	}	

}
