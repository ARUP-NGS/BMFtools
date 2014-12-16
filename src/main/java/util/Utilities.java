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
import java.util.logging.Logger;

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

import main.java.bmftools.BMFToolsMain;

public class Utilities {
	String OutputStr = "";
	String ErrorStr = "";
	LinkedTreeMap<String, Object> Output = new LinkedTreeMap<String, Object>();

	/*
	 * A generalized call to the system.
	 * 
	 * @param call List of strings for a function call
	 * 
	 * @param metadata Object for holding extra data, included for future
	 * extending classes.
	 * 
	 * @return Object Default: LinkedTreeMap.
	 */
	private void Call(List<String> call, Object metadata) {

		// Resetting output and error strings for each new Call
		this.OutputStr = "";
		this.ErrorStr = "";
		try {
			String line;
			// OutputStream stdin = null;
			// InputStream stderr = null;
			// InputStream stdout = null;

			// launch the command and grab stdin/stdout and stderr
			ProcessBuilder pb = new ProcessBuilder(call);
			Process process = pb.start();
			BufferedReader stdout = new BufferedReader(new InputStreamReader(
					process.getInputStream()));
			BufferedReader stderr = new BufferedReader(new InputStreamReader(
					process.getErrorStream()));

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
		 * @param filename Path to file for reading in.
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
		 * First token becomes key, second value. Additional tokens are ignored.
		 */
		public LinkedTreeMap<String, String> parseConfig(String filename)
				throws IOException {
			LinkedTreeMap<String, String> config = new LinkedTreeMap<String, String>();
			String[] lines = readLines(filename);
			for (String line : lines) {
				if (!line.startsWith("#"))
					config.put(line.split("=")[0], line.split("=")[1]);
			}
			return config;
		}
	}

	public LinkedTreeMap<String, Object> ParseRunJson(String jsonPath)
			throws JsonParseException, IOException, ThisIsMadnessException {
		Scanner scanner = new Scanner(new File(jsonPath));
		String jsonStr = scanner.useDelimiter("\\Z").next();
		Gson gson = new Gson();
		JsonReader reader = new JsonReader(new StringReader(jsonStr));
		reader.setLenient(true);
		LinkedTreeMap<String, Object> map = new LinkedTreeMap<String, Object>();
		ValueErrorThrower Pitcher = new ValueErrorThrower();
		scanner.close();
		try {
			map = gson.fromJson(reader,
					new TypeToken<LinkedTreeMap<String, Object>>() {
					}.getType());
		} catch (JsonParseException e) {
			e.printStackTrace();
			Pitcher.ValueError("The JSON was not successfully parsed. Check your path and file for validity.");
		}
		System.out.println("Json string: " + jsonStr);
		System.out.println("Map: " + map);
		return map;
	}

	/**
	 * Execute the given system command in its own process, and wait until the
	 * process has completed to return. If the exit value of the process is not
	 * zero, an OperationFailedException is thrown
	 * 
	 * @param command
	 * @throws OperationFailedException
	 */
	protected void executeCommand(final String command)
			throws ThisIsMadnessException {
		executeCommand(command, false);
	}

	/**
	 * Execute the given system command in its own process, and wait until the
	 * process has completed to return. If the exit value of the process is not
	 * zero, an ThisIsMadnessException in thrown TODO: Add logging
	 * 
	 * @param command
	 * @param permitNonZero
	 *            If true, tolerate nonzero exit values from subordinate
	 *            processes
	 * @throws ThisIsMadnessException
	 */
	protected void executeCommand(final String command,
			final boolean permitNonZero) throws ThisIsMadnessException {
		Runtime r = Runtime.getRuntime();
		final Process p;

		try {
			p = r.exec(command);

			// Weirdly, processes that emits tons of data to their error stream
			// can cause some kind of
			// system hang if the data isn't read. Since BWA and samtools both
			// have the potential to do this
			// we by default capture the error stream here and write it to
			// System.err to avoid hangs. s
			final Thread errConsumer = new StringPipeHandler(
					p.getErrorStream(), System.err);
			errConsumer.start();
			// Apparently, same goes for std out, if we dont capture and
			// redirect it we can encounter a hang
			final Thread outConsumer = new StringPipeHandler(
					p.getInputStream(), System.out);
			outConsumer.start();

			// If runtime is going down, destroy the process so it won't become
			// orphaned
			Runtime.getRuntime().addShutdownHook(new Thread() {
				public void run() {
					// System.err.println("Invoking shutdown thread, destroying task with command : "
					// + command);
					p.destroy();
					errConsumer.interrupt();
					outConsumer.interrupt();
				}
			});

			try {
				int exitVal = p.waitFor();
				if (exitVal != 0) {
					if (permitNonZero == false) {
						throw new ThisIsMadnessException(
								"Task terminated with nonzero exit value : "
										+ System.err.toString()
										+ " command was: " + command, this);
					} else {
						// Logger.getLogger(BMFTools.Pipeline.primaryLoggerName).info("Task terminated with nonzero exit value: "
						// + System.err.toString() + " command was: " +
						// command);
						// Logger.getLogger(Pipeline.primaryLoggerName).info("Settings: Nonzero exit status permitted. Continuing Pipeline.");
					}
				} else {
					// Logger.getLogger(Pipeline.primaryLoggerName).info("Task completed successfully: "
					// + command);
				}
			} catch (InterruptedException e) {
				throw new ThisIsMadnessException("Task was interrupted : "
						+ System.err.toString() + "\n"
						+ e.getLocalizedMessage(), this);
			}

		} catch (IOException e1) {
			throw new ThisIsMadnessException(
					"Task encountered an IO exception : "
							+ System.err.toString() + "\n"
							+ e1.getLocalizedMessage(), this);
		}
	}

	public class ValueErrorThrower {
		public void ValueError(String message) throws ThisIsMadnessException {
			String fp = "............................................________"
					+ "\n";
			fp += "....................................,.-'\"...................``~.,"
					+ "\n";
			fp += ".............................,.-\"...................................-.,"
					+ "\n";
			fp += ".........................,/...............................................,"
					+ "\n";
			fp += ".....................,?......................................................,"
					+ "\n";
			fp += ".................../...........................................................,}"
					+ "\n";
			fp += "................./......................................................,:`^`..}"
					+ "\n";
			fp += ".............../...................................................,:\"........./"
					+ "\n";
			fp += "..............?.....__.........................................:`.........../"
					+ "\n";
			fp += "............./__.(.....\"~-,_..............................,:`........../"
					+ "\n";
			fp += ".........../(_....\"~,_........\"~,_....................,:`........_/"
					+ "\n";
			fp += "..........{.._$;_......\"=,_.......\"-,_.......,.-~-,},.~\";/....}"
					+ "\n";
			fp += "...........((.....*~_.......\"=-._......\";,,./`..../\"............../"
					+ "\n";
			fp += "...,,,___.`~,......\"~.,....................`.....}............../"
					+ "\n";
			fp += "............(....`=-,,.......`........................(......;_,,-"
					+ "\n";
			fp += "............/.`~,......`-...................................../"
					+ "\n";
			fp += ".............`~.*-,.....................................|,./.....,__"
					+ "\n";
			fp += ",,_..........}.>-._...................................|..............`=~-,"
					+ "\n";
			fp += ".....`=~-,__......`,................................."
					+ "\n";
			fp += "...................`=~-,,.,..............................."
					+ "\n";
			fp += "................................`:,,...........................`..............__"
					+ "\n";
			fp += ".....................................`=-,...................,%`>--==``"
					+ "\n";
			fp += "........................................_..........._,-%.......`"
					+ "\n";
			fp += "..................................., " + "\n";
			System.out.println(fp);
			throw new ThisIsMadnessException(message);
		}

	}

}
