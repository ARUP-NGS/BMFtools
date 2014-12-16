package main.java.bmftools;

import java.io.IOException;
import java.util.ArrayList;
import java.util.logging.*;

import com.google.gson.internal.LinkedTreeMap;

import main.java.util.StringPipeHandler;
import main.java.util.ThisIsMadnessException;
import main.java.util.Utilities;
import main.java.util.Utilities.ValueErrorThrower;

public class BMFAnalysis {
	public static final String BEDFILE = "bedfile";

	public String bedfile = "default";
	public String LoggerName = "default";
	public boolean i7_index = false;
	public boolean paired = true;
	LinkedTreeMap<String, Object> FastqSettings = null;
	LinkedTreeMap<String, Object> GlobalSettings = null;

	public BMFAnalysis(String LoggerName,
			LinkedTreeMap<String, Object> runProtocol,
			ArrayList<String> FastqAL) throws ThisIsMadnessException {
		Utilities deepMagic = new Utilities();
		ValueErrorThrower VET = deepMagic.new ValueErrorThrower();
		LinkedTreeMap<String, Object> GlobalSettings = (LinkedTreeMap<String, Object>) runProtocol.get("Config");
		String bedfile = ((String) GlobalSettings.get("bedfile"));
		LinkedTreeMap<String, Object> FastqSettings = (LinkedTreeMap<String, Object>) GlobalSettings
				.get("fastqConfig");
		this.i7_index = Boolean.parseBoolean((String) FastqSettings
				.get("i7_index"));
		this.i7_index =
		this.paired = Boolean.parseBoolean((String) GlobalSettings
				.get("paired"));
		if (bedfile == null)
			VET.ValueError("Bed file required.");

		String homing = (String) FastqSettings.get("homing");
		if (!i7_index && homing == null) {
			VET.ValueError("Neither i7 index nor a homing sequence were provided. Abort mission!");
		}
		System.out.println(Boolean.toString(i7_index));

		LinkedTreeMap<String, Object> BamSettings = (LinkedTreeMap<String, Object>) GlobalSettings
				.get("bamConfig");

		int expectedFastqCount = 0;
		if (paired)
			expectedFastqCount += 2;
		else
			expectedFastqCount++;
		if (i7_index)
			expectedFastqCount += 1;
		if (FastqAL.size() != expectedFastqCount) {
			VET.ValueError("Expected " + expectedFastqCount
					+ " fastq files but received " + FastqAL.size() + "!");
		}

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
	 * zero, an ThisIsMadnessException in thrown
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
						Logger.getLogger(this.LoggerName).info(
								"Task terminated with nonzero exit value: "
										+ System.err.toString()
										+ " command was: " + command);
						Logger.getLogger(this.LoggerName)
								.info("Settings: Nonzero exit status permitted. Continuing Pipeline.");
					}
				} else {
					Logger.getLogger(this.LoggerName).info(
							"Task completed successfully: " + command);
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

	public void setBedfile(String bedfile) {
		this.bedfile = bedfile;
	}

	public String getBedfile() {
		return this.bedfile;
	}

	public void setLoggerName(String LoggerName) {
		this.LoggerName = LoggerName;
	}

	public String getLogger() {
		return this.LoggerName;
	}
}