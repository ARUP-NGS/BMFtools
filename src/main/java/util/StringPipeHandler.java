package main.java.util;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.PrintStream;

/**
 * --- Taken from ARUP-NGS Pipeline, 12/11/14 ---
 * 
 * This class uses a thread to read text data from an input stream and write it to an output stream.
 * It's used in Pipeline to read the data emitted to stdout (or stderr) by a process and write it
 * to a file. Without this running as a separate thread, buffers used to store data from stdout will
 * fill up and the process generating the data may hang. 
 *  For binary data, use a BinaryPipeHandler
 * @author brendan
 *
 */
public class StringPipeHandler extends Thread {

		InputStream inpStr;
		PrintStream stream;
		
		public StringPipeHandler(InputStream inpStr, OutputStream stream) {
			this.inpStr = inpStr;
			this.stream = new PrintStream(stream);
		}
		
		public StringPipeHandler(InputStream inpStr, PrintStream stream) {
			this.inpStr = inpStr;
			this.stream = stream;
		}

		public void run() {
			try {
				InputStreamReader inpStrd = new InputStreamReader(inpStr);
				BufferedReader buffRd = new BufferedReader(inpStrd);
				String line = null;
				int count = 0;
				while((line = buffRd.readLine()) != null) {
					if (stream != null) {
						stream.println(line);

						count++;
					}
				}
				buffRd.close();

			} catch(Exception e) {
				System.out.println(e);
			}

		}
}
