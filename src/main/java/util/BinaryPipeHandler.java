package main.java.util;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

/**
 * 
 * --- Taken from ARUP-NGS Pipeline, 12/16/14 ---
 * 
 * This class uses a thread to read binary data from an input stream and write it to an output stream
 * Its used in Pipeline to read the data emitted to stdout (or stderr) by a process and write it
 * to a file. Without this running as a separate thread, buffers used to store data from stdout will
 * fill up and the process generating the data may hang. 
 * @author brendan
 *
 */
public class BinaryPipeHandler extends Thread {

	InputStream inpStr;
	OutputStream writer;

	public BinaryPipeHandler(InputStream inpStr, OutputStream writer) {
		this.inpStr = inpStr;
		this.writer = writer;
		this.setPriority(Thread.MAX_PRIORITY);
	}

	public void run() {
			try {
				//Attempt to read data in big chunks, this is a lot faster than
				//doing things one byte at a time
				//If an application writes binary data to stdout extremely quickly
				//then we may need to increase the buffer size, but this appears to be OK for now
				byte[] data = new byte[1024*1024]; //One megabyte
				int bytesRead = inpStr.read(data);
				int count = 0;
				while(bytesRead >= 0) {
					writer.write(data, 0, bytesRead);
					bytesRead = inpStr.read(data);
					count++;
				}
				writer.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace(System.err);
			}

	}
}
