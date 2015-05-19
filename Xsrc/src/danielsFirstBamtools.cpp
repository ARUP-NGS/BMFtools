//Hey this is a C++ file!

#include <iostream>
#include <string>
#include "api/BamReader.h"
#include "api/BamWriter.h"

using namespace std;
using namespace BamTools;

int main(int argc, char* argv[])
{   
    string inputFilename;
    string outputFilename;
    float AF;
    float minAF;
    if(argc == 3) {
        inputFilename = argv[1];
        outputFilename = argv[2];
    }
    else if(argc > 3) {
        cerr << "Too many arguments" << endl;
        return 1;
    }
    else {
        printf("Description: Calculates the mean aligned fraction for a BAM.\n\n");
        printf("Usage: \n");
        printf("First positional argument: input BAM")
        printf("Second positional argument: output BAM")

    }
    BamReader reader;
    if(!reader.Open(inputFilename)) {
        cerr << "Could not open input BAM" << endl;
        return 1;
    }
    BamWriter writer;
    const SamHeader header = reader.GetHeader();
    const RefVector references = reader.GetReferenceData();
    if(!writer.Open(outputFilename, header, references)) {
        cerr << "Could not open output BAM file for writing." << endl;
        return 1;
    }
    BamAlignment rec;
    while (reader.GetNextAlignment(rec)) {
        
    }
    reader.Close();
    return 0;
}
