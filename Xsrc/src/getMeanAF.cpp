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
    BamReader reader;
    if(argc == 2) {
        inputFilename = argv[1];
    }
    else if(argc > 2) {
        cerr << "Too many arguments" << endl;
        return 1;
    }
    else {
        printf("Description: Calculates the mean aligned fraction for a BAM.\n\n");
        printf("Usage: \n");
        printf("First positional argument: input BAM\n");
        printf("That is all.\n");
        return 1;
    }
    if(!reader.Open(inputFilename)) {
        cerr << "Could not open input BAM" << endl;
        return 1;
    }
    const SamHeader header = reader.GetHeader();
    const RefVector references = reader.GetReferenceData();
    BamAlignment rec;
    AFSum = 0.;
    numReads = 0;

    while (reader.GetNextAlignment(rec)) {
        if(rec.IsMapped()) {
            rec.GetTag("AF", AF);
            numReads++;
            AFSum += AF;
        }
    }
    meanAF = AFSum / numReads;
    printf("Mean Aligned Fraction for BAM %s: %f\n", inputFilename.c_str(), meanAF);
    
    reader.Close();
    return 0;
}
