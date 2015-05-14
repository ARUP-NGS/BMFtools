//Hey this is a C++ file!

#include <iostream>
#include <string>
#include <cmath>
#include "api/BamReader.h"
#include "api/BamWriter.h"

using namespace std;
using namespace BamTools;

int main(int argc, char* argv[])
{
    //MaxInsert is the largest insert included in the calculation of mean insert size.
    const long MaxInsert = 2000;
    string inputFilename;
    float meanAF, AF, AFSum, fracMapped, meanInsert;
    long numReads, numMapped, InsertSum, NonzeroInserts;
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
    }
    BamReader reader;
    if(!reader.Open(inputFilename)) {
        cerr << "Could not open input BAM" << endl;
        return 1;
    }
    const SamHeader header = reader.GetHeader();
    const RefVector references = reader.GetReferenceData();
    BamAlignment rec;
    AFSum = 0.;
    numReads = numMapped = NonzeroInserts = InsertSum = 0;

    while (reader.GetNextAlignment(rec)) {
        numReads++;
        if(rec.IsMapped()) {
            rec.GetTag("AF", AF);
            numMapped++;
            AFSum += AF;
            if(rec.InsertSize != 0 && abs(rec.InsertSize) < MaxInsert) {
                NonzeroInserts++;
                InsertSum += abs(rec.InsertSize);
            }
        }

    }
    meanAF = AFSum / numMapped;
    fracMapped = float(numMapped) / numReads;
    meanInsert = float(InsertSum) / NonzeroInserts;
    printf("Mean Aligned Fraction %s=%f\n", inputFilename.c_str(), meanAF);
    printf("Reads mapped=%i\n", numMapped);
    printf("Reads total=%i\n", numReads);
    printf("Fraction aligned=%f\n", fracMapped);
    printf("Mean insert size=%f\n", meanInsert);
    reader.Close();
    return 0;
}
