//Hey this is a C++ file!

#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <numeric>
#include "api/BamReader.h"
#include "api/BamWriter.h"

using namespace BamTools;

int main(int argc, char* argv[])
{
    //MaxInsert is the largest insert included in the calculation of mean insert size.
    const long MaxInsert = 2000;
    std::string inputFilename;
    std::vector<float> AFs;
    std::vector<float> InsertVec;
    float meanAF, AF, AFSum, fracMapped, meanInsert;
    long numReads, numMapped, InsertSum, NonzeroInserts;
    if(argc == 2) {
        inputFilename = argv[1];
    }
    else if(argc > 2) {
        std::cerr << "Too many arguments" << std::endl;
        return 1;
    }
    else {
        printf("Description: Calculates the mean aligned fraction for a BAM.\n\n");
        printf("Usage: \n");
        printf("First positional argument: input BAM\n");
        printf("That is all.\n");
        return 1;
    }
    BamReader reader;
    if(!reader.Open(inputFilename)) {
        std::cerr << "Could not open input BAM" << std::endl;
        return 1;
    }
    const SamHeader header = reader.GetHeader();
    const RefVector references = reader.GetReferenceData();
    BamAlignment rec;
    numReads = NonzeroInserts = InsertSum = 0;

    while (reader.GetNextAlignment(rec)) {
        numReads++;
        if(rec.IsMapped()) {
            rec.GetTag("AF", AF);
            AFs.push_back(AF);
            if(rec.InsertSize != 0 && abs(rec.InsertSize) < MaxInsert) {
                NonzeroInserts++;
                InsertVec.push_back(abs(rec.InsertSize));
            }
        }

    }
    NonzeroInserts = long(InsertVec.size());
    AFSum = std::accumulate(AFs.begin(), AFs.end(), 0.0);
    InsertSum = std::accumulate(InsertVec.begin(), InsertVec.end(), 0);
    numMapped = long(AFs.size());
    meanAF = AFSum / numMapped;
    fracMapped = float(numMapped) / numReads;
    meanInsert = float(InsertSum) / NonzeroInserts;
    printf("#Filename=%s\n", inputFilename.c_str());
    printf("Mean Aligned Fraction=%f\n", meanAF);
    printf("Reads mapped=%i\n", numMapped);
    printf("Reads total=%i\n", numReads);
    printf("Fraction aligned=%f\n", fracMapped);
    printf("Mean insert size=%f\n", meanInsert);
    reader.Close();
    return 0;
}
