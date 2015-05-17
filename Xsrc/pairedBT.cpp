//Hey this is a C++ file that tags ur BAMzzz!!!11!oneone!

#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <numeric>
#include <exception>
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "kseq.hpp"
#include <sys/stat.h>
#include <fcntl.h>

using namespace BamTools;

int main(int argc, char* argv[])
{
    //MaxInsert is the largest insert included in the calculation of mean insert size.
    const long MaxInsert = 2000;
    long recCount = 0;
    long pairCount = 0;
    long kread1 = 0;
    long kread2 = 0;
    long ifLen1, ifLen2, fq1, fq2;
    kseq seq1, seq2;
    std::string inputBam, inputFq1, inputFq2;
    FunctorZlib gzr1, gzr2;

    if(argc == 3) {
        inputBam = argv[1];
        inputFq1 = argv[2];
        inputFq2 = "default";
        throw std::runtime_error("Function not yet written to support single-end tagging.");
    }
    else if(argc > 4) {
        std::cerr << "Too many arguments" << std::endl;
        return 1;
    }
    else if(argc <= 2) {
        printf("Description: adds the appropriate barcode-related tags to a BAM file from a pair of fastq files.\n\n");
        printf("BAM must be name-sorted with only one primary mapping read per template accepted.");
        printf("Usage: \n");
        printf("First positional argument: input BAM.\n");
        printf("Second positional argument: fastq file.\n");
        printf("Third (optional) positional argument: second fastq file.\n");
        printf("If two fastq files are specified, data is treated as paired-end.\n");
        return 1;
    }
    FunctorRead r1, r2;
    inputBam = argv[1];
    inputFq1 = argv[2];
    inputFq2 = argv[3];
    ifLen1 = inputFq1.length();
    fq1 = open(inputFq1.c_str(), O_RDONLY);
    ifLen2 = inputFq2.length();
    fq2 = open(inputFq2.c_str(), O_RDONLY);
    kstream<int, FunctorRead> ks1(fq1, r1);
    kstream<int, FunctorRead> ks2(fq2, r2);
    BamReader reader;
    if(!reader.Open(inputBam)) {
        std::cerr << "Could not open input BAM" << std::endl;
        return 1;
    }
    const SamHeader header = reader.GetHeader();
    const RefVector references = reader.GetReferenceData();
    BamAlignment rec1, rec2, rec;
    kseq* kseq1_ptr = &seq1; 
    kseq* kseq2_ptr = &seq2; 

    while (reader.GetNextAlignment(rec)) {
        if(!rec.IsPrimaryAlignment()){
            continue;
        }
        if(rec.IsFirstMate()) {
            rec1 = rec;
            kread1 = ks1.read(seq1);
            if(kread1 <= 0) {
                throw std::runtime_error("Fewer fastq records than BAMs - these files do not match!");
            }
            continue;
        }
        if(rec.IsSecondMate()) {
            rec2 = rec;
            kread2 = ks2.read(seq2);
            if(kread2 <= 0) {
                throw std::runtime_error("Fewer fastq records than BAMs - these files do not match!");
            }
        }
        if(rec1.Name != rec2.Name) {
            throw std::runtime_error("Read names should match!!");
        }
        if(!kseq1_ptr) {
            throw std::runtime_error("Read 1 fastq object is null!");
        }
        if(!kseq2_ptr) {
            throw std::runtime_error("Read 2 fastq object is null!");
        }
    }
    printf("#Bam=%s\n", inputBam.c_str());
    printf("#Fq1=%s\n", inputFq1.c_str());
    if(inputFq2.compare("default") == 0) {
        printf("#Fq2=%s\n", inputFq2.c_str());
    }
    reader.Close();
    close(fq2);
    close(fq1);
    return 0;
}
