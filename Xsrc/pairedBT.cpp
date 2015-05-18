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
#include <sstream>
#include <sys/stat.h>
#include <fcntl.h>

using namespace BamTools;

std::vector<std::string> splitString(std::string inStr, std::string delimiter) {
    std::vector<std::string> list;
    list.clear();
    size_t pos, dlen;
    dlen = delimiter.length();
    while( (pos = inStr.find(delimiter)) != std::string::npos) {
        list.push_back(inStr.substr(0, pos));
        inStr.erase(0, pos + dlen);
    }
    return list;
}

std::vector<std::string> splitByChar(std::string inStr, char delimiter) {
    std::vector<std::string> tokens;
    std::stringstream ss;
    ss.str(inStr);
    std::string item;
    while (getline(ss, item, delimiter)) {
        tokens.push_back(item);
    }
    return tokens;
}

std::string CoorString(RefVector refs, BamAlignment rec1, BamAlignment rec2) {
    std::stringstream stringbuilder;
    std::string rec1Str, rec2Str;
    if(rec1.IsMapped()) {
        stringbuilder << refs[rec1.RefID].RefName << ":" << rec1.Position;
        rec1Str = stringbuilder.str();
        if(!rec2.IsMapped()) {
            return rec1Str + ",None";
        }
        stringbuilder.clear(); //clear any bits set
        stringbuilder.str(""); //clear any characters set
        stringbuilder << refs[rec2.RefID].RefName << ":" << rec2.Position;
        rec2Str = stringbuilder.str();
        if(rec1.RefID < rec2.RefID) {
            return rec1Str + "," + rec2Str;
        }
        else if(rec2.RefID < rec1.RefID) {
            return rec2Str + "," + rec1Str;
        }
        else if(rec1.Position < rec2.Position) {
            return rec1Str + "," + rec2Str;
        }
        else {
            return rec2Str + "," + rec1Str;
        }
    }
    else{
        rec1Str = "None";
        if(!rec2.IsMapped()) {
            return "None,None";
        }
        stringbuilder << refs[rec2.RefID].RefName << ":" << rec2.Position;
        return stringbuilder.str() + ",None";
    }
}

float SoftClippedFraction(BamAlignment rec) {
    int cigarLen = rec.CigarData.size();
    int recLen = 0;
    int mappedBases = 0;
    CigarOp operation;
    for(int i=0; i++; i < cigarLen) {
        operation = rec.CigarData[i];
        if(operation.Type == 'S'){
            mappedBases += operation.Length;
        }
        recLen += operation.Length;
    }
    return float(mappedBases) / recLen;
}

float AlignedFraction(BamAlignment rec) {
    int cigarLen = rec.CigarData.size();
    int recLen = 0;
    int mappedBases = 0;
    CigarOp operation;
    for(int i=0; i++; i < cigarLen) {
        operation = rec.CigarData[i];
        if(operation.Type == 'M'){
            mappedBases += operation.Length;
        }
        recLen += operation.Length;
    }
    return float(mappedBases) / recLen;
}

int main(int argc, char* argv[])
{
    long recCount = 0;
    long pairCount = 0;
    long kread1 = 0;
    long kread2 = 0;
    long tagsLen;

    // Variables that will be rewritten as command-line argument settable.
    bool markQCFlag = true;
    int head = 4;  // Number of bases to grab off the "head" of R1 and R2 to pad the molecular barcode

    gzFile fq1, fq2;
    kseq seq1, seq2;
    std::string inputBam, outputBam, inputFq1, inputFq2, head1, head2, BarcodeString;
    FunctorZlib gzr1, gzr2;

    if(argc == 3) {
        inputBam = argv[1];
        inputFq1 = argv[2];
        inputFq2 = "default";
        throw std::runtime_error("Function not yet written to support single-end tagging.");
    }
    else if(argc > 5) {
        std::cerr << "Too many arguments" << std::endl;
        return 1;
    }
    else if(argc <= 2) {
        printf("Description: adds the appropriate barcode-related tags to a BAM file from a pair of fastq files.\n\n");
        printf("BAM must be name-sorted with only one primary mapping read per template accepted.\n");
        printf("Usage: \n");
        printf("First positional argument: input BAM.\n");
        printf("Second positional argument: fastq file.\n");
        printf("Third (optional) positional argument: second fastq file.\n");
        printf("If two fastq files are specified, data is treated as paired-end.\n");
        return 1;
    }
    inputBam = argv[1];
    inputFq1 = argv[2];
    inputFq2 = argv[3];
    if(argc == 4) {
        outputBam = splitString(inputBam, ".")[0] + ".tagged.bam";
    }
    else {
        outputBam = argv[4];
    }
    fq1 = gzopen(inputFq1.c_str(), "r");
    fq2 = gzopen(inputFq2.c_str(), "r");
    kstream<gzFile, FunctorZlib> ks1(fq1, gzr1);
    kstream<gzFile, FunctorZlib> ks2(fq2, gzr2);
    BamReader reader;
    if(!reader.Open(inputBam)) {
        std::cerr << "Could not open input BAM" << std::endl;
        return 1;
    }
    std::cout << "Opened bam reader" << std::endl;
    const SamHeader header = reader.GetHeader();
    const RefVector references = reader.GetReferenceData();
    BamAlignment rec1, rec2, rec;
    BamWriter writer;
    if(!writer.Open(outputBam, header, references) ) {
        std::cerr << "ERROR: could not open " + outputBam + " for writing. Abort mission!" << std::endl;
        throw std::runtime_error("Could not open " + outputBam + " for writing.");
    }
    std::cout << "Opened bam writer" << std::endl;

    while (reader.GetNextAlignment(rec)) {
        recCount++;
        std::vector<std::string> TagStrings, TagPair;
        int FM, ND1, ND2, FP;
        float NF1, NF2;
        std::string FAStr1, PVStr1, FAStr2, PVStr2;
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
            pairCount++;
        }
        if(rec1.Name != rec2.Name) {
            printf("Read 1 name: %s. Read 2 name: %s.\n", rec1.Name.c_str(), rec2.Name.c_str());
            throw std::runtime_error("Read names should match!!");
        }
        if(!&seq1) {
            throw std::runtime_error("Read 1 fastq object is null!");
        }
        if(!&seq2) {
            throw std::runtime_error("Read 2 fastq object is null!");
        }
        // Get the information needed for barcode tagging from read 1.
        TagStrings = splitByChar(std::string(seq1.comment), '|');
        tagsLen = TagStrings.size();
        BarcodeString = "";
        for(int i = 0; i < tagsLen; i++) {
            TagPair = splitByChar(TagStrings[i], '=');
            //printf("TagPair[0]: %s\n", TagPair[0].c_str());
            // Get family size
            if(TagPair[0].compare("FM") == 0) {
                FM = atoi(TagPair[1].c_str());
            }
            if(TagPair[0].compare("ND") == 0) {
                ND1 = atoi(TagPair[1].c_str());
                NF1 = float(ND1) / FM;
            }
            // Check for Barcode QC fail
            if(TagPair[0].compare("FP") == 0) {
                if(TagPair[1].find("Fail") != std::string::npos) {
                    FP = 0; // Fails barcode QC
                    if(markQCFlag) {
                        rec1.SetIsFailedQC(true);
                        rec2.SetIsFailedQC(true);
                    }
                }
                else{
                    FP = 1; // Passes barcode QC
                }
            }
            if(TagPair[0].compare("FA") == 0) {
                FAStr1 = TagPair[1];
            }
            if(TagPair[0].compare("PV") == 0) {
                PVStr1 = TagPair[1];
            }
            if(TagPair[0].compare("BS") == 0) {
                BarcodeString = TagPair[1];
            }
        }
        TagStrings = splitByChar(std::string(seq2.comment), '|');
        for(int i = 0; i < tagsLen; i++) {
            TagPair = splitByChar(TagStrings[i], '=');
            if(TagPair[0].compare("ND") == 0) {
                ND2 = atoi(TagPair[1].c_str());
                NF2 = float(ND2) / FM;
            }
            if(TagPair[0].compare("FA") == 0) {
                FAStr2 = TagPair[1];
            }
            if(TagPair[0].compare("PV") == 0) {
                PVStr2 = TagPair[1];
            }
        }

        std::string CoordinateString = CoorString(references, rec1, rec2);
        //Add these tags.
        //Read 1 tags
        rec1.AddTag("FA", "Z", FAStr1);
        rec1.AddTag("PV", "Z", PVStr1);
        rec1.AddTag("ND", "i", ND1);
        rec1.AddTag("NF", "f", NF1);
        rec1.AddTag("AF", "f", AlignedFraction(rec1));
        rec1.AddTag("SF", "f", SoftClippedFraction(rec1));

        //Read 2 tags
        rec2.AddTag("FA", "Z", FAStr2);
        rec2.AddTag("PV", "Z", PVStr2);
        rec2.AddTag("ND", "i", ND2);
        rec2.AddTag("NF", "f", NF2);
        rec2.AddTag("AF", "f", AlignedFraction(rec2));
        rec2.AddTag("SF", "f", SoftClippedFraction(rec2));

        // Shared tags
        rec1.AddTag("CS", "Z", CoordinateString);
        rec2.AddTag("CS", "Z", CoordinateString);
        rec1.AddTag("BS", "Z", BarcodeString);
        rec2.AddTag("BS", "Z", BarcodeString);
        rec1.AddTag("FM", "i", FM);
        rec2.AddTag("FM", "i", FM);
        rec1.AddTag("FP", "i", FP);
        rec2.AddTag("FP", "i", FP);
        
        writer.SaveAlignment(rec1);
        writer.SaveAlignment(rec2);

    }
    writer.Close();
    reader.Close();
    gzclose(fq2);
    gzclose(fq1);
    printf("#Bam=%s\n", inputBam.c_str());
    printf("#Fq1=%s\n", inputFq1.c_str());
    if(inputFq2.compare("default") == 0) {
        printf("#Fq2=%s\n", inputFq2.c_str());
    }
    printf("#Outbam=%s\n", outputBam.c_str());
    printf("#Records processed=%i\n", recCount);
    printf("#Pairs processed=%i\n", pairCount);
    return 0;
}
