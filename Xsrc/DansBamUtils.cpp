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
#include <assert.h>

using namespace BamTools;

const int MaxInsert = 1000;
class AlignmentLayoutPos {
    /** Contains a nucleotide, 
     */
    private:
        char Operation; //cigar operation for this base.
        int RefID;
        int Position; // Set to -1 for an insertion or soft-clipping
        int ReadPosition; // Set to -1 for a deletion
        std::string base; // Set to "D" for a deletion.
        int quality; // Set to -1 for a deletion.

    public:
        AlignmentLayoutPos(std::string, char, int, int, int, int); // Base, ref id, Position (ref), Position (read), Quality
        AlignmentLayoutPos();
        void setAttributes(std::string, char, int, int, int, int);
        int getReferenceID() {return RefID;}
        int getReadBasePos() {return ReadPosition;} // Returns which base in a read this tuple is.
        std::string getBase() {return base;}
        int getPos() {return Position;}
        int getQuality() {return quality;}
        // More of this check could be filled out, but let's get this compiling.
        void check() {
            if(Operation == 'D' && base != "D") {
                throw std::runtime_error("base for AlignmentLayoutPos isn't what it should be...");
            }
            if(base.size() != 1){
                throw std::runtime_error("base is the wrong size for AlignmentLayoutPos. It should be precisely one character.");
            }
        }
};

void AlignmentLayoutPos::setAttributes(std::string baseArg, char opArg, int RefIDArg, int posArg, int readPosArg, int qualArg) {
    Operation = opArg;
    base = baseArg;
    RefID = RefIDArg;
    Position = posArg;
    ReadPosition = readPosArg;
    quality = qualArg;
}

AlignmentLayoutPos::AlignmentLayoutPos(std::string baseArg, char opArg, int RefIDArg, int posArg, int readPosArg, int qualArg) {
    Operation = opArg;
    base = baseArg;
    RefID = RefIDArg;
    Position = posArg;
    ReadPosition = readPosArg;
    quality = qualArg;
}

AlignmentLayoutPos::AlignmentLayoutPos() {
    //Empty constructor
}

class AlignmentLayoutOp {
    /** 
     * Contains the information needed from a CigarOp object and a read to
     * make a list of AlignmentLayoutPos*/
    private:
        char Operation; // Character definition of cigar operation.
        int length; // Length of cigar operation
        int pos; // 0-based genomic start position of the cigar operation.
        int readPos; // 0-based read start position of the cigar operation.
        int RefID; // Reference ID.
        std::string seq;
        std::vector<int> quality; // Store the information as integers. For this purpose, use the PV tags if available.
        std::vector<AlignmentLayoutPos> layoutPositions;

    public:
        AlignmentLayoutOp(std::string, int, int, int, std::string, char); //Constructor
        AlignmentLayoutOp();
        void setAttributes(std::string, int, int, int, std::string, char); // For setting values if declared before constructed.
        bool isIns() {return Operation == 'I';}
        bool isDel() {return Operation == 'D';}
        bool isMap() {return Operation == 'M';}
        bool isSoftClipped() {return Operation == 'S';}
        char getOperation() {return Operation;}
        int getLength() {return length;}
        int getPos() {return pos;}
        std::vector<int> getQuality() {return quality;}
        std::string getSequence() {return seq;}
        std::vector<AlignmentLayoutPos> getLayoutPositions() {
            //assert length == layoutPositions.size(); // If layout positions have changed, make sure things are still sane.
            return layoutPositions;
        }
        void replacePos(AlignmentLayoutPos lPos, int index) {
            layoutPositions[index] = lPos;
        }
};

//Constructor meant to be empty
AlignmentLayoutOp::AlignmentLayoutOp() {
    //empty on purpose.
}

void AlignmentLayoutOp::setAttributes(std::string cigarSeq, int RefIDArg, int startPos, int readStart, std::string cigarQual,
                                 char cigarOperation) {
    // Convert the quality string to ints.
    //qualities = std::vector<int>;
    length = cigarQual.size();
    for(int k = 0;k < length; k++) {
        quality.push_back((int) cigarQual.at(k));
        //qualities.push_back((int) cigarQual.at(k));
    }
    //quality = qualities;
    Operation = cigarOperation;
    pos = startPos;
    readPos = readStart;
    seq = cigarSeq;
    RefID = RefIDArg;
    //Build AlignmentLayoutPos object vector!
    std::vector<AlignmentLayoutPos> layoutPosVector;
    // Create AlignmentLayoutPos objects for each.
    AlignmentLayoutPos ALP = AlignmentLayoutPos();
    for(int k = 0; k < length; k++) {
        switch (Operation){
        case 'I':
            ALP.setAttributes(seq.substr(k, k + 1), Operation, RefID, -1, readPos + k, quality[k]);
            break;
        case 'S':
            ALP.setAttributes(seq.substr(k, k + 1), Operation, RefID, -1, readPos + k, quality[k]);
            break;
        case 'D':
            ALP.setAttributes("D", Operation, RefID, pos + k, -1, -1);
            break;
        case 'M':
            ALP.setAttributes(seq.substr(k, k + 1), Operation, RefID, pos + k, readPos + k, quality[k]);
            break;
        default:
            throw std::runtime_error("Sorry, unsupported cigar character. Email me and I'll change this.");
        }
        layoutPosVector.push_back(ALP);
    }
    layoutPositions = layoutPosVector;
}

/**
 * Function to turn a bam record into a list of AlignmentLayoutOp objects.
std::vector<AlignmentLayoutOp> GetAlignmentLayoutOps(BamAlignment rec) {
    int cigarLen;
    int StartPosition = rec.Position; // Needed to keep read indices in line with reference.
} */

// Warning... This constructor should only be used for reads without barcodes.
AlignmentLayoutOp::AlignmentLayoutOp(std::string cigarSeq, int RefIDArg, int startPos, int readStart, std::string cigarQual,
                                     char cigarOperation) {
    // Convert the quality string to ints.
    //qualities = std::vector<int>;
    length = cigarQual.size();
    for(int k = 0;k < length; k++) {
        quality.push_back((int) cigarQual.at(k));
        //qualities.push_back((int) cigarQual.at(k));
    }
    //quality = qualities;
    Operation = cigarOperation;
    pos = startPos;
    readPos = readStart;
    seq = cigarSeq;
    RefID = RefIDArg;
    //Build AlignmentLayoutPos object vector!
    std::vector<AlignmentLayoutPos> layoutPosVector;
    // Create AlignmentLayoutPos objects for each.
    AlignmentLayoutPos ALP = AlignmentLayoutPos();
    for(int k = 0; k < length; k++) {
        switch (Operation){
        case 'I':
            ALP.setAttributes(seq.substr(k, k + 1), Operation, RefID, -1, readPos + k, quality[k]);
            break;
        case 'S':
            ALP.setAttributes(seq.substr(k, k + 1), Operation, RefID, -1, readPos + k, quality[k]);
            break;
        case 'D':
            ALP.setAttributes("D", Operation, RefID, pos + k, -1, -1);
            break;
        case 'M':
            ALP.setAttributes(seq.substr(k, k + 1), Operation, RefID, pos + k, readPos + k, quality[k]);
            break;
        default:
            throw std::runtime_error("Sorry, unsupported cigar character. Email me and I'll change this.");
        }
        layoutPosVector.push_back(ALP);
    }
    layoutPositions = layoutPosVector;
}


std::vector<AlignmentLayoutOp> GetAlignmentLayoutOps(BamAlignment rec) {
    int cigarLen;
    int StartPosition = rec.Position; // Needed to keep read indices in line with reference.
    int cumCigarSum = 0;
    // These next variables are ones I hope to remove once I figure out what I need
    // to initialize AlignmentLayoutOp.
    std::vector<AlignmentLayoutOp> operations;
    AlignmentLayoutOp TmpOp;
    CigarOp TmpCigarOp;
    cigarLen = rec.CigarData.size();
    std::string PVString;
    for(int i = 0; i < cigarLen; i++) {
        TmpCigarOp = rec.CigarData[i];
        if(!rec.HasTag("PV")) {
            operations.push_back(
                AlignmentLayoutOp(rec.QueryBases.substr(cumCigarSum, cumCigarSum + TmpCigarOp.Length), // Read sequence
                                  rec.RefID, // Reference ID
                                  StartPosition + cumCigarSum, // Genomic coords for start of cigar
                                  cumCigarSum, // Read coords for start of cigar
                                  rec.Qualities.substr(cumCigarSum, cumCigarSum + TmpCigarOp.Length),
                                  TmpCigarOp.Type));
        }
        else {
            int PVArray[cigarLen];
            for(int k = 0; k < cigarLen; k++) {
                
            }
        }
    }
}

std::vector<BamAlignment> MergePairedAlignments(BamAlignment rec1, BamAlignment rec2) {
    // If the two overlap, make one new read that's all perfect for it. :)
    std::vector<BamAlignment> returnList;
    int startpos;
    if(rec1.RefID != rec2.RefID || abs(rec1.Position - rec2.Position) > MaxInsert || rec1.IsReverseStrand() == rec2.IsReverseStrand()){
        // Reads aren't anywhere near each other or are on the same strand. Return original alignments after adding an attempted merge tag.
        //rec1.AddTag("MP", "A", "F");
        //rec2.AddTag("MP", "A", "F");
        returnList.push_back(rec1);
        returnList.push_back(rec2);
        return returnList;
    }
    if(rec1.Position <= rec2.Position) {
        startpos = rec1.Position;
    }
}

int main(){
    return 0;
}
