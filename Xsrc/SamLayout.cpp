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
#include <google/sparse_hash_map>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include "SamLayout.hpp"

using namespace BamTools;
using google::sparse_hash_map;
using namespace boost::algorithm;

typedef sparse_hash_map<int64_t, int64_t> hash_t;
typedef sparse_hash_map<const char *, int64_t> str_hash_t;

const int MaxInsert = 1000;

// LayoutPos implementation
	int LayoutPos::getReferenceID() {return RefID;}
	void LayoutPos::setReferenceID(int newRef) {RefID = newRef;}

	int LayoutPos::getReadBasePos() {return ReadPosition;} // Returns which base in a read this tuple is.
	void LayoutPos::setReadBasePos(int newReadBasePos) {ReadPosition = newReadBasePos;}

	char LayoutPos::getBase() {return base;}
	void LayoutPos::setBase(char newChar){base = newChar;}

	char LayoutPos::getOperation() {return Operation;}
	void LayoutPos::setOperation(char newChar) {Operation = newChar;}

	int LayoutPos::getPos() {return Position;}
	void LayoutPos::setPos(int newPos) {Position = newPos;}

	int LayoutPos::getQuality() {return quality;}
	void LayoutPos::setQuality(int newQual) {quality = newQual;}
// More of this check could be filled out, but let's just get this compiling.
	void LayoutPos::check() {
        if(Operation == 'D' && base != 'D') {
            std::cerr <<"base for LayoutPos isn't what it should be... Fixing!" << std::endl;
            base = 'D';
        }
    }
    bool LayoutPos::incrementRefPos(){
     	return !(Operation == 'S' || Operation == 'I');
    }
    bool LayoutPos::incrementReadPos(){
     	return !(Operation == 'D');
    }

void LayoutPos::setAttributes(char baseArg, char opArg, int RefIDArg, int posArg, int readPosArg, int qualArg, int strand) {
    Operation = opArg;
    base = baseArg;
    RefID = RefIDArg;
    Position = posArg;
    ReadPosition = readPosArg;
    quality = qualArg;
    strandedness = strand;
}

LayoutPos::LayoutPos(char baseArg, char opArg, int RefIDArg, int posArg, int readPosArg, int qualArg, int strand) {
    setAttributes(baseArg, opArg, RefIDArg, posArg, readPosArg, qualArg, strand);
}

LayoutPos::LayoutPos(){
    //Empty constructor
}

bool LayoutOp::isIns() {return Operation == 'I';}
bool LayoutOp::isDel() {return Operation == 'D';}
bool LayoutOp::isMap() {return Operation == 'M';}
bool LayoutOp::isSoftClipped() {return Operation == 'S';}

char LayoutOp::getOperation() {return Operation;}
void LayoutOp::setOperation(char newOp) {Operation = newOp;}

int LayoutOp::getLength() {return layoutPositions.size();}

int LayoutOp::getReadPos() {return readPos;}
void LayoutOp::setReadPos(int newReadPos) {readPos = newReadPos;}

int LayoutOp::getPos() {return pos;}
void LayoutOp::setPos(int newPos) {pos = newPos;}

int LayoutOp::getRef() {return RefID;}
void LayoutOp::setRef(int newRef) {RefID = newRef;}

int LayoutOp::getStrandedness() {return strandedness;}
void LayoutOp::setStrandedness(int newStrandedness) {strandedness = newStrandedness;}

std::string LayoutOp::getSequence() {return seq;}
void LayoutOp::setSequence(std::string newSeq) {seq = newSeq;}

std::vector<int> LayoutOp::getQuality() {return quality;}
void LayoutOp::setQuality(std::vector<int> quals) {quality = quals;}
void LayoutOp::clearQuality() {quality.clear();}

void LayoutOp::mergeOp(LayoutOp ALO) {
	int newPos = ALO.getPos() > pos ? pos : ALO.getPos();
}
std::vector<LayoutPos> LayoutOp::getLayoutPositions() {
    return layoutPositions;
}
void LayoutOp::replacePos(LayoutPos lPos, int index) {
    assert(lPos.getOperation() == getOperation());
    layoutPositions[index] = lPos;
}
bool LayoutOp::incrementRefPos(){
	return !(Operation == 'S' || Operation == 'I');
}
bool LayoutOp::incrementReadPos(){
	return !(Operation == 'D');
}
void LayoutOp::updatePositions(){
	int cumSum = 0;
	for(int i = 0; i < layoutPositions.size(); i++){
		if(layoutPositions[i].incrementRefPos()){cumSum++;}
		layoutPositions[i].setPos(layoutPositions[i].incrementRefPos() ? pos + cumSum : -1);
	}
}
void LayoutOp::updateReadPositions() {
	int cumSum = 0;
	for(int i = 0; i < layoutPositions.size(); i++){
		cumSum = (layoutPositions[i].getOperation() != 'D') ? cumSum + 1 : cumSum; //Increments if the operation is anything but a deletion.
		layoutPositions[i].setReadBasePos(layoutPositions[i].incrementReadPos() ? readPos + cumSum : -1); // Sets read position to -1 for a deletion or updates the position within the read.
	}
}
void LayoutOp::update(){
	updatePositions();
	updateReadPositions();
}

LayoutOp::LayoutOp(){
    //Empty constructor
}

void LayoutOp::setAttributes(std::string cigarSeq, std::vector<int> qualitySlice, char cigarOperation, int RefIDArg, int startPos, int readStart, int strand)
{
    // Convert the quality string to ints.
    //qualities = std::vector<int>;
    quality = qualitySlice;
    Operation = cigarOperation;
    pos = startPos;
    readPos = readStart;
    seq = cigarSeq;
    RefID = RefIDArg;
    strandedness = strand;
    //Build LayoutPos object vector!
    std::vector<LayoutPos> layoutPosVector;
    // Create LayoutPos objects for each.
    LayoutPos ALP = LayoutPos();
    for(int k = 0; k < qualitySlice.size(); k++) {
        switch (Operation){
        case 'I':
            ALP.setAttributes(seq.at(k), Operation, RefID, -1, readPos + k, quality[k], strandedness);
            break;
        case 'S':
            ALP.setAttributes(seq.at(k), Operation, RefID, -1, readPos + k, quality[k], strandedness);
            break;
        case 'D':
            ALP.setAttributes('D', Operation, RefID, pos + k, -1, -1, strandedness);
            break;
        case 'M':
            ALP.setAttributes(seq.at(k), Operation, RefID, pos + k, readPos + k, quality[k], strandedness);
            break;
        default:
            throw std::runtime_error("Sorry, unsupported cigar character. Email me and I'll change this.");
        }
        layoutPosVector.push_back(ALP);
    }
    layoutPositions = layoutPosVector;
}

LayoutOp::LayoutOp(std::string cigarSeq, std::vector<int> qualitySlice, char cigarOperation, int RefIDArg, int startPos, int readStart,
		int strand) {
    setAttributes(cigarSeq, qualitySlice, cigarOperation, RefIDArg, startPos, readStart, strand);
}

void LayoutOp::setAttributes(std::string cigarSeq, std::string cigarQual, char cigarOperation, int RefIDArg, int startPos, int readStart, int strand) {
    // Convert the quality string to ints.
    //qualities = std::vector<int>;
    for(int k = 0;k < cigarQual.size(); k++) {
        quality.push_back((int) cigarQual.at(k) - 33); // Subtract 33 to convert from ASCII value to phred score.
    }
    strandedness = strand;
    Operation = cigarOperation;
    pos = startPos;
    readPos = readStart;
    seq = cigarSeq;
    RefID = RefIDArg;
    //Build LayoutPos object vector!
    std::vector<LayoutPos> layoutPosVector;
    // Create LayoutPos objects for each.
    LayoutPos ALP = LayoutPos();
    for(int k = 0; k < quality.size(); k++) {
        switch (Operation){
        case 'I':
            ALP.setAttributes(seq.at(k), Operation, RefID, -1, readPos + k, quality[k], strandedness);
            break;
        case 'S':
            ALP.setAttributes(seq.at(k), Operation, RefID, -1, readPos + k, quality[k], strandedness);
            break;
        case 'D':
            ALP.setAttributes('D', Operation, RefID, pos + k, -1, -1, strandedness);
            break;
        case 'M':
            ALP.setAttributes(seq.at(k), Operation, RefID, pos + k, readPos + k, quality[k], strandedness);
            break;
        default:
            throw std::runtime_error("Sorry, unsupported cigar character. Email me and I'll change this.");
        }
        layoutPosVector.push_back(ALP);
    }
    layoutPositions = layoutPosVector;
}

// Warning... This constructor should only be used for reads without barcodes.
LayoutOp::LayoutOp(std::string cigarSeq, std::string cigarQual, char cigarOperation, int RefIDArg, int startPos, int readStart,
		int strand)
{
    setAttributes(cigarSeq, cigarQual, cigarOperation, RefIDArg, startPos, readStart, strand);
}


std::vector<int> PhredVectorFromString(std::string PVString){
	std::vector<int> returnVec;
	std::vector<std::string> stringVec;
	boost::split(stringVec, PVString, is_any_of(","));
	for(int i = 0; i < returnVec.size(); i++){
		returnVec.push_back(std::atoi(stringVec[i].c_str()));
	}
	return returnVec;
}

std::vector<std::string> IntToStringVec(std::vector<int> TemplateVec){
	std::vector<std::string> returnVec;
	for(int i = 0; i < TemplateVec.size(); i++){
		returnVec.push_back(boost::lexical_cast<std::string>(TemplateVec[i]));
	}
	return returnVec;
}

std::string PhredStringFromVector(std::vector<int> PVArray){
	// Gets the phred value (PV) string from such a vector. (114,133,160,...)
	return boost::algorithm::join(IntToStringVec(PVArray), ",");
}


std::vector<LayoutOp> GetLayoutOps(BamAlignment rec) {
    int cigarLen;
    int StartPosition = rec.Position; // Needed to keep read indices in line with reference.
    int cumCigarSum = 0;
    // These next variables are ones I hope to remove once I figure out what I need
    // to initialize LayoutOp.
    std::vector<LayoutOp> operations;
    LayoutOp TmpOp;
    CigarOp TmpCigarOp;
    cigarLen = rec.CigarData.size();
    std::string PVString;
    for(int i = 0; i < cigarLen; i++) {
        TmpCigarOp = rec.CigarData[i];
        if(!rec.HasTag("PV")) {
        	// Gets quality scores from ASCII phred scores rather than PV numbers.
            operations.push_back(
                LayoutOp(rec.QueryBases.substr(cumCigarSum, cumCigarSum + TmpCigarOp.Length), // Read sequence
                                  rec.Qualities.substr(cumCigarSum, cumCigarSum + TmpCigarOp.Length), // Read qualities
                                  TmpCigarOp.Type, // Cigar operation
                                  rec.RefID, // Reference ID
                                  StartPosition + cumCigarSum, // Genomic coords for start of cigar
                                  cumCigarSum, // Read coords for start of cigar
								  rec.IsReverseStrand() ? -1 : 1
                                  ));
            cumCigarSum += TmpCigarOp.Length;
        }
        else {
            std::string PVString;
            rec.GetTag("PV", PVString);
            std::vector<int> PhredSubVector;
            std::vector<int> PhredVector = PhredVectorFromString(PVString);
            for(int k = cumCigarSum; k < cigarLen + cumCigarSum; k++)  {
                PhredSubVector.push_back(PhredVector[k]);
            }
            LayoutOp tmpOp = LayoutOp(rec.QueryBases.substr(cumCigarSum, cumCigarSum + TmpCigarOp.Length), //Read sequence
                    PhredSubVector, TmpCigarOp.Type, rec.RefID, StartPosition + cumCigarSum,
                    cumCigarSum, rec.IsReverseStrand() ? -1 : 1);
            operations.push_back(tmpOp);
            cumCigarSum += TmpCigarOp.Length;
        }
    }
    return operations;
}

class AlignmentLayout {
    private:
        int length; // Length of the read.
        int pos; // 0-based genomic start position of the cigar operation.
        int RefID;
        std::string seq;
        std::string Name; // Read name
        std::string TagData; // Tags for the read, used for finally recreating a Bam record.
        std::vector<int> quality; // Quality scores as integers. PV tags if available, normal phred-encoded otherwise.
        std::vector<LayoutOp> operations;
        int strandedness; // -1 for reverse, 1 for forward, 0 for a merged read pair.
        bool pairMerged;

    public:
        AlignmentLayout(BamAlignment rec);

        int getRef() {return RefID;}
        void setRef(int reference) {RefID = reference;}
        int getPos() {return pos;}
        void setPos(int position) {pos = position;}

        bool isMerged() {return pairMerged;}
        void setIsMerged(bool newMergedValue) {pairMerged = newMergedValue;} // Used to denote that read 1 and read 2 in a pair have been merged into a single read.

        int getLen() {return length;}
        void updateLen() {length = seq.size();}
        std::string getSeq() {return seq;}
        void setSeq(std::string newSeq) {
        	seq = newSeq;
        }
        void updateOpCoords() {
        	pos = operations[0].getPos();
        	int cumPos = 0;
        	for(int i = 0; i < operations.size(); i++){
        		operations[i].setPos(pos + cumPos);
        		cumPos += operations[i].getLength();
        		operations[i].setReadPos(cumPos);
        		operations[i].update();
        	}
        }
        void addOperation(LayoutOp ALO) {
        	if(RefID != ALO.getRef()){
        		std::cerr << "Can not add LayoutOp to read, as it is mapped to a different reference." << std::endl;
        		return;
        	}
        	length += ALO.incrementReadPos() ? ALO.getLength(): 0;
        	if(ALO.getPos() -1 == length + pos){
        		//LayoutOp is just after current vector of operations.
        		operations.push_back(ALO);
        		seq += ALO.getSequence();
        		quality.reserve(quality.size() + ALO.getQuality().size());
        		quality.insert(quality.end(), ALO.getQuality().begin(), ALO.getQuality().end()); // Appends ALO's qualities to the appropriate vector
        		// ALO.clearQuality(); not sure if I need this line.
        		if(ALO.incrementRefPos()) {pos += ALO.getLength();}
        		updateOpCoords();
        	}
        	else if (ALO.getPos() + 1 == pos){
        		//LayoutOp is just before the current vector of operations.
        		length += ALO.incrementReadPos() ? ALO.getLength(): 0;
        		operations.insert(operations.begin(), ALO);
        		seq = ALO.getSequence() + seq;
        		if(ALO.incrementRefPos()) {pos -= ALO.getLength();}
        		updateOpCoords();
        	}
        	else {
        		std::cerr << "I can't figure out what you're doing. ALO doesn't border the layout to which it is being added." << std::endl;
        	}
        	for(int i = 0; i < operations.size() - 1; i++){
        		if(operations[i].getOperation() == operations[i + 1].getOperation()){
        			operations[i].mergeOp(operations[i + 1]);
        			operations.erase(operations.begin() + i + 1);
        		}
        	}
        }
        std::vector<LayoutOp> getOps() {return operations;}
        void setOpts(std::vector<LayoutOp> newVec) {operations = newVec;}
        int getStrandedness() {return strandedness;}
        void setStrandedness(int newStrandedness) {strandedness = newStrandedness;}
        std::vector<LayoutPos> getLayoutPositions(){
            std::vector<LayoutPos> positions;
            std::vector<LayoutPos> tmpPositions;
            for(int i = 0; i < operations.size(); i++){
                tmpPositions = operations[i].getLayoutPositions();
                for(int ii = 0; ii < tmpPositions.size(); ii++) {
                    positions.push_back(tmpPositions[i]);
                }
            }
            return positions;
        }
        int getAlignedLen() {
            int len = 0;
            for(int i = 0; i < operations.size(); i++){
                if(operations[i].getOperation() == 'M') {
                    len += operations[i].getLength();
                }
            }
            return len;
        }

        BamAlignment toBam(){
            BamAlignment returnBam;
            returnBam.Name = Name;
            returnBam.QueryBases = seq;
            returnBam.RefID = RefID;
            returnBam.Position = pos;
            returnBam.SetIsReverseStrand(strandedness);
            if(pairMerged){
            	returnBam.AddTag("MP", "A", 'T');
            }
            returnBam.AddTag("PV", "Z", PhredStringFromVector(quality));
            std::vector<std::string> Tags;
            boost::split(Tags, returnBam.TagData,is_any_of("\t"));
            for(int i = 0; i < Tags.size(); i++){
            	std::vector<std::string> TagElements;
            	boost::split(TagElements, Tags[i], is_any_of(":"));
            	if(TagElements[0].compare("PV") == 0 || TagElements[0].compare("MP") == 0 ){
            		continue;
            	}
            	switch (TagElements[1].at(0)){
				    case 'A':
            		    returnBam.AddTag(TagElements[0], "A", (char) TagElements[2].at(0));
            		    break;
				    case 'Z':
					    returnBam.AddTag(TagElements[0], "Z", TagElements[2]);
					    break;
				    case 'i':
				    	returnBam.AddTag(TagElements[0], "i", atoi(TagElements[2].c_str()));
				    	break;
				    case 'f':
				    	returnBam.AddTag(TagElements[0], "f", atof(TagElements[2].c_str()));
				    	break;
				    default:
					    continue; // Don't sweat it.
            	}
            }
            return returnBam;
        }
};

AlignmentLayout::AlignmentLayout(BamAlignment rec) {
	if(rec.HasTag("MP")){
		char MPChar;
		rec.GetTag("MP", MPChar);
		pairMerged = (MPChar == 'T') ? true : false;
	}
	else {
		pairMerged = false;
	}
    Name = rec.Name;
    length = rec.Length;
    seq = rec.QueryBases;
    RefID = rec.RefID;
    pos = rec.Position;
    TagData = rec.TagData;
    std::string PVString;
    rec.GetTag("PV", PVString);
    quality = PhredVectorFromString(PVString);
    strandedness = rec.IsReverseStrand() ? -1 : 1; // Replace that if/then chain with ternary
    operations = GetLayoutOps(rec);
}


std::vector<BamAlignment> MergePairedAlignments(BamAlignment rec1, BamAlignment rec2) {
    // If the two overlap, make one new read that's all perfect for it. :)
    std::vector<BamAlignment> returnList;
    int startpos;
    if(rec1.RefID != rec2.RefID || abs(rec1.Position - rec2.Position) > MaxInsert || rec1.IsReverseStrand() == rec2.IsReverseStrand()){
        // Reads aren't anywhere near each other or are on the same strand. Return original alignments after adding an attempted merge tag.
        rec1.AddTag("MP", "A", 'F');
        rec2.AddTag("MP", "A", 'F');
        returnList.push_back(rec1);
        returnList.push_back(rec2);
        return returnList;
    }
    if(rec1.Position <= rec2.Position) {
        startpos = rec1.Position;
    }
    else {
        startpos = rec2.Position;
    }
    assert(false); // This function hasn't been finished. Make it fail assertion in case it ends up being used.
    return returnList;
}

int main(){
    return 0;
}
