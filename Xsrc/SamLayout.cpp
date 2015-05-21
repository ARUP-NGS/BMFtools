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

std::string LayoutPos::__str__(){
    /*
     * String representation of a Layout Position.
     * RefID, position, read position, base call, cigar operation,
     * and base quality.
     */
    std::stringstream ss;
    ss.str("");
    ss.clear();
    ss << "RefID:" << RefID << ";Position:"<< Position;;
    ss << ";readPosition:" << ReadPosition << ";base:" << base;
    ss << ";CigarOperation:" << Operation << ";Quality:" << quality;
    return ss.str();
}

void LayoutPos::setAttributes(char baseArg, char opArg, int RefIDArg, int posArg,
                              int readPosArg, int qualArg, int strand, int FA) {
    Operation = opArg;
    base = baseArg;
    RefID = RefIDArg;
    Position = posArg;
    ReadPosition = readPosArg;
    quality = qualArg;
    strandedness = strand;
    agreement = FA;
}

LayoutPos::LayoutPos(char baseArg, char opArg, int RefIDArg, int posArg, int readPosArg, int qualArg, int strand, int FA) {
    setAttributes(baseArg, opArg, RefIDArg, posArg, readPosArg, qualArg, strand, FA);
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

std::vector<std::string> LayoutOp::getPositionStrs() {
    std::vector<std::string> returnVec;
    for(int i = 0; i < getLength(); i++){
        returnVec.push_back(layoutPositions[i].__str__());
    }
    return returnVec;
}

std::vector<char> LayoutOp::getBaseCalls() {
    std::vector<char> returnVec;
    for(int i = 0; i < layoutPositions.size(); i++) {
        returnVec.push_back(layoutPositions[i].getBase());
    }
    return returnVec;
}


std::vector<int> LayoutOp::getReadPositions() {
    std::vector<int> returnVec;
    for(int i = 0; i < layoutPositions.size(); i++) {
        returnVec.push_back(layoutPositions[i].getReadBasePos());
    }
    return returnVec;
}


std::vector<int> LayoutOp::getPositions() {
    std::vector<int> returnVec;
    for(int i = 0; i < layoutPositions.size(); i++) {
        returnVec.push_back(layoutPositions[i].getPos());
    }
    return returnVec;
}

/*
 * Returns a string representation of the LayoutOp object.
 * Should be useful for debugging and knowing what we're working with.
 */
std::string LayoutOp::__str__(){
    std::stringstream ss;
    ss.str("");
    ss.clear();
    ss << "RefID:" << RefID << ";Position:"<< pos;
    ss << ";readPosition:" << readPos << ";Sequence:" << seq;
    ss << ";CigarOperation:" << Operation << ";Quality:" << PhredStringFromVector(quality);
    return ss.str();
}

LayoutOp::LayoutOp(){
    //Empty constructor
}

void LayoutOp::setAttributes(std::string cigarSeq, std::vector<int> qualitySlice, std::vector<int> agreementSlice,
        char cigarOperation, int RefIDArg, int startPos, int readStart, int strand)
{
    // Convert the quality string to ints.
    //qualities = std::vector<int>;
    quality = qualitySlice;
    agreement = agreementSlice;
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
            ALP.setAttributes(seq.at(k), Operation, RefID, -1, readPos + k, quality[k], strandedness, agreement[k]);
            break;
        case 'S':
            ALP.setAttributes(seq.at(k), Operation, RefID, -1, readPos + k, quality[k], strandedness, agreement[k]);
            break;
        case 'D':
            ALP.setAttributes('D', Operation, RefID, pos + k, -1, -1, strandedness, agreement[k]);
            break;
        case 'M':
            ALP.setAttributes(seq.at(k), Operation, RefID, pos + k, readPos + k, quality[k], strandedness, agreement[k]);
            break;
        default:
            throw std::runtime_error("Sorry, unsupported cigar character. Email me and I'll change this.");
        }
        layoutPosVector.push_back(ALP);
    }
    layoutPositions = layoutPosVector;
}

LayoutOp::LayoutOp(std::string cigarSeq, std::vector<int> qualitySlice, std::vector<int> agreementSlice,
        char cigarOperation, int RefIDArg, int startPos, int readStart, int strand) {
    setAttributes(cigarSeq, qualitySlice, agreementSlice, cigarOperation, RefIDArg, startPos, readStart, strand);
}

//This constructor should only be used for reads that have not been barcoded!
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
            ALP.setAttributes(seq.at(k), Operation, RefID, -1, readPos + k, quality[k], strandedness, 1);
            break;
        case 'S':
            ALP.setAttributes(seq.at(k), Operation, RefID, -1, readPos + k, quality[k], strandedness, 1);
            break;
        case 'D':
            ALP.setAttributes('D', Operation, RefID, pos + k, -1, -1, strandedness, 1);
            break;
        case 'M':
            ALP.setAttributes(seq.at(k), Operation, RefID, pos + k, readPos + k, quality[k], strandedness, 1);
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
    boost::split(stringVec, PVString, boost::is_any_of(","));
    for(int k = 0; k < stringVec.size(); k++){
    	std::cerr << "String item here is " << stringVec[k] << std::endl;
    }
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
	std::cerr << "Starting GetLayoutOps." << std::endl;
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
	std::cerr << "About to start creating cigar operations.." << std::endl;
    for(int i = 0; i < cigarLen; i++) {
    	std::cerr << "Now working with cigar operation " << i + 1 << " for this record." << std::endl;
        TmpCigarOp = rec.CigarData[i];
        if(!rec.HasTag("PV")) {
        	std::cerr << "This record has no PV tag" << std::endl;
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
        	std::cerr << "This record has a PV tag!" << std::endl;
            std::string PVString;
            rec.GetTag("PV", PVString);
            std::string FAString;
            rec.GetTag("FA", FAString);
            std::vector<int> PhredSubVector;
            std::vector<int> AgreementSubVector;
        	std::cerr << "About to get Phred/Agreement Vectors" << std::endl;
            std::vector<int> PhredVector = PhredVectorFromString(PVString);
            std::vector<int> AgreementVector = PhredVectorFromString(FAString);
        	std::cerr << "About to get my slice from  Phred/Agreement Vectors" << std::endl;
        	for(int i = 0; i < PhredVector.size(); i++){
        		std::cerr << "Phred value at position " << i + 1 << " is " << PhredVector[i] << std::endl;
        	}
            for(int k = cumCigarSum; k < cigarLen + cumCigarSum; k++)  {
            	std::cerr << "Here is the PhredVector string: " << PVString << std::endl;
                PhredSubVector.push_back(PhredVector[k]);
                AgreementSubVector.push_back(AgreementVector[k]);
            }
        	std::cerr << "About to initialize the new LayoutOp object" << std::endl;
            LayoutOp tmpOp = LayoutOp(rec.QueryBases.substr(cumCigarSum, cumCigarSum + TmpCigarOp.Length), //Read sequence
                    PhredSubVector, AgreementSubVector, TmpCigarOp.Type, rec.RefID, StartPosition + cumCigarSum,
                    cumCigarSum, rec.IsReverseStrand() ? -1 : 1);
        	std::cerr << "Successfully initialized the new LayoutOp object" << std::endl;
            operations.push_back(tmpOp);
            cumCigarSum += TmpCigarOp.Length;
        }
    }
    return operations;
}

int AlnLayout::getRef() {return RefID;}
void AlnLayout::setRef(int reference) {RefID = reference;}
int AlnLayout::getPos() {return pos;}
void AlnLayout::setPos(int position) {pos = position;}

bool AlnLayout::isMerged() {return pairMerged;}
void AlnLayout::setIsMerged(bool newMergedValue) {pairMerged = newMergedValue;} // Used to denote that read 1 and read 2 in a pair have been merged into a single read.

int AlnLayout::getLen() {return length;}
void AlnLayout::updateLen() {length = seq.size();}
std::string AlnLayout::getSeq() {return seq;}
void AlnLayout::setSeq(std::string newSeq) {
    seq = newSeq;
}

void AlnLayout::updateOpCoords() {
    pos = operations[0].getPos();
    int cumPos = 0;
    for(int i = 0; i < operations.size(); i++){
        operations[i].setPos(pos + cumPos);
        cumPos += operations[i].getLength();
        operations[i].setReadPos(cumPos);
        operations[i].update();
    }
}

void AlnLayout::addOperation(LayoutOp ALO) {
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
std::vector<LayoutOp> AlnLayout::getOps() {return operations;}
void AlnLayout::setOpts(std::vector<LayoutOp> newVec) {operations = newVec;}
int AlnLayout::getStrandedness() {return strandedness;}
void AlnLayout::setStrandedness(int newStrandedness) {strandedness = newStrandedness;}
/*
 * Returns a std::vector of all the LayoutPos objects
 * contained in the AlnLayout.
 */
std::vector<LayoutPos> AlnLayout::getLayoutPositions(){
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

std::string CigarDataToString(std::vector<CigarOp> cigarVec){
    std::string returnStr = "";
    std::stringstream ss;
    for(int i = 0; i < cigarVec.size(); i++) {
        ss << cigarVec[i].Length << cigarVec[i].Type;
        returnStr += ss.str();
        ss.str(std::string());
        ss.clear();
    }
    return returnStr;
}

/*
 * Turns a BAM record into a string.
 */
std::string BamToString(BamAlignment rec, RefVector vec){
    std::stringstream ss;
    std::string returnStr;
    std::string TagData = getBamTagString(rec);
    ss << rec.Name + "\t" << rec.AlignmentFlag << "\t" + vec[rec.RefID].RefName
    << "\t" << rec.Position << "\t" << rec.MapQuality << "\t"
    << CigarDataToString(rec.CigarData) << "\t" << vec[rec.MateRefID].RefName
    << "\t" << rec.MatePosition << "\t" << rec.InsertSize << "\t"
    << rec.QueryBases << "\t" << rec.Qualities << "\t" << TagData;
    // Gets first 11 fields and the entire tag list string.
    return ss.str();
}

/*
 * Turns a TagData object into something human readable.
 */

/*
 * Returns a vector of all of the reference bases
 * contained in the LayoutPos objects in the order
 * in which they are found.
 * a -1 means that it is either soft-clipped or there
 * is an insertion.
 */

std::vector<int> AlnLayout::getPositions() {
    std::vector<int> returnVec;
    std::vector<LayoutPos> layouts = getLayoutPositions();
    for(int i = 0; i < layouts.size(); i++){
        returnVec.push_back(layouts[i].getPos());
    }
    return returnVec;
}

/*
 * Gets the aligned length of the read. (Number of "M" bases.)
 */
int AlnLayout::getAlignedLen() {
    int len = 0;
    for(int i = 0; i < operations.size(); i++){
        if(operations[i].getOperation() == 'M') {
            len += operations[i].getLength();
        }
    }
    return len;
}

std::vector<CigarOp> AlnLayout::makeCigar(){
    std::vector<CigarOp> returnVec;
    for(int i = 0; i < operations.size(); i++){
        returnVec.push_back(CigarOp(operations[i].getOperation(),
                                    operations[i].getLength()));
    }
    return returnVec;
}


BamAlignment AlnLayout::toBam(){
    BamAlignment returnBam;
    returnBam.Name = Name;
    returnBam.QueryBases = seq;
    returnBam.RefID = RefID;
    returnBam.Position = pos;
    returnBam.SetIsReverseStrand(strandedness); // Keep the strandedness of the original read 1 for convention's sake.
    returnBam.CigarData = makeCigar(); // Create a std::Vector<CigarOp> object from the LayoutOperation objects.
    if(pairMerged){
        returnBam.AddTag("MP", "A", 'T'); // Add Merged Pair tag.
        returnBam.SetIsFirstMate(false); // Set that it is neither first nor second, as it is both mates
        returnBam.SetIsSecondMate(false);
    }
    returnBam.AddTag("PV", "Z", PhredStringFromVector(quality));
    returnBam.AddTag("FA", "Z", PhredStringFromVector(agreement));
    std::vector<std::string> Tags;
    boost::split(Tags, returnBam.TagData, boost::is_any_of("\t"));
    for(int i = 0; i < Tags.size(); i++){
        std::vector<std::string> TagElements;
        boost::split(TagElements, Tags[i], boost::is_any_of(":"));
        if(TagElements[0].compare("PV") == 0 || TagElements[0].compare("MP") == 0 || TagElements[0].compare("FA") == 0){
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

std::string AlnLayout::toBamStr(RefVector references){
    return BamToString(toBam(), references);
}

template <typename T>
inline std::vector<std::string> vecToStrVec(T inVec){
	std::vector<std::string>returnVec;
	for(int i = 0; i < inVec.size(); i++){
		returnVec.push_back(inVec[i].__str__());
	}
	return returnVec;
}
/**
std::vector<std::string> PosVecToStr(std::vector<LayoutPos> inVec){
	std::vector<std::string>returnVec;
	for(int i = 0; i < inVec.size(); i++){
		returnVec.push_back(inVec[i].__str__());
	}
	return returnVec;
} */

std::string AlnLayout::__str__(){
    std::string returnStr;
    std::vector<std::string> posVecStrs = vecToStrVec(getLayoutPositions());
    std::vector<std::string> opVecStrs = vecToStrVec(getOps());
    returnStr = boost::algorithm::join(posVecStrs, ",") + "|" + boost::algorithm::join(opVecStrs, ",");
    return returnStr;
}

/*
 * Creates a vector of tag strings from a BamAlignment object
 * because the BamAlignment strings aren't actually the strings
 * that would be created for a SAM file.
 */

std::vector<std::string> getBamTagVector(BamAlignment rec){
    std::vector<std::string> TagStrings;
    int tmpInt;
    float tmpFloat;
    std::string tmpString;
    if(rec.HasTag("AF")){
        rec.GetTag("AF", tmpFloat);
        TagStrings.push_back("AF:f:" + std::to_string(tmpFloat));
    }
    if(rec.HasTag("BS")){
        rec.GetTag("BS", tmpString);
        TagStrings.push_back("BS:Z:" + tmpString);
    }
    if(rec.HasTag("FM")){
        rec.GetTag("FM", tmpInt);
        TagStrings.push_back("FM:i:" + std::to_string(tmpInt));
    }
    if(rec.HasTag("FP")){
        rec.GetTag("FP", tmpInt);
        TagStrings.push_back("FP:i:" + std::to_string(tmpInt));
    }
    if(rec.HasTag("MQ")){
        rec.GetTag("MQ", tmpInt);
        TagStrings.push_back("MQ:i:" + std::to_string(tmpInt));
    }
    if(rec.HasTag("ND")){
        rec.GetTag("ND", tmpInt);
        TagStrings.push_back("ND:i:" + std::to_string(tmpInt));
    }
    if(rec.HasTag("NF")){
        rec.GetTag("NF", tmpFloat);
        TagStrings.push_back("NF:f:" + std::to_string(tmpInt));
    }
    if(rec.HasTag("NM")){
        rec.GetTag("NM", tmpInt);
        TagStrings.push_back("NM:i:" + std::to_string(tmpInt));
    }
    if(rec.HasTag("PV")){
        rec.GetTag("PV", tmpString);
        TagStrings.push_back("PV:Z:" + tmpString);
    }
    if(rec.HasTag("FA")){
        rec.GetTag("FA", tmpString);
        TagStrings.push_back("FA:Z:" + tmpString);
    }
    if(rec.HasTag("SV")){
        rec.GetTag("SV", tmpString);
        TagStrings.push_back("SV:Z:" + tmpString);
    }
    if(rec.HasTag("RP")){
        rec.GetTag("RP", tmpString);
        TagStrings.push_back("RP:Z:" + tmpString);
    }
    if(rec.HasTag("SC")){
        rec.GetTag("SC", tmpString);
        TagStrings.push_back("SC:Z:" + tmpString);
    }
    if(rec.HasTag("YA")){
        rec.GetTag("YA", tmpString);
        TagStrings.push_back("YA:Z:" + tmpString);
    }
    if(rec.HasTag("YO")){
        rec.GetTag("YO", tmpString);
        TagStrings.push_back("YO:Z:" + tmpString);
    }
    if(rec.HasTag("SF")){
        rec.GetTag("SF", tmpFloat);
        TagStrings.push_back("SF:f:" + std::to_string(tmpFloat));
    }
    if(rec.HasTag("X0")){
        rec.GetTag("X0", tmpInt);
        TagStrings.push_back("X0:i:" + std::to_string(tmpInt));
    }
    if(rec.HasTag("X1")){
        rec.GetTag("X1", tmpInt);
        TagStrings.push_back("X1:i:" + std::to_string(tmpInt));
    }
    if(rec.HasTag("XM")){
        rec.GetTag("XM", tmpInt);
        TagStrings.push_back("XM:i:" + std::to_string(tmpInt));
    }
    if(rec.HasTag("YX")){
        rec.GetTag("YX", tmpInt);
        TagStrings.push_back("YX:i:" + std::to_string(tmpInt));
    }
    if(rec.HasTag("YM")){
        rec.GetTag("YM", tmpInt);
        TagStrings.push_back("YM:i:" + std::to_string(tmpInt));
    }
    if(rec.HasTag("YQ")){
        rec.GetTag("YQ", tmpInt);
        TagStrings.push_back("YQ:i:" + std::to_string(tmpInt));
    }
    if(rec.HasTag("YR")){
        rec.GetTag("YR", tmpInt);
        TagStrings.push_back("YR:i:" + std::to_string(tmpInt));
    }
    if(rec.HasTag("MP")){
        char op;
        rec.GetTag("MP", op);
        TagStrings.push_back("MP:A:" + op);
    }
    return TagStrings;
}

std::string getBamTagString(BamAlignment rec){
	std::vector<std::string> bamtags = getBamTagVector(rec);
    return boost::algorithm::join(bamtags, "\t");
}


AlnLayout::AlnLayout(BamAlignment rec) {
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
    TagData = getBamTagString(rec);
    std::string PVString;
    rec.GetTag("PV", PVString);
    std::string FAString;
    rec.GetTag("FA", FAString);
    quality = PhredVectorFromString(PVString);
    agreement = PhredVectorFromString(FAString);
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

int main(int argc, char* argv[]){
    std::string inputBam, outputBam;
    if(argc == 3) {
        inputBam = argv[1];
        outputBam = argv[2];
    }
    else if(argc == 2){
        inputBam = argv[1];
        std::vector<std::string> splitOut;
        boost::split(splitOut, inputBam, boost::is_any_of("."));
        outputBam = "";
        for(int i = 0; i < splitOut.size() - 1; i++){
            outputBam += splitOut[i] + ".";
        }
        outputBam += "out.bam";
    }
    else if(argc == 1){
        printf("Usage: I don't know. Position arguments: InputBam, OutputBam.\n");
        printf("Output bam optional - replaces .bam with .out.bam in the filename.\n");
    }
    std::cerr << "Output BAM :" << outputBam;
    BamReader reader;
    if(!reader.Open(inputBam)) {
        std::cerr << "Could not open input BAM" << std::endl;
        return 1;
    }
    std::cerr << "Opened bam reader" << std::endl;
    const SamHeader header = reader.GetHeader();
    const RefVector references = reader.GetReferenceData();
    BamAlignment rec;
    BamWriter writer;
    if(!writer.Open(outputBam, header, references) ) {
        std::cerr << "ERROR: could not open " + outputBam + " for writing. Abort mission!" << std::endl;
        throw std::runtime_error("Could not open " + outputBam + " for writing.");
    }
    std::cerr << "Opened bam writer" << std::endl;
    while (reader.GetNextAlignment(rec)) {
        //std::cout << BamToString(rec, references) << std::endl; // This checked to see if BAM formats were being parsed and reformatted correctly.
    	std::cerr << "I'm about to make an AlnLayout object." << std::endl;
    	AlnLayout newAlnLayout = AlnLayout(rec);
    	std::cerr << "I'm about to get the first operation from that object." << std::endl;
    	LayoutOp newOp = newAlnLayout.getOps()[0];
    	std::cout << AlnLayout(rec).getOps()[0].__str__() << std::endl;
    }
    reader.Close();
    writer.Close();
    return 0;
}
