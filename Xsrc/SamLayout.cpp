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

int LayoutPos::getAgreement() {return agreement;}

bool LayoutPos::getIsMerged() {
    return isMerged;
}

void LayoutPos::setIsMerged(bool newMerged) {isMerged = newMerged;}

void LayoutPos::setAgreement(int newAgreement) {agreement = newAgreement;}
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
    ss << ";IsMerged:" << isMerged;
    return ss.str();
}

void LayoutPos::setAttributes(char baseArg, char opArg, int RefIDArg, int posArg,
                              int readPosArg, int qualArg, int strand, int FA) {
    Operation = opArg;
    base = baseArg;
    RefID = RefIDArg;
    Position = posArg;
    ReadPosition = readPosArg;
    quality = (base != 'N') ? qualArg: 0; // Turn all quality scores for N bases to 0
    strandedness = strand;
    agreement = FA;
    isMerged = false;
}

LayoutPos::LayoutPos(char baseArg, char opArg, int RefIDArg, int posArg, int readPosArg, int qualArg, int strand, int FA) {
    setAttributes(baseArg, opArg, RefIDArg, posArg, readPosArg, qualArg, strand, FA);
}

int LayoutPos::getStrandedness(){
    return strandedness;
}

void LayoutPos::setStrandedness(int newStrandedness) {strandedness = newStrandedness;}

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

std::string LayoutOp::getSequence() {
    updateSequence();
    return seq;
}
void LayoutOp::setSequence(std::string newSeq) {seq = newSeq;}

void LayoutOp::update(){
    updateAgreement();
    updateQuality();
    updateSequence();
}

void LayoutOp::setAgreement(std::vector<int> newAgreements){
    if(newAgreements.size() != layoutPositions.size()){
        throw std::runtime_error("Length of new agreements vector is not the same as the number of l");
    }
    agreement = newAgreements;
    int i = 0;
    for(LayoutPos& position: layoutPositions){
        position.setAgreement(newAgreements[i]);
        i++;
    }
}

void LayoutOp::updateAgreement(){
    std::vector<int> tmpAgreements;
    for(auto& lPos : layoutPositions){
        tmpAgreements.push_back(lPos.getAgreement());
    }
    setAgreement(tmpAgreements);
}
void LayoutOp::updateQuality(){
    std::vector<int> tmpQualities;
    for(auto& lPos : layoutPositions){
        tmpQualities.push_back(lPos.getQuality());
    }
    setQuality(tmpQualities);
}
void LayoutOp::updateSequence(){
    std::string returnStr;
    for(LayoutPos& lPos : layoutPositions){
        returnStr += lPos.getBase();
    }
    setSequence(returnStr);
}

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
    char base;
    int baseQuality, agreeCount, Pos, tmpRPos;
    /* char baseArg, char opArg, int RefIDArg, int posArg,
                              int readPosArg, int qualArg, int strand, int FA*/
    LayoutPos ALP;
    for(int k = 0; k < quality.size(); k++) {
        //std::cerr << "Now making ALP # " << k + 1 << std::endl;
        switch(Operation){
        case 'S':
            baseQuality = quality[k];
            agreeCount = agreement[k];
            Pos = -1;
            tmpRPos = readPos + k;
            base = seq.at(k);
            break;
        case 'D':
            Pos = pos + k;
            tmpRPos = -1;
            agreeCount = -1;
            base = 'D';
            baseQuality = -1;
            break;
        case 'I':
            base = seq.at(k);
            baseQuality = quality[k];
            agreeCount = agreement[k];
            Pos = -1;
            tmpRPos = readPos + k;
            break;
        case 'M':
            base = seq.at(k);
            baseQuality = quality[k];
            agreeCount = agreement[k];
            Pos = pos + k;
            tmpRPos = readPos + k;
            break;
        default:
            throw std::runtime_error("Sorry, unsupported cigar character. Email me and I'll change this.");
        }
        ALP = LayoutPos(base, Operation, RefID, Pos, tmpRPos, baseQuality, strandedness, agreeCount);
        //std::cerr << ALP.__str__() << " is ALP string." << std::endl;
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
            ALP.setAttributes(seq.at(k), Operation, RefID, -1 * pos, readPos + k, quality[k], strandedness, 1);
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


std::vector<int> IntVectorFromString(std::string PVString){
    std::vector<int> returnVec;
    std::vector<std::string> stringVec;
    boost::split(stringVec, PVString, boost::is_any_of(","));
    int vecSize = stringVec.size();
    for(int i = 0; i < vecSize; i++){
        returnVec.push_back(std::stoi(stringVec[i]));
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

std::vector<int> sliceVector(std::vector<int> inVec, int start, int end){
    return std::vector<int>(&inVec[start], &inVec[end]);
}
/*
template <typename T>
std::vector<T> sliceVector(std::vector<T> inVec, int start, int end){
    return std::vector<T>(&inVec[start], &inVec[end]);
} */

std::vector<int> InitializeVector(int length, int value){
    std::vector<int> returnVec;
    for(int i = 0; i < length; i++){
        returnVec.push_back(value);
    }
    return returnVec;
}

std::string InitializeString(int length, char character){
    std::stringstream ss;
    for(int i = 0; i < length; i++){
        ss >> character;
    }
    return ss.str();
}


std::vector<LayoutOp> GetLayoutOps(BamAlignment rec) {
    int cigarLen;
    int StartPosition = rec.Position; // Needed to keep read indices in line with reference.
    int cigOffset = 0;
    // These next variables are ones I hope to remove once I figure out what I need
    // to initialize LayoutOp.
    std::vector<LayoutOp> operations;
    LayoutOp TmpOp;
    CigarOp TmpCigarOp;
    cigarLen = rec.CigarData.size();
    std::string PVString;
    for(int i = 0; i < cigarLen; i++) {
        TmpCigarOp = rec.CigarData[i];
        if(rec.HasTag("PV")) {
            std::string PVString;
            rec.GetTag("PV", PVString);
            std::string FAString;
            rec.GetTag("FA", FAString);
            std::string querySub;
            std::vector<int> PhredVector, AgreementVector, PhredSubVector, AgreementSubVector;
            if(TmpCigarOp.Type != 'D'){
                PhredVector = IntVectorFromString(PVString);
                AgreementVector = IntVectorFromString(FAString);
                PhredSubVector = sliceVector(PhredVector, cigOffset, cigOffset + TmpCigarOp.Length);
                AgreementSubVector = sliceVector(AgreementVector, cigOffset, cigOffset + TmpCigarOp.Length);
                querySub = rec.QueryBases.substr(cigOffset, TmpCigarOp.Length);
                assert(querySub.size() == TmpCigarOp.Length);
            }
            else {
                PhredSubVector = InitializeVector(TmpCigarOp.Length, -1);
                AgreementSubVector = InitializeVector(TmpCigarOp.Length, -1);
                querySub = InitializeString(TmpCigarOp.Length, 'D');
            }
            LayoutOp tmpOp = LayoutOp(querySub, //Read sequence
                    PhredSubVector, AgreementSubVector, TmpCigarOp.Type, rec.RefID, StartPosition + cigOffset,
                    cigOffset, rec.IsReverseStrand() ? -1 : 1);
            operations.push_back(tmpOp);
            if(tmpOp.getOperation() != 'D'){
                cigOffset += TmpCigarOp.Length;
            }
        }
        else {
            std::cerr << "This record has no PV tag" << std::endl;
            // Gets quality scores from ASCII phred scores rather than PV numbers.
            LayoutOp tmpOp = LayoutOp(rec.QueryBases.substr(cigOffset, TmpCigarOp.Length), // Read sequence
                    rec.Qualities.substr(cigOffset, TmpCigarOp.Length), // Read qualities
                    TmpCigarOp.Type, // Cigar operation
                    rec.RefID, // Reference ID
                    StartPosition + cigOffset, // Genomic coords for start of cigar
                    cigOffset, // Read coords for start of cigar
                    rec.IsReverseStrand() ? -1 : 1
                    );
            operations.push_back(tmpOp);
            if(tmpOp.getOperation() != 'D'){
                cigOffset += TmpCigarOp.Length;
            }

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
void AlnLayout::updateLen() {length = getOps().size();}
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
 * Changed: Now outputs the incremented SAM position,
 * not the 0-based BAM.
 */
std::string BamToString(BamAlignment rec, RefVector vec){
    std::stringstream ss;
    std::string returnStr;
    std::string TagData = getBamTagString(rec);
    ss << rec.Name + "\t" << rec.AlignmentFlag << "\t" + vec[rec.RefID].RefName
    << "\t" << rec.Position + 1 << "\t" << rec.MapQuality << "\t"
    << CigarDataToString(rec.CigarData) << "\t" << vec[rec.MateRefID].RefName
    << "\t" << rec.MatePosition + 1 << "\t" << rec.InsertSize << "\t"
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

std::string IntVecToPhred33(std::vector<int> inVec){
    std::string returnStr;
    for(auto &i: inVec){
        if(i <= 93){
        returnStr += char(i + 33);
        }
        else {
            returnStr += char(126);
        }
    }
    return returnStr;
}


BamAlignment AlnLayout::toBam(BamAlignment templateAln){
    BamAlignment returnBam = BamAlignment(templateAln);
    returnBam.Name = Name;
    returnBam.QueryBases = seq;
    returnBam.RefID = RefID;
    returnBam.Position = pos;
    returnBam.SetIsReverseStrand(strandedness); // Keep the strandedness of the original read 1 for convention's sake.
    returnBam.CigarData = makeCigar(); // Create a std::Vector<CigarOp> object from the LayoutOperation objects.
    returnBam.Qualities = IntVecToPhred33(quality);
    returnBam.SetIsMateReverseStrand(templateAln.IsMateReverseStrand());
    /*
    if(pairMerged){
        std::cerr << "Please make this merged bam record" << std::endl;
    	tmpTag = 'T';
        if(returnBam.HasTag("MP")){
            returnBam.EditTag("MP", "A", tmpTag); // Add Merged Pair tag.
        }
        else {
            returnBam.AddTag("MP", "A", tmpTag); // Add Merged Pair tag.
        }
        returnBam.SetIsFirstMate(false); // Set that it is neither first nor second, as it is both mates
        returnBam.SetIsSecondMate(false);
    }
    else {
    	tmpTag = 'F';
        std::cerr << "Please make this bam record" << std::endl;
        if(returnBam.HasTag("MP")){
        	std::cerr << "It has the tag!" << std::endl;
            returnBam.EditTag("MP", "A", tmpTag);
        }
        else {
        	std::cerr << "It does not have the tag!" << std::endl;
        	returnBam.AddTag("MP", "A", tmpTag);
        }
    } */
    returnBam.EditTag("PV", "Z", PhredStringFromVector(quality));
    returnBam.EditTag("FA", "Z", PhredStringFromVector(agreement));
    returnBam.EditTag("DB", "Z", PhredStringFromVector(mergedPositions));
    std::vector<std::string> Tags = getBamTagVector(templateAln);
    for(int i = 0; i < Tags.size(); i++){
        std::vector<std::string> TagElements;
        boost::split(TagElements, Tags[i], boost::is_any_of(":"));
        //std::cerr << "Now adding tag " << TagElements[0] << std::endl;
        if(TagElements[0].compare("PV") == 0 || TagElements[0].compare("MP") == 0 || TagElements[0].compare("FA") == 0){
            continue;
        }
        switch (TagElements[1].at(0)){
                case 'A':
                   returnBam.EditTag(TagElements[0], "A", (char) TagElements[2].at(0));
                int tmpInt;
                float tmpFloat;
                break;
                    case 'Z':
                        returnBam.EditTag(TagElements[0], "Z", TagElements[2]);
                        break;
                    case 'i':
                        tmpInt = std::stoi(TagElements[2]);
                        returnBam.EditTag(TagElements[0], "i", tmpInt);
                        break;
                    case 'f':
                        tmpFloat = std::stof(TagElements[2]);
                        returnBam.EditTag(TagElements[0], "f", tmpFloat);
                        break;
                    default:
                        throw std::runtime_error("Hey, what is this format?");
        }
    }
    // std::cerr << "Now returning the bam with name " << returnBam.Name << "and is it R1? " << returnBam.IsFirstMate() << std::endl;
    return returnBam;
}


std::string AlnLayout::toBamStr(BamAlignment templateBam,
                                RefVector references){
    return BamToString(toBam(templateBam), references);
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
    returnStr = boost::algorithm::join(posVecStrs, "&") + "|" + boost::algorithm::join(opVecStrs, "&");
    return returnStr;
}

void AlnLayout::updateSequence() {

}

void AlnLayout::update(){

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

int AlnLayout::getFirstAlignedBase(BamAlignment templateRec){
    return getastart(toBam(templateRec));
}

int AlnLayout::getFirstAlignedBase(){
    int cigSum = 0;
    for(LayoutOp& op : operations){
        if(op.getOperation() == 'S'){
            return -1 * op.getPos(); // Returns the first reference base to which it is aligned
        }
        else if(op.getOperation() == 'M'){
            return op.getPos();
        }
        else {
        	continue;
        }
    }
    throw std::runtime_error("Working under the assumption that at least one cigar operation is mapped or softclipped.");
}

int AlnLayout::getLastMappedBase(){
    int cigSum = 0;
    int decrementor = operations.size();
    LayoutOp op;
    for(int decrementor = operations.size() - 1; decrementor > -1; decrementor--){
        op = operations[decrementor];
        if(op.isMap()){
            return op.getPos() + op.getLength() - 1; // Add length of the cigar element, minus one because the position is already set for the first base in the set.
        }
    }
    throw std::runtime_error("Working under the assumption that at least one cigar operation is mapped.");
}

/*
 * Returns -1 for none, otherwise the first
 * position "mapped" to a reference, whether
 * soft-clipped or actually considered matched.
 */

int AlnLayout::getFirstMatchingRef(){
    for(LayoutOp& op: operations){
        if(op.getOperation() == 'M'){
            return op.getPos();
        }
        if(op.getOperation() == 'S'){
            return -1 * op.getPos();
        }
    }
    return -1;
}

int getastart(BamAlignment rec){
    if(!rec.IsMapped()) {
        return -1;
    }
    int cigSum = 0;
    for(auto& cig: rec.CigarData){
        if(cig.Type == 'M'){
            return cigSum;
        }
        cigSum += cig.Length;
    }
    return cigSum;
}

/*
 * Empty constructor
 */

AlnLayout::AlnLayout(){
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
    quality = IntVectorFromString(PVString);
    agreement = IntVectorFromString(FAString);
    strandedness = rec.IsReverseStrand() ? -1 : 1;
    mateStrandedness = rec.IsMateReverseStrand() ? -1 : 1;
    operations = GetLayoutOps(rec);
    firstAlignedBase = getFirstAlignedBase(rec);
    //std::cerr << "Finished initializing an AlnLayout from a bam rec" << std::endl;
}

std::string AlnLayout::getTagString() {return TagData;}
void AlnLayout::setTagString(std::string newTagStr) {
    TagData = newTagStr;
}

std::string AlnLayout::getName() {return Name;}
void AlnLayout::setName(std::string newName) {
    Name = newName;
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

std::string CigarOpToStr(CigarOp CigarData){
    std::stringstream ss;
    ss >> CigarData.Length >> CigarData.Type;
    return ss.str();
}


std::string CigarDataToStr(std::vector<CigarOp> CigarData){
    std::string returnStr;
    for(CigarOp& op: CigarData){
        returnStr += CigarOpToStr(op);
    }
    return returnStr;
}

int firstIntersect(AlnLayout R1Layout, AlnLayout R2Layout){
    int R2First;
    R2First = R2Layout.getFirstAlignedBase();
    for(LayoutPos& pos: R1Layout.getLayoutPositions()){
        if(std::abs(pos.getPos()) == R2First){ // Both reads are mapped here in some fashion
            return R2First;
        }
        if(std::abs(pos.getPos()) < R2First) {
            throw std::runtime_error("Somehow read 2's first base got missed when checking for use of ");
        }
    }
    std::cerr << "No intersection found between pairs..." << std::endl;
    return -1; // No match found...
}

LayoutPos mergePositions(LayoutPos pos1, LayoutPos pos2){
    int newQual, newAgreement, newStrandedness;
    char newBase;
    if(pos1.getOperation() != pos2.getOperation()){
        throw std::runtime_error("Looks like these two have different operations. Read pairs disagree, give up!");
    }
    if(pos1.getPos() != pos2.getPos()){
        throw std::runtime_error("Looks like these two have different positions. Read pairs disagree, give up!");
    }
    if(pos1.getBase() == pos2.getBase()){
        newQual = pos1.getQuality() + pos2.getQuality();
        newAgreement = pos1.getAgreement() + pos2.getAgreement();
        newBase = pos1.getBase();
    }
    else {
        if(pos1.getQuality() > pos2.getQuality()){
            newQual = pos1.getQuality() - pos2.getQuality();
            newBase = pos1.getBase();
            newAgreement = pos1.getAgreement();
        }
        else {
            newQual = pos2.getQuality() - pos1.getQuality();
            newBase = pos2.getBase();
            newAgreement = pos2.getAgreement();
        }
    }
    if(newQual < 3){
        newBase = 'N';
    }
    newStrandedness = 0;
    LayoutPos returnLP = LayoutPos(newBase, pos1.getReferenceID(), pos1.getOperation(), pos1.getPos(), -1, newQual, newStrandedness, newAgreement);
    returnLP.setIsMerged(true);
    return returnLP;
}

/*
 *
 */

std::vector<LayoutPos> MergeLayouts(AlnLayout R1Layout, AlnLayout R2Layout){
    // Find the one which starts earlier on the contig
    int firstR2Base, cigOffset, opOffset;
    std::vector<LayoutPos> R1PosV, R2PosV, NewLayouts;
    LayoutPos tmpPos;
    if(R1Layout.getFirstMatchingRef() > R2Layout.getFirstMatchingRef()){
        std::swap(R1Layout, R2Layout); // Start laying out the first one
    }
    firstR2Base = R2Layout.getFirstAlignedBase();
    if(R1Layout.getLastMappedBase() < firstR2Base){
    	return NewLayouts;
    }
    // Go LayoutOp by LayoutOp until there's some overlap
    R1PosV = R1Layout.getLayoutPositions();
    R2PosV = R2Layout.getLayoutPositions();
    for(opOffset = 0; opOffset < R1PosV.size(); opOffset++){
        if(R1PosV[opOffset].getPos() == firstR2Base){
            break;
        }
        NewLayouts.push_back(R1PosV[opOffset]);
    }
    int countOffset = 0;
    while(opOffset < R1PosV.size()){
        try {
            tmpPos = mergePositions(R1PosV[opOffset], R2PosV[countOffset]);
        }
        catch(const std::exception& ex){
            std::vector<LayoutPos> emptyVec;
            return emptyVec; // If this returns an empty vector, then we just couldn't successfully merge the two in the pair.
        }
        tmpPos.setReadBasePos(opOffset);
        NewLayouts.push_back(tmpPos);
        countOffset++;
        opOffset++;
    }
    for(int i = countOffset;i < R2PosV.size();i++){
        // Append all the layout position objects from the part that don't overlap.
        NewLayouts.push_back(R2PosV[i]);
    }
    return NewLayouts;
}

AlnLayout::AlnLayout(std::vector<LayoutPos> lPositions, bool merged=true){
    std::string cigarOpSeq;
    std::vector<int> sliceQual, sliceAgreement;
    std::vector<LayoutOp> tmpOperations;
    int cigOffset = 0;
    char workingOp = lPositions[0].getOperation();
    int workingOpCount = 0;
    Name = "MergedLayoutPositionListsInto an AlnLayoutObject. Change me!";
    seq = "";
    TagData = "";
    pos = std::abs(lPositions[0].getPos());
    RefID = lPositions[0].getReferenceID();
    cigarOpSeq = "";
    for(int i = 0; i < lPositions.size(); i++){
        if(lPositions[i].getIsMerged()){
            mergedPositions.push_back(i);
        }
        switch(lPositions[i].getOperation()){
            case 'I':
                quality.push_back(lPositions[i].getQuality());
                seq += lPositions[i].getBase();
                cigarOpSeq += lPositions[i].getBase();
                agreement.push_back(lPositions[i].getAgreement());
                break;
            case 'D':
                break;
            case 'S':
                quality.push_back(lPositions[i].getQuality());
                seq += lPositions[i].getBase();
                cigarOpSeq += lPositions[i].getBase();
                agreement.push_back(lPositions[i].getAgreement());
                break;
            case 'M':
                quality.push_back(lPositions[i].getQuality());
                seq += lPositions[i].getBase();
                cigarOpSeq += lPositions[i].getBase();
                agreement.push_back(lPositions[i].getAgreement());
                break;
        }
        if(lPositions[i].getOperation() == workingOp){
            workingOpCount++;
            continue;
        }
        else { // Now make a LayoutOp from this list of LayoutPos objects.
            sliceQual = sliceVector(quality, cigOffset, cigOffset + cigarOpSeq.size());
            sliceAgreement = sliceVector(agreement, cigOffset, cigOffset + cigarOpSeq.size());
            switch(workingOp) {
            case 'M':
                operations.push_back(LayoutOp(cigarOpSeq, sliceQual, sliceAgreement, workingOp, lPositions[0].getReferenceID(), std::abs(lPositions[cigOffset].getPos()) + cigOffset, cigOffset, lPositions[cigOffset].getStrandedness()));
                break;
            case 'D':
                operations.push_back(LayoutOp(InitializeString(cigarOpSeq.size(), 'D'), InitializeVector(cigarOpSeq.size(), -1), InitializeVector(cigarOpSeq.size(), -1),
                                              workingOp, lPositions[i].getReferenceID(),
                                              -1, -1, lPositions[0].getStrandedness()));
                break;
            case 'I':
                operations.push_back(LayoutOp(cigarOpSeq, sliceQual, sliceAgreement, workingOp, lPositions[cigOffset].getReferenceID(), -1, cigOffset, lPositions[cigOffset].getStrandedness()));
                break;
            case 'S':
                operations.push_back(LayoutOp(cigarOpSeq, sliceQual, sliceAgreement, workingOp, lPositions[cigOffset].getReferenceID(), std::abs(lPositions[cigOffset].getPos()), cigOffset, lPositions[cigOffset].getStrandedness()));
                break;
            default:
            	int cigarCharEq = (int) workingOp;
            	std::cerr << "Cigar character: " << workingOp << ". ASCII character code: " << cigarCharEq << std::endl;
                throw std::runtime_error("Unsupported cigar character.");
            }
            // Reset temporary variables
            workingOpCount = 0;
            cigarOpSeq = "";
            cigOffset += cigarOpSeq.size();
        }
    }
    pairMerged = true;
    updateLen(); // Set the new length for the AlnLayout object to be equal to the number of cigar operations
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
        boost::split(splitOut, outputBam, boost::is_any_of("/"));
        outputBam = splitOut[splitOut.size() - 1]; // Make it a relative rather than absolute path
    }
    else if(argc == 1){
        printf("Usage: I don't know. Position arguments: InputBam, OutputBam.\n");
        printf("Output bam optional - replaces .bam with .out.bam in the filename.\n");
    }
    std::cerr << "Output BAM :" << outputBam << std::endl;
    BamReader reader;
    bool openSuccess = reader.Open(inputBam);
    if(!openSuccess) {
        std::cerr << "Could not open input BAM" << std::endl;
        return 1;
    }
    std::cerr << "Opened bam reader. " << std::endl;
    const SamHeader header = reader.GetHeader();
    const RefVector references = reader.GetReferenceData();
    BamAlignment rec, rec1, rec2;
    BamWriter writer;
    bool writerSuccess = writer.Open(outputBam, header, references);
    int tlen;
    if(!writerSuccess) {
        std::cerr << "ERROR: could not open " + outputBam + " for writing. Abort mission!" << std::endl;
        throw std::runtime_error("Could not open " + outputBam + " for writing.");
    }
    std::cerr << "Opened bam writer" << std::endl;
    AlnLayout newAln;
    int ncount = 0;
    std::string tagTrue = "True";
    std::string tagFalse = "False";
    while (reader.GetNextAlignment(rec)) {
        ncount++;
        if(rec.IsFirstMate()){
            rec1 = rec;
            continue;
        }
        else if(rec.IsSecondMate()){
            rec2 = rec;
        }
        //std::cout << BamToString(rec, references) << std::endl; // This checked to see if BAM formats were being parsed and reformatted correctly.
        //if(rec.CigarData.size() == 0){
        if(!rec1.IsMapped() || !rec2.IsMapped()){
            writer.SaveAlignment(rec1);
            writer.SaveAlignment(rec2);
            continue;
        }
        if(rec1.RefID != rec2.RefID){
            writer.SaveAlignment(rec1); // If pairs aren't on the same contig, don't merge them.
            writer.SaveAlignment(rec2);
            continue;
        }
        tlen = std::abs(rec1.Position - rec2.Position);
        if(tlen > MaxInsert){
            writer.SaveAlignment(rec1);
            writer.SaveAlignment(rec2);
            continue;
        }
        if(rec1.Name != rec2.Name) {
            throw std::runtime_error("Read names are not the same for a pair!");
        }
        std::vector<LayoutPos> mergedPositions;
        AlnLayout AlnL1 = AlnLayout(rec1);
        AlnLayout AlnL2 = AlnLayout(rec2);
        mergedPositions = MergeLayouts(AlnL1, AlnL2);
        /*
        catch(const std::exception& ex) {
        	std::cerr << "Exception " << ex.what() << std::endl;
            std::cerr << "Failed to merge reads of name because of an exception: " << rec.Name << std::endl;
            writer.SaveAlignment(rec1);
            writer.SaveAlignment(rec2);
            continue;
        } */
        if(mergedPositions.size() == 0){ // Length is 0 if merging failed.
            //std::cerr << "Failed to merge reads of name " << rec.Name << " for not being mergable." <<std::endl;
        	rec1.AddTag("MP", "Z", tagFalse);
        	rec2.AddTag("MP", "Z", tagFalse);
            writer.SaveAlignment(rec1);
            writer.SaveAlignment(rec2);
            continue;
        }
        newAln = AlnLayout(mergedPositions, true);
        //std::cerr << "newAln string repr: " << newAln.__str__() << std::endl;
        newAln.setName(AlnL1.getName() + "PairMerged");
        newAln.setStrandedness(0);
        rec = newAln.toBam(rec);

        rec.AddTag("MP", "Z", tagTrue);
        //std::cerr << "Made bam record from newAln!" << std::endl;
        writer.SaveAlignment(rec); // Save merged record
    }
    bool closeSuccess = reader.Close();
    writer.Close();
    return 0;
}
