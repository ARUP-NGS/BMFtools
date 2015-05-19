//header file for pairedBT.cpp
//It tags ur BAMs!
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

std::vector<std::string> splitString(std::string inStr, std::string delimiter); // Splits a string into a string vector with a string delimiter. Does not work!
std::vector<std::string> splitByChar(std::string inStr, char delimiter); // Splits a string into a string vector with a character delimiter. Does work!
std::string CoorString(RefVector refs, BamAlignment rec1, BamAlignment rec2); // Returns a unique string for each set of coordinate starts.
float SoftClippedFraction(BamAlignment rec); // Returns the fraction of bases in read softclipped.
float AlignedFraction(BamAlignment rec); // Returns the fraction of bases aligned.
