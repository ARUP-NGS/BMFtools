//Hey this is a C++ file!

#include <iostream>
#include "api/BamReader.h"
#include "api/BamWriter.h"

using namespace std;
using namespace BamTools;

int main(int argc, char* argv[])
{
    string inputFilename;
    if(argc == 2) {
        inputFilename = argv[1];
    }
    else {
        cerr << "Too many arguments" << endl;
        return 1;
    }
    BamReader reader;
    if(!reader.Open(inputFilename)) {
        cerr << "Could not open input BAM" << endl;
        return 1;
    }
    const BamIndex::IndexType type = BamIndex::STANDARD;
    if(! reader.CreateIndex(type)) {
        cerr << "Create index failed!" << endl;
        return 1;
    }
    cerr << "Index created successfully!" << endl;
    reader.Close();
    return 0;
}
