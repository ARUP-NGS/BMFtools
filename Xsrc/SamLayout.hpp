class LayoutPos {
    /** Contains information regarding an alignment layout with nucleotide resolution.
     */
    private:
        char Operation; //cigar operation for this base. Set to "Y" for uninitialized.
        int RefID;
        int Position; // Set to -1 for an insertion or soft-clipping
        int ReadPosition; // Set to -1 for a deletion
        char base; // Set to "D" for a deletion.
        int quality; // Set to -1 for a deletion.
        int strandedness; // -1 for reverse, 1 for forward, 0 for neither.

    public:
        LayoutPos(char, char, int, int, int, int, int); // Base, Op, ref id, Position (ref), Position (read), Quality, Strand
        LayoutPos();
        void setAttributes(char, char, int, int, int, int, int);
        int getReferenceID();
        void setReferenceID(int);
        int getReadBasePos(); // Returns which base in a read this tuple is.
        void setReadBasePos(int);
        char getBase();
        void setBase(char);
        char getOperation();
        void setOperation(char);
        int getPos();
        void setPos(int);
        int getQuality();
        void setQuality(int);
        void check(); // More of this check could be filled out, but let's just get this compiling.
        bool incrementRefPos();
        bool incrementReadPos();
};


class LayoutOp {
    /**
     * Contains the information needed from a CigarOp object and a read to
     * make a list of LayoutPos. Contrasting with LayoutPos
     * objects, this has resolution equal to a cigar operation.*/
    private:
        char Operation; // Character definition of cigar operation.
        int pos; // 0-based genomic start position of the cigar operation.
        int readPos; // 0-based read start position of the cigar operation.
        int RefID; // Reference ID.
        std::string seq;
        std::vector<int> quality; // Store the information as integers. For this purpose, use the PV tags if available.
        std::vector<LayoutPos> layoutPositions;
        int strandedness;

    public:
        LayoutOp(std::string, std::string, char, int, int, int, int); //Constructor
        LayoutOp();
        LayoutOp(std::string, std::vector<int>, char, int, int, int, int); //For PV reads.

        void setAttributes(std::string, std::string, char, int, int, int, int); // For setting values if declared before constructed.
        void setAttributes(std::string, std::vector<int>, char, int, int, int, int);

        bool isIns();
        bool isDel();
        bool isMap();
        bool isSoftClipped();

        char getOperation();
        void setOperation(char);

        int getLength();

        int getReadPos();
        void setReadPos(int);

        int getPos();
        void setPos(int);

        int getRef();
        void setRef(int);

        int getStrandedness();
        void setStrandedness(int);

        std::string getSequence();
        void setSequence(std::string);

        std::vector<int> getQuality();
        void setQuality(std::vector<int>);
        void clearQuality();
        void mergeOp(LayoutOp);
        std::vector<LayoutPos> getLayoutPositions();
        void replacePos(LayoutPos, int);
        bool incrementRefPos(); // Returns true if the operation requires the incrementing of the reference position in the LayoutOp
        bool incrementReadPos(); // Returns true if the operation requires the increment
        void updatePositions();
        void updateReadPositions();
        void update();
};
