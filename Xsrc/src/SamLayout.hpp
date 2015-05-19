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
