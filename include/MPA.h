/*
 * ArrayLayout header!
 */

typedef struct ArrayLayoutPos{
    int pos;
    int readPos;
    int quality;
    int agreement;
    char operation;
    char base;
    char mergeAgreed;
} ArrayLayoutPos_t;

typedef struct ArrayLayout {
	ArrayLayoutPos_t * layouts;
	int length;
} ArrayLayout_t;

ArrayLayoutPos_t cMergeLayoutPositions(ArrayLayoutPos_t L1, ArrayLayoutPos_t L2);
int getFirstAlignedRefPos(ArrayLayout_t layout);
ArrayLayout_t MergeLayouts(ArrayLayout_t L1, ArrayLayout_t L2);
