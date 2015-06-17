/*
 * ArrayLayout header!
 */
#include "stdint.h"

typedef struct ArrayLayoutPos{
    int pos;
    uint16_t readPos;
    int quality;
    uint16_t agreement;
    char operation;
    char base;
    char mergeAgreed;
} ArrayLayoutPos_t;

typedef struct ArrayLayout {
	ArrayLayoutPos_t * layouts;
	int length;
} ArrayLayout_t;

typedef struct MergeRet {
	ArrayLayout_t Layout;
	char Success;
} MergeRet_t;

ArrayLayoutPos_t cMergeLayoutPositions(ArrayLayoutPos_t L1, ArrayLayoutPos_t L2);
int getFirstAlignedRefPos(ArrayLayout_t layout);
ArrayLayout_t MergeLayouts(ArrayLayout_t L1, ArrayLayout_t L2);
//ArrayLayout_t MergeWithPassFail(ArrayLayout_t AL1, ArrayLayout_t AL2);
