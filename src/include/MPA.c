#include "MPA.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

ArrayLayoutPos_t cMergeLayoutPositions(ArrayLayoutPos_t L1, ArrayLayoutPos_t L2){
    if(L1.base == L2.base){
        if(L1.operation == L2.operation || L2.operation == 83){
            L1.agreement += L2.agreement;
            L1.quality += L2.quality;
            L1.mergeAgreed = 2;
        }
        else if(L2.base == 78){
            return L1;
        }
        else {
            L1.base = 78; // Set base to N
            L1.mergeAgreed = 0;
            L1.agreement = 0;
        }
    }
    else {
        if(L1.operation == L2.operation){
            if(L1.quality > L2.quality){
                L1.base = L2.base;
                L1.quality = L2.quality - L1.quality;
                L1.mergeAgreed = 0;
            }
            else {
                L1.quality -= L2.quality;
                L1.mergeAgreed = 0;
            }
        }
        else {
            L1.mergeAgreed = -1;
        }
    }
    return L1;
}

int getFirstAlignedRefPos(ArrayLayout_t layout){
    for(int i = 0; i < layout.length; i++){
        if(layout.layouts[i].operation == 77){
            // If operation is "M"
            return layout.layouts[i].pos - i;
        }
    }
    return -1;
}


char MergeOverlappedLayouts(ArrayLayout_t AL1, ArrayLayout_t AL2){
    int start1, start2, offset;
    start1 = getFirstAlignedRefPos(AL1);
    start2 = getFirstAlignedRefPos(AL2);
    /*
    //You'll want to make sure that AL1's start position is less than
    //AL2's since I'm trying to avoid this kind of moving.
    if(start1 > start2){
        // Switch the order to make AL1 the first read.
        ArrayLayout_t tmpLayout = AL2;
        AL2 = AL1;
        AL1 = tmpLayout;
        int tmpInt = start2;
        start2 = start1;
        start1 = tmpInt;
    } */
    offset = start2 - start1;
    // Now figuring out precisely how long this is going to have to be
    for(int i = 0; i < offset; i++){
        if(AL1.layouts[i].operation == 68){
            // Correct for "D" operations
            offset += 1;
        }
    }
    // AL1.length = offset + AL2.length;
    //printf("Reallocating");
    //AL1.layouts = (ArrayLayoutPos_t*)realloc(AL1.layouts, (AL1.length) * sizeof(ArrayLayoutPos_t));
    // Since we're updating AL1 in place, no need to copy entries before the overlap.

    // Merge positions where there is an overlap.
    for(int i = offset; i < AL1.length; i++){
        assert(i - offset < AL2.length);
        assert(i > 0);
        AL1.layouts[i] = cMergeLayoutPositions(AL1.layouts[i], AL2.layouts[i - offset]);
    }
    //LayoutOffset_t retValue;
    //retValue.Layout = AL1;
    //retValue.offset = offset;
    return offset;  // Note: Only merges the overlapping positions!
}
/*
ArrayLayout_t MergeLayouts(ArrayLayout_t AL1, ArrayLayout_t AL2){
    int start1, start2, offset;
    start1 = getFirstAlignedRefPos(AL1);
    start2 = getFirstAlignedRefPos(AL2);
    if(start1 > start2){
        // Switch the order to make AL1 the first read.
        ArrayLayout_t tmpLayout = AL2;
        AL2 = AL1;
        AL1 = tmpLayout;
        int tmpInt = start2;
        start2 = start1;
        start1 = tmpInt;
    }
    offset = start2 - start1;
    AL1.length = offset + AL2.length;
    AL1.layouts = (ArrayLayoutPos_t*)realloc(AL1.layouts, (AL1.length) * sizeof(ArrayLayoutPos_t));
    // Since we're updating AL1 in place, no need to copy entries before the overlap.

    // Merge positions where there is an overlap.
    for(int i = offset; i < AL1.length; i++){
        AL1.layouts[i] = cMergeLayoutPositions(AL1.layouts[i], AL2.layouts[i - offset]);
    }
    // Copy over the entries after the overlap.
    for(int i = AL1.length; i < AL1.length + offset; i++){
        AL1.layouts[i] = AL2.layouts[i - AL1.length + offset];
    }
    return AL1;
}
*/
