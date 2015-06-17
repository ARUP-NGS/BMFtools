#include "MPA.h"
#include <stdlib.h>

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
			L1.agreement = 0;
			L1.quality = 0;
			L1.mergeAgreed = 0;
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
	ArrayLayout_t RetLayout;
	RetLayout.length = offset + AL2.length;
	RetLayout.layouts = (ArrayLayoutPos_t*)malloc((offset + AL2.length) * sizeof(ArrayLayoutPos_t));
	for(int i = 0; i < offset; i++){
		// Copy the entries before the overlap.
		RetLayout.layouts[i] = AL1.layouts[i];
	}
	for(int i = offset; i < AL1.length; i++){
		RetLayout.layouts[i] = cMergeLayoutPositions(AL1.layouts[i], AL2.layouts[i - offset]);
	}
	for(int i = AL1.length; i < AL1.length + offset; i++){
		RetLayout.layouts[i] = AL2.layouts[i - AL1.length + offset];
	}
	free(AL1.layouts);
	AL1 = RetLayout;
	return AL1;
}
