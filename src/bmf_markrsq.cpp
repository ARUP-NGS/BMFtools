/*
 * OUTLINE
 * 1. Mark/prepare for sorting.
 *   1. Handle supplementary/secondary
 *     1. How?
 *       1. Make a stack for reads with a given read name
 *         1. Add tags for SU/MU/LM to all reads in the set such that they have the same keys
 *         2. Add a tag to read 1 and read 2 for each (if any) of its supplemental alignments to know what to wait for.
 *           1. Add a tag to the supplementals for their actual position/contig, flag as unmapped, and move to their primaries.
 *         3. Output to stdout
 * 2. Pass to samtools sort and sort by coordinate.
 * 3. Load in a buffer of reads
 *   1. Fill a large stack of buffered reads.
 *   2. Build a hashmap for r1/r2 combinations.
 *   3. Once all reads needed for an alignment signature set have been loaded, collapse them, putting the supplementals in a separate table.
 *     1. If a read name set is not collapsed and there are supplementals, unset the unmapped flag and change the coordinates back to what they should be.
 *     2. Otherwise, ignore the supplementals because they'll be realigned.
 */
#include "bmf_markrsq.h"

namespace BMF {

}
