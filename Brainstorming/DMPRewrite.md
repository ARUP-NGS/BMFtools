## Approaches
1. fqmarksplit
  1. Description
    1. Select a prefix length. (1-4)
    2. Open 4**(Prefix length) handles
    3. Split the fastq while adding the barcode sequence to the comment field, writing the record (on a single line) to temporary files or mmap'd files.
    4. Call Heng Li's fork of GNU sort within C, piped (within C) to the meta-analysis with the CEPHES astrophysics library's igamcll on each base call with a #pragma omp parallel for after preparing the file handles for input/output, whether all in RAM as mmap or having it actually be written to disk.
    5. Apply the CEPHES astrophysics library's igamcll to each base call with the appropriate C code (see BCFastq.pyx for an example).
    6. Bring all of the demultiplexed sets back to a final set of fastqs.
  2. Strengths
    1. Makes parallelizing the sort and DMP easier.
  3. Weaknesses
    1. I/O, unless written to mmap'd files. (Ideally, should be written so that the DMP can take streaming input, and have the sort write to stdout)
    2. Merging the final files back together - unless the right answer is to make a shell call to bwa which takes <(cat ${FileListR1}) <(cat ${FileListR2}) for R1/R2.
2. Read many, write once.
  1. Description
    1. Pick a prefix length.
    2. Make a list of each prefix, then send off workers which read the whole file but skip any reads that don't have its prefix.
    3. Place those reads into a hashmap, then demultiplex those sets of records, writing to a merged fastq file.
    4. Bring all of the demultiplexed sets back to a final set of fastqs.
  2. Strengths
    1. Writing to disk: only writes the final, demultiplexed fastqs.
    2. RAM efficient - since only these small fractions of the input fastq are worked with at a given point of time, there isn't one gigantic hashmap that has to fit into RAM.
  3. Weaknesses
    1. Reading from disk: reads the file 4**(Prefix length) times rather than one.
3. Post-alignment
  1. Description
    1. Align reads after adding barcode/pass fail to the fastq comment field.
    2. Use bwa -C to put these comments into the BAM.
    3. Create a tmp_stack_t object (see bam_rmdup.c). Move through the stack, looking for reads within some mismatch limit of the barcode present.
    4. Write the flattened reads to an output stream.
  2. Strengths
    1. Reads with the same/similar barcode have been grouped together somewhat. Since sort is already needed for BAM random access, a lot of the sort overhead is dealt with.
  3. weaknesses
    1. BWA overhead - if family sizes are large, we're aligning 20x as many reads.
    2. Computational complexity of handling the barcode set comparisons, especially as read depth increases.
