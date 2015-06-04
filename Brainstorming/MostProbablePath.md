## Fastq Recalibration
0. Prep
    1. Large, barcoded, deep family truth set for motif-causing errors
        1. Example: Phage Lambda will get us a lot of this context-specific error rate information, but not the whole larger space.
        2. Sesquicentamers!
1. Normalization:
    1. Run a barcoded PhiX each run
    2. Calculate Cycle by Cycle and Base by Base Error rate
    3. Re-calibrate each quality score for the sample.
2. Creating new reads:
    1. Slide down each read for context, using a context hashmap based on the output of (0.).
    2. Should probably have a probability for each potential nucleotide.
3. Demultiplexing families:
    1. bmftools dmp
4. Flattening read pairs.
3. Post-processing
    1. Cutadapt ("N" them all)
    2. Quality trimming? (e.g., Sickle?)
    3. Inverted repeats?
