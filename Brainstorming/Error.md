#Sequencing Error Characterization and Correction (MPS)
## Sources Of Error and Bias
1. Chemical
    1. PCR-induced errors
        1. Substitution (nucleotide).
        2. Indels.
        3. Template-switching
    2. Library Prep
        1. Ligation issues
        2. Other sources?
    3. Capture/Amplicon
        1. Capture
            1. Off-target (e.g., regions with homology issues)
            2. Bias
        2. Amplicon
            1. Mispriming
            2. Other sources?
    4. Correcting for strandedness of number of molecules in terms of de-duplicating.
2. Optical/Base calling
    1. Duplicates
        1. Barcodes eliminates the issues with PCR duplicates but not optical ones. We should be careful aboutthis.
    2. Substitutions (optical problems)
    3. Indels (optical problems)
3. Mapping errors
    1. False negatives/positives due to incorrect mapping location
        1. SNVs (I'm guessing primarily from a homology issue not properly accounted for?)
        2. Small indels (similar to SNVs, but aggravated in the context of tandem repeats)
            1. Shannon entropy filter available, but it is important to call in low complexity regions, too.
            2. Hopefully realignment (GATK/ABRA/bamleftalign/?) solves that.
        3. Large indels
            1. Read simply doesn't go to correct location.
                1. Mapping subreads should be able to help, but it looks like local assembly is the best way to look for this.
4. Additional variables
    1. Cluster size
    2. Fraction of clusters passing filters
    3. Number of clusters
5. Attributes to build SVM or other classifier on for fastq quality score and base recalibration
    1. "Correct" base at location.
    2. Cycle number of machine.
    3. Number of clusters.
    4. Fraction of clusters passing filters.
    5. Whether base is supported by the other read in the pair.
    6. Preceding sequence context.
    7. Perhaps succeeding context.
    8. Template switching/ concatemers.
