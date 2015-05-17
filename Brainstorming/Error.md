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
3. Additional variables
    1. Cluster size
    2. Fraction of clusters passing filters
    3. Number of clusters
4. Attributes to build SVM or other classifier on
    1. "Correct" base at location.
    2. Cycle number of machine.
    3. Number of clusters.
    4. Fraction of clusters passing filters.
    5. Whether base is supported by the other read in the pair.
    6. Preceding sequence context.
    7. Perhaps succeeding context.
    8. Template switching/ concatemers.
