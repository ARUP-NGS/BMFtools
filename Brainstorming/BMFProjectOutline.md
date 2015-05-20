#Barcoded Molecular Families Outline
##What can be separated?

1. Error Characterization/Correction:
    1. BMF Paper
        1. Family collapsing (complete)
            1. Purpose: SNVs and small indels, easing local assembly by correcting errors.
            2. Needed for accurate quantification of observed frequencies.
        2. Read pair collapsing (in progress)
            1. Purpose: SNVs and small indels
            2. Needed for accurate quantification of observed frequencies.
        3. Optical duplicate identification
            1. Needed for accurate quantification of observed frequencies.
        3. Polymerase error characterization
            1. Purpose: SNVs and small indels, easing local assembly by correcting errors.
            2. Look for bias - certain regions or sequences amplify better?
        4. Sequencer error characterization
            1. Purpose: SNVs and small indels, easing local assembly by correcting errors.
            2. Additionally, will permit us to lower sequencing power needed for specific tasks, guide experiments.

    2. Error Characterization/Correction project (FQSR)
        1. Use molecular barcodes with deep family coverage to get a set of "true" positives and negatives.
        2. Build a hash table/database listing error rates by contexts, as well as including extra information to be used by an SVM to look for contributions of error.
        3. With a given chemistry, a run output XML, and the error rates of the sequencer by cycle and by "correct" base.
        4. Process a fastq file to update quality score information and change a base if necessary.
        5. Ultimately, to improve the accuracy/sensitivity/specificity of the BMF methods.

2. Complex variants:

    1. Large indels (size >= 8 ?)
        1. Local assembly
            1. Tests on the HDx standards have not been successful with either scalpel or freebayes.
            2. Representative kmers/double hash map + local assembly.
                1. Fall back to an alignment should the kmers fail to be sufficiently unique markers.

    2. Structural variants (translocations)
        1. Mapping-based
            2. Subreads (e.g., first 30 bases + full reads)
        2. Assembly-based
            1. Try method in 2:1:2.
            2. Perhaps it would work with very large deletions, as well.

    3. Copy Number Aberrations/Alterations
        1. First: What kind of data do we need?
            1. Amplicon-based, or tons of probes with high on-target. Natera uses 23000 potentially heterozygous sites - I think we need a different assay.
        2. Combination of allelic imbalance and read depth.
        3. Note: This requires an accurate estimate of the tumor fraction.
            1. How do we go about that with a heterogeneous tumor population?
