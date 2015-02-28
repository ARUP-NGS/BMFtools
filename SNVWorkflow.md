#SNV Analysis Workflow

1. Demultiplexing families
    a. Current use
        i. FM: Number of reads in family
        ii. FA: Number of reads in family agreeing on the consensus nucleotide for each base.
        iii. FP: Pass/Fail for barcode processing.
        iv. PV: Summed phred scores for flattened families.
    b. Questions
        i. Keep or remove families of size 1? Currently they are kept but filtered out from inclusion in downstream calls.
        ii. Keep the "summed" phred scores in the description? Compress to base 85?
        iii. Should we remove from demultiplexing bases with a phred score under 10? [Messy in practice and probably marginal effect]
    c. Note: Can filter out read pairs here. What criteria?

2. Processing fastq records
    a. cutadapt

3. Alignment
    a. Currently, just bwa mem with primarily default settings.
    b. Cython wrapper of mem for indels and splicing.
    c. Further work to be done here? Realigning particular reads afterwards, perhaps?

3. Mark BAM
    a. Note: Can filter out read pairs here. What criteria?
    b. Speed - can we write BAM records to file faster by batching them rather than writing one at a time?
    c. TODO: Add SV tags (and the new tags) to BAM tagging.
    d. Indel realignment (optimize abra parameters...)

4. Pileups
    a. Control which reads are included in the call.
        i. Minimum MQ
        ii. Minimum BQ
            1. Minimum "summed" BQ
            2. Minimum phred-encoded BQ
        iii. Extra tags
            1. Include MDC, LI, ORB, ORU read pairs in SNV calls?
        iv. Minimum FA for a base to be included in a call.
        v. Minimum FA fraction (FA / FM) for a base to be included in a call.

5. Actual call (INFO fields)
    a. Minimum mean MQ (redundant with above)
    b. Minimum mean BQ (redundant with above)
    c. Requiring a certain number of "duplex" reads.
    d. Max or min "mean" position in read or maximum SD. Should probably normalize SD to a coherent statistical metric.
    e. Strand Bias
    f. CSE-based notations.
    g. Check to see if reads 1 and 2 share significant sequence but map different places?
