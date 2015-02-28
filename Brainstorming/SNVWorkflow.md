#SNV Analysis Workflow

1. Demultiplexing families:
    1. Current use:
        1. FM: Number of reads in family
        2. FA: Number of reads in family agreeing on the consensus nucleotide for each base.
        3. FP: Pass/Fail for barcode processing.
        4. PV: Summed phred scores for flattened families.
    2. Questions:
        1. Keep or remove families of size 1? Currently they are kept but filtered out from inclusion in downstream calls.
        2. Keep the "summed" phred scores in the description? Compress to base 85?
        3. Should we remove from demultiplexing bases with a phred score under 10? [Messy in practice and probably marginal effect]
    3. Note: Can filter out read pairs here. What criteria?

2. Processing fastq records:
    1. cutadapt

3. Alignment:
    1. Currently, just bwa mem with primarily default settings.
    2. Cython wrapper of mem for indels and splicing.
    3. Further work to be done here? Realigning particular reads afterwards, perhaps?

3. Mark BAM:
    1. Note: Can filter out read pairs here. What criteria?
    2. Speed - can we write BAM records to file faster by batching them rather than writing one at a time?
    3. TODO: Add SV tags (and the new tags) to BAM tagging.
    4. Indel realignment (optimize abra parameters...)

4. Pileups:
    1. Control which reads are included in the call:
        1. Minimum MQ
        2. Minimum BQ:
            1. Minimum "summed" BQ
            2. Minimum phred-encoded BQ
        3. Extra tags:
            1. Include MDC, LI, ORB, ORU read pairs in SNV calls?
        4. Minimum FA for a base to be included in a call.
        5. Minimum FA fraction (FA / FM) for a base to be included in a call.

5. Actual call (INFO fields):
    1. Minimum mean MQ (redundant with above)
    2. Minimum mean BQ (redundant with above)
    3. Requiring a certain number of "duplex" reads.
    4. Max or min "mean" position in read or maximum SD. Should probably normalize SD to a coherent statistical metric.
    5. Strand Bias
    6. CSE-based notations.
    7. Check to see if reads 1 and 2 share significant sequence but map different places?
