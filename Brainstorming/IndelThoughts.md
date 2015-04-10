#Indel Calling
##"Local Reassembly"
1. For each window of size "w"
    1. Grab some number of kmers from the reference (not low-complexity, for example, perhaps, or just some at the start and finish)
    2. Align those kmers to the hg19 reference. Ignore those with MQ < minMQ. This will give us a net with which to pull down reads which are relevant to any portion of it and be flexible regarding reads which might have an indel or other event in the way preventing it from aligning to that section of the genome.
    3. Create a custom reference of those kmers padded with Ns. (make sure kmers used are not seriously redundant or else everything will get MQ 0)
    4. Align full fastq to said reference with bowtie (most flexible with permitting or not permitting mismatches). If it matches with some certain criteria being satisfied to any of these tiny kmers, that's fine. IME, bowtie is much, much faster than using grep, which we could also do.
    5. Find reads where the contig column is not "*" (Pipe bwa's output to awk '$3 != "*"')
    6. Feed those reads to fermi, velvet, or ???
    7. Make sure that multiple assemblies are permitted so that we can look at multiple genotypes. (I don't know enough about using them)
    8. Somehow reconnect those assemblies back to the region of interest. (Align straight? bwasw would likely be flexible enough to permit the necessary gap extensions. Old-fashioned BLAST could work, but it does use heuristics and accepts sub-optimal alignments.)
    9. Either align original sets of reads to this assembly or use that alignment to make a .vcf style record which could then be fed into a graph-based aligner.
    10. ???
    11. PROFIT
    12. Get off my lawn!
2. Thoughts on speed
    1. It's super fast to align to these mini-references. Maybe we won't take much of a hit by doing this.
    2. It would be an (almost) reference-agnostic local assembly-based indel caller. Bigger windows --> detecting bigger indels.
3. Broader applicability
    1. This could also be used for translocation detection and breakpoint reconstruction.
