# TODO

## Software Development Housekeeping
1. Unit tests:
	1. MawCluster
		1. BCBam unit tests
		2. Pileup/SV/SNV unit tests?
			1. These modules need to be cleaned up a bit.
	2. utilBMF
		1. Finish QC unit tests.
		2. Identify essential HTSUtils unit tests.
		3. MergePairedAlignments unit test
			1. Contingent on the tool being finished.
2. PEP8 Compliance
	1. Write a script that runs pep8 on the BMFTools folder and ignores the errors permitted by our README.md file.
	2. Check the outputs of that script every ... how often?

## Analysis
1. Experimental Design:
    1. RESOLVED
	    1. Number of Ns needed for a given library diversity.
        2. We just decided "As many as we can get." Somewhere 12 <= N <= 16.
	2. What FM (Family Member size) do we need to naively get the specificity we want?
		1. We'll get this from the Phage data, hopefully?
	3. (After VC re-write) - what criteria do we need to make a correct call with the needed specificity?
		1. Combination of the [buccal/buffy] experiment and the tumor/normal cell line.
2. Analysis Framework
	1. RESOLVED
        1. Hellinger Distance - written.
	2. Compare R1/R2 Hellinger's.
	2. Test CycleByNucleotide script.
	

## Backend Software Development
1. Probabilistic Merging:
	1. Add read kickout for > 10% disagreement from consensus for fastq processing.
	2. Finish MPA
		1. With Unit test
2. "Variant Caller"
	1. Fix it up to skeleton level. (Good enough for tier 1 SNVs)
	2. Re-write it to be faster/better. [Stripped-down C struct flavor]
4. Non-SNV mutations
	1. Binning
	2. Assembly
	3. SW back
## 
