#Error Profiling
##Metacode
### n is an integer in (4,12)
1. For each motif of size n:
    1. For each successive consensus nucleotide:
        1. Check for an indel (slide over f bases, where 1 >= f >=5?):
            1. If so, mark the error as that type, make counts and frequencies for it
            2. If not, proceed to next step
        2. Check for single nucleotide substitution calls. (Easier and quicker, though it's important to check for indels in order to make sure we're actually characterizing SNV errors separately from indels.

#CNA
##Rough outline
1. Pass 1:
    1. Establish background (counts/ratios) for capture with normal samples. 
    2. Establish variance/variability/noise
    3. Make putative calls solely based on counts.
2. Pass 2:
    1. Haplotype the sample (freebayes?) for potential het positions.
    2. Find allelic imbalances (part of this: bounds for normal allelic imbalances with our deep coverage.):
        1. Nifty for us to show if we're better at allelic imbalance discovery with our barcoding than people without. I assume it'll be a bit cleaner.
    3. If clearly outside the bounds of normal allelic imbalance/sequencing bias/etc. in such a way that it supports the putative CNA call, then mark it as an official call.
3. Pass 3:
    1. For regions where allelic imbalance was found and a raw count was not significantly different, I'd be concerned. How might you handle that?:
        1. Make a putative call, not unlike step 3 from pass 1.
