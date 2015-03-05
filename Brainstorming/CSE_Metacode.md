#Error Profiling
##Metacode
### n is an integer in (4,12)
1. For each motif of size n:
    1. For each successive consensus nucleotide:
        1. Check for an indel (slide over f bases, where 1 >= f >=5?):
            1. If so, mark the error as that type, make counts and frequencies for it
            2. If not, proceed to next step
        2. Check for single nucleotide substitution calls. (Easier and quicker, though it's important to check for indels in order to make sure we're actually characterizing SNV errors separately from indels.

