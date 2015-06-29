
'''
Unexpected procotols...

To get the number of non-reference mapped bases from a VCF, do these two things! :)

```bash
grep -v ^# $Infile| cut -f8 | grep -v 'INDEL' | cut -d';' -f2 | cut -d'=' -f2 | cut -d',' -f3,4 | sort | uniq -c
```
sum(int(k[0]) * sum(map(int, k[1].split(","))) for k in [i.strip().split() for i in open("NoBMF.nonref.hist", "r").read().split("\n") if i != ''] if sum(map(int, k[1].split(","))) < 100)
'''
