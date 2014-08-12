#!/mounts/anaconda/bin/python
import vcf

def MPileup(inputBAM,outputBCF="default",bedfile,ref):
    from subprocess import call;
    if(outputBCF=="default"):
        if(len(inputBAM.split('.')) >= 4):
            outputBCF = '.'.join(inputBAM.split('.')[0:-3]) + ".fullMP.vcf";
        else:
            outputBCF = '.'.join(inputBAM.split('.')[0:-1]) + ".fullMP.vcf";
    commandStr = "samtools mpileup -f {} -F 0.0001 -I -S -g -D -R -q 10 -Q 30 -l {} {} > {}".format(ref,bedfile,inputBAM,outputBCF)
    call(commandStr, shell=True)
    return outputBCF