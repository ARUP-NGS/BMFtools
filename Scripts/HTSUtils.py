import logging


def align_bowtie2(R1, R2, ref, opts, outsam):
    import subprocess
    output = open(
        outsam,
        'w', 0)
    if(opts == ""):
        opts = '--threads 4 '
    if('--reorder' not in opts):
        opts += '--reorder '
    if('--mm' not in opts):
        opts += ' --mm '
    opt_concat = ' '.join(opts.split())
    command_str = ('bowtie2 {} --local --very'.format(opt_concat) +
                   '-sensitive-local -x {} -1 {} -2 {}'.format(ref, R1, R2))
    printlog(command_str)
    # command_list=command_str.split(' ')
    subprocess.call(command_str, stdout=output, shell=True)
    output.close()
    return(command_str)


def align_bwa(R1, R2, ref, opts, outsam):
    """Aligns a set of paired-end
    reads to a reference
    with provided options using bwa mem.
    Defaults to 4 threads, silent alignment, listing
    supplementary alignments, and
    writing each reads' alignment,
    regardless of mapping quality.
    """
    import subprocess
    if(opts == ""):
        opts = '-t 4 -v 1 -Y -T 0'
    output = open(outsam, 'w', 0)
    opt_concat = ' '.join(opts.split())

    command_str = 'bwa mem {} {} {} {}'.format(opt_concat, ref, R1, R2)
    # command_list = command_str.split(' ')
    printlog(command_str)
    subprocess.call(command_str, stdout=output, shell=True)
    output.close()
    return outsam, command_str


def align_bwa_se(reads, ref, opts, outsam):
    """Aligns a set of reads to a reference
    with provided options. Defaults to
    4 threads, silent alignment, listing
    supplementary alignments, and
    writing each reads' alignment,
    regardless of mapping quality.
    """
    import subprocess
    if(opts == ""):
        opts = '-t 4 -v 1 -Y -T 0'
    output = open(outsam, 'w', 0)
    opt_concat = ' '.join(opts.split())
    command_str = 'bwa mem {} {} {}'.format(opt_concat, ref, reads)
    # command_list = command_str.split(' ')
    printlog(command_str)
    subprocess.call(command_str, stdout=output, shell=True)
    output.close()
    return outsam, command_str


def align_snap(R1, R2, ref, opts, outbam):
    import subprocess
    opt_concat = ""
    for i, opt_it in enumerate(opts.split()):
        opt_concat += opt_it + " "
    command_str = "snap paired {} {} {} -o {} {}".format(
        ref,
        R1,
        R2,
        outbam,
        opt_concat)
    printlog(command_str)
    subprocess.call(command_str)
    return(command_str)


def CustomRefBowtiePaired(mergedFq, ref, output="default"):
    from subprocess import call
    if(output == "default"):
        output = mergedFq.split('.')[0] + '.mergingFamilies.sam'
    command_list = [
        "bowtie",
        "--threads",
        "4",
        "-S",
        "-n",
        "3",
        "-l",
        "24",
        "--norc",
        ref,
        mergedFq,
        output]
    call(command_list)
    return output


# This function is to handle StopIterations with a little elegance
def has_elements(iterable):
    from itertools import tee
    iterable, any_check = tee(iterable)
    try:
        any_check.next()
        return True, iterable
    except StopIteration:
        return False, iterable


def indexBowtie(fasta):
    from subprocess import call
    call('bowtie-build {0} {0}'.format(fasta), shell=True)
    return


def mergeBam(samList, memoryStr="-XmX16",
             MergeJar="/mounts/bin/picard-tools/MergeSamFiles.jar",
             outBam="default"):
    if(outBam == "default"):
        outBam = '.'.join(samList[0].split('.')[0:-1]) + '.merged.bam'
    from subprocess import call
    cStr = ("java -jar " + MergeJar + " " + memoryStr + " I=" +
            " I=".join(samList) + " O=" + outBam + " MSD=True " +
            "AS=True SO=coordinate"
            )
    printlog("About to merge bams. Command string: " + cStr)
    call(cStr, shell=True)
    return outBam


def sam_sort(insam, outsam):
    # Skip header and funnel all reads into a temp file
    import random
    output = open(outsam, 'w+')
    tmpname = 'temp{}.txt'.format(random.randint(0, 200))
    tmp = open(tmpname, 'w', 0)
    import subprocess
    command_str = str('grep -v "@SQ\|@PG\|VN:\|@HD" {}'.format(insam))
    printlog(command_str)
    subprocess.call(command_str, stdout=tmp, shell=True)
    tmp.close()
    # Save the header to the outsam
    command_header = 'grep "@SQ\|@PG\|@HD" {}'.format(insam)
    subprocess.call(command_header, stdout=output, shell=True)
    # sort the reads by query name
    tmp = open(tmpname, 'r')
    command_str1 = str('sort -k1,1 -t " " {}'.format(tmpname))
    printlog(command_str1)
    subprocess.call(command_str1, stdout=output, shell=True)
    output.close()
    tmp.close()
    subprocess.call('rm {}'.format(tmpname), shell=True)
    both_cmds = command_str + "\n" + command_str1
    return(both_cmds)


def printlog(string):
    print(string)
    logging.info(string)
    return


def parseConfigXML(string):
    import xmltodict
    orderedDict = xmltodict.parse(open(string, "r"))
    return orderedDict


def parseConfig(string):
    parsedLines = [l.strip() for l in open(string, "r").readlines()]
    ConfigDict = {}
    for line in parsedLines:
        ConfigDict[line.split("=").strip()[0]] = line.split("=")[1].strip()
    return ConfigDict

def FacePalm(string):
    str = ("............................................________ "
    "\n....................................,.-'\"...................``~. "
    "\n.............................,.-\"...................................\"-., "
    "\n.........................,/...............................................\":, "
    "\n.....................,?......................................................, "
    "\n.................../...........................................................,} "
    "\n................./......................................................,:`^`..} "
    "\n.............../...................................................,:\"........./ "
    "\n..............?.....__.........................................:`.........../ "
    "\n............./__.(.....\"~-,_..............................,:`........../ "
    "\n.........../(_....\"~,_........\"~,_....................,:`........_/ "
    "\n..........{.._$;_......\"=,_.......\"-,_.......,.-~-,},.~\";/....} "
    "\n...........((.....*~_.......\"=-._......\";,,./`..../\"............../ "
    "\n...,,,___.`~,......\"~.,....................`.....}............../ "
    "\n............(....`=-,,.......`........................(......;_,,-\" "
    "\n............/.`~,......`-...................................../ "
    "\n.............`~.*-,.....................................|,./.....,__ "
    "\n,,_..........}.>-._...................................|..............`=~-, "
    "\n.....`=~-,__......`,................................. "
    "\n...................`=~-,,.,............................... "
    "\n................................`:,,...........................`..............__ "
    "\n.....................................`=-,...................,%`>--==`` "
    "\n........................................_..........._,-%.......` "
    "\n...................................,")
    print(str)
    raise IllegalArgumentError(string)


class IllegalArgumentError(ValueError):
    pass
