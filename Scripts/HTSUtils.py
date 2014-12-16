import logging
import shlex
import subprocess


class Configurations:
    """
    Holds both the config '=' dictionary
    and the config xml orderedDict
    """
    def __init__(self, configFile, configXml):
        self.config = parseConfig(configFile)
        self.configOrderedDict = parseConfigXML(configXml)

    def GetValueForKey(self, key, multiple=False):
        if(key in self.config.keys()):
            return self.config[key]
        else:
            return get_recursively(
                self.configOrderedDict, key, multiple=multiple)


def align_bowtie2(R1, R2, ref, opts, outsam):
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
    subprocess.check_call(command_str, stdout=output, shell=True)
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
    if(opts == ""):
        opts = '-t 4 -v 1 -Y -T 0 -U 5'
    output = open(outsam, 'w', 0)
    opt_concat = ' '.join(opts.split())

    command_str = 'bwa mem {} {} {} {}'.format(opt_concat, ref, R1, R2)
    # command_list = command_str.split(' ')
    printlog(command_str)
    subprocess.check_call(shlex.split(command_str), stdout=output)
    output.close()
    return outsam


def align_bwa_se(reads, ref, opts, outsam):
    """Aligns a set of reads to a reference
    with provided options. Defaults to
    4 threads, silent alignment, listing
    supplementary alignments, and
    writing each reads' alignment,
    regardless of mapping quality.
    """
    if(opts == ""):
        opts = '-t 4 -v 1 -Y -T 0'
    output = open(outsam, 'w', 0)
    opt_concat = ' '.join(opts.split())
    command_str = 'bwa mem {} {} {}'.format(opt_concat, ref, reads)
    # command_list = command_str.split(' ')
    subprocess.check_call(shlex.split(command_str), stdout=output, shell=False)
    output.close()
    return outsam, command_str


def align_snap(R1, R2, ref, opts, outbam):
    opt_concat = " ".join(opts.split())
    command_str = "snap paired {} {} {} -o {} {}".format(
        ref,
        R1,
        R2,
        outbam,
        opt_concat)
    printlog(command_str)
    subprocess.check_call(shlex.split(command_str), shell=False)
    return(command_str)


def PipedShellCall(commandStr):
    import uuid
    PipedShellCallFilename = "PipedShellCall{}.sh".format(
        str(uuid.uuid4().get_hex().upper()[0:8]))
    printlog("Command string: {}".format(commandStr))
    open(PipedShellCallFilename, "w").write(commandStr)
    subprocess.check_call(['bash', PipedShellCallFilename])
    return commandStr


def CustomRefBowtiePaired(mergedFq,
                          ref,
                          output="default",
                          barLen="default",
                          bowtiePath="bowtie",
                          mismatchLimit="default"):
    if(output == "default"):
        output = mergedFq.split('.')[0] + '.mergingFamilies.sam'
    if(barLen == "default"):
        FacePalm("Barcode length must be set. Abort mission!")
    if(bowtiePath == "bowtie"):
        printlog("Defaulting to bowtie for path to executable.")
    elif("2" in bowtiePath.split("/")[-1]):
        FacePalm("Do not use bowtie2!")
    command_list = [
        "bowtie",
        "--threads",
        "4",
        "-S",
        "-n",
        str(mismatchLimit),
        "-l",
        str(barLen),
        "--norc",
        ref,
        mergedFq,
        output]
    subprocess.check_call(command_list)
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
    subprocess.check_call('bowtie-build {0} {0}'.format(fasta), shell=True)
    return


def mergeBam(samList, memoryStr="-XmX16",
             MergeJar="/mounts/bin/picard-tools/MergeSamFiles.jar",
             outBam="default"):
    if(outBam == "default"):
        outBam = '.'.join(samList[0].split('.')[0:-1]) + '.merged.bam'
    cStr = ("java -jar " + MergeJar + " " + memoryStr + " I=" +
            " I=".join(samList) + " O=" + outBam + " MSD=True " +
            "AS=True SO=coordinate"
            )
    printlog("About to merge bams. Command string: " + cStr)
    subprocess.check_call(shlex.split(cStr))
    return outBam


def sam_sort(insam, outsam):
    # Skip header and funnel all reads into a temp file
    import random
    output = open(outsam, 'w+')
    tmpname = 'temp{}.txt'.format(random.randint(0, 200))
    tmp = open(tmpname, 'w', 0)
    command_str = str('grep -v "@SQ\|@PG\|VN:\|@HD" {}'.format(insam))
    printlog(command_str)
    subprocess.check_call(shlex.split(command_str), stdout=tmp)
    tmp.close()
    # Save the header to the outsam
    command_header = 'grep "@SQ\|@PG\|@HD" {}'.format(insam)
    subprocess.check_call(shlex.split(command_header), stdout=output)
    # sort the reads by query name
    tmp = open(tmpname, 'r')
    command_str1 = str('sort -k1,1 -t " " {}'.format(tmpname))
    printlog(command_str1)
    subprocess.check_call(shlex.split(command_str1), stdout=output)
    output.close()
    tmp.close()
    subprocess.check_call(['rm', tmpname])
    printlog("Command 1 for sam sort: " + command_str)
    printlog("Command 2 for sam sort: " + command_str1)
    return(outsam)


def get_recursively(search_dict, field, multiple=False):
    """
    Takes a dict with nested lists and dicts,
    and searches all dicts for a key of the field
    provided.
    If multiple is set, it returns a list.
    Otherwise, the first found value is returned.
    """
    fields_found = []

    for key, value in search_dict.iteritems():

        if key == field:
            if(multiple is False):
                return value
            fields_found.append(value)

        elif isinstance(value, dict):
            results = get_recursively(value, field)
            for result in results:
                fields_found.append(result)

        elif isinstance(value, list):
            for item in value:
                if isinstance(item, dict):
                    more_results = get_recursively(item, field)
                    for another_result in more_results:
                        fields_found.append(another_result)

    return fields_found


def printlog(string):
    print(string.replace(
        "\t", "\\t").replace("\n", "\\n"))
    logging.info(string.replace(
        "\t", "\\t").replace("\n", "\\n"))
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
    Str = ("............................................________ "
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
    print(Str)
    raise IllegalArgumentError(string)


class IllegalArgumentError(ValueError):
    pass


class ThisIsMadness(Exception):
    pass
