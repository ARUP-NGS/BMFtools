import logging
import shlex
import subprocess
import pybedtools
from pybedtools import BedTool as BT


class Configurations:
    """
    Holds the config json dict
    """
    def __init__(self, configJSON):
        self.config = parseConfigJSON(configJSON)

    def GetValueForKey(self, key, multiple=False):
        if(key in self.config.keys()):
            return self.config[key]
        else:
            return get_recursively(
                self.config, key, multiple=multiple)


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
        opts = '-t 4 -v 1 -Y -M -T 0'
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


def PipedShellCall(commandStr, delete=True):
    import uuid
    PipedShellCallFilename = "PipedShellCall{}.sh".format(
        str(uuid.uuid4().get_hex().upper()[0:8]))
    printlog("Command string: {}".format(commandStr))
    open(PipedShellCallFilename, "w").write(commandStr)
    subprocess.check_call(['bash', PipedShellCallFilename])
    if(delete is True):
        subprocess.check_call(shlex.split(
            "rm {}".format(PipedShellCallFilename)))
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


def CoorSortAndIndexBam(inbam, prefix="MetasyntacticVar",
                        outbam="default",
                        uuid="true",
                        threads="4"):
    # If uuid is either a boolean true or is a string containing true,
    # then a random string is generated for the output
    if(str(uuid).lower() == "true"):
        import uuid
        prefix += str(uuid.uuid4().get_hex().upper()[0:8])
    if(outbam == "default"):
        outbam = '.'.join(inbam.split('.')[0:-1]) + '.CoorSort.bam'
    CommandStr = ("samtools sort -T {} -O bam -o {}".format(prefix, outbam) +
                  " -@ {} {}".format(threads, inbam))
    printlog("About to call sort command: {}".format(CommandStr))
    subprocess.check_call(shlex.split(CommandStr))
    printlog("Now indexing.")
    subprocess.check_call(shlex.split("samtools index {}".format(outbam)))
    return outbam


def NameSort(inbam, outbam="default", prefix="MetasyntacticVar",
             uuid="true", threads="4"):
    # If uuid is either a boolean true or is a string containing true,
    # then a random string is generated for the output
    if(str(uuid).lower() == "true"):
        import uuid
        prefix += str(uuid.uuid4().get_hex().upper()[0:8])
    if(outbam == "default"):
        outbam = '.'.join(inbam.split('.')[0:-1]) + '.NameSort.bam'
    CommandStr = ("samtools sort -T {} -O bam -o {}".format(prefix, outbam) +
                  " -@ {} -n {}".format(threads, inbam))
    printlog("About to call sort command: {}".format(CommandStr))
    subprocess.check_call(shlex.split(CommandStr))
    printlog("Namesort successful, sorted bam available at: {}".format(outbam))
    return outbam


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


def GetBamTagsDictionary(samRecord):
    TagsDict = {}
    for pair in samRecord.tags:
        TagsDict[pair[0]] = pair[1]
    return TagsDict


def GetBamTag(samRecord, tag):
    try:
        return GetBamTagsDictionary(samRecord)[tag]
    except KeyError:
        raise ThisIsMadness("Alignment record missing requested tag.")


def ReadContainedInBed(samRecord, bedRef="default"):
    """
    Checks to see if a samRecord is contained in a bedfile.
    """
    if(isinstance(bedRef, BT) is False):
        raise ThisIsMadness("bedRef must be a pybedtools bedtool object.")
    ersatzBF = BT(PysamToChrDict[samRecord.reference_id],
                  str(samRecord.reference_start),
                  str(samRecord.reference_end),
                  from_string=True)
    if(len(ersatzBF.intersect(bedRef)) > 0):
        return True
    else:
        return False


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


def printlog(string, level=logging.INFO):
    Logger = logging.getLogger("Primarylogger")
    if(level == logging.DEBUG):
        Logger.debug(string.replace(
            "\t", "\\t").replace("\n", "\\n"))
    elif(level == logging.INFO):
        Logger.info(string.replace(
            "\t", "\\t").replace("\n", "\\n"))
    elif(level == logging.WARNING):
        Logger.warning(string.replace(
            "\t", "\\t").replace("\n", "\\n"))
    else:
        Logger.critical(string.replace(
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


def parseConfigJSON(path):
    import json
    json_file = open(path)
    json_str = json_file.read()
    json_data = json.loads(json_str)
    return json_data


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


"""
This dictionary is required by some tools in the package to deal with the fact
that pysam doesn't store contig names, only which contig number it is,
based on the order found in the sam header.
"""

PysamToChrDict = {}
for i in range(22):
    PysamToChrDict[i] = str(i + 1)
PysamToChrDict[22] = "X"
PysamToChrDict[23] = "Y"
PysamToChrDict[24] = "MT"
PysamToChrDict[25] = "GL000207.1"
PysamToChrDict[26] = "GL000226.1"
PysamToChrDict[27] = "GL000229.1"
PysamToChrDict[28] = "GL000231.1"
PysamToChrDict[29] = "GL000210.1"
PysamToChrDict[30] = "GL000239.1"
PysamToChrDict[31] = "GL000235.1"
PysamToChrDict[32] = "GL000201.1"
PysamToChrDict[33] = "GL000247.1"
PysamToChrDict[34] = "GL000245.1"
PysamToChrDict[35] = "GL000197.1"
PysamToChrDict[36] = "GL000203.1"
PysamToChrDict[37] = "GL000246.1"
PysamToChrDict[38] = "GL000249.1"
PysamToChrDict[39] = "GL000196.1"
PysamToChrDict[40] = "GL000248.1"
PysamToChrDict[41] = "GL000244.1"
PysamToChrDict[42] = "GL000238.1"
PysamToChrDict[43] = "GL000202.1"
PysamToChrDict[44] = "GL000234.1"
PysamToChrDict[45] = "GL000232.1"
PysamToChrDict[46] = "GL000206.1"
PysamToChrDict[47] = "GL000240.1"
PysamToChrDict[48] = "GL000236.1"
PysamToChrDict[49] = "GL000241.1"
PysamToChrDict[50] = "GL000243.1"
PysamToChrDict[51] = "GL000242.1"
PysamToChrDict[52] = "GL000230.1"
PysamToChrDict[53] = "GL000237.1"
PysamToChrDict[54] = "GL000233.1"
PysamToChrDict[55] = "GL000204.1"
PysamToChrDict[56] = "GL000198.1"
PysamToChrDict[57] = "GL000208.1"
PysamToChrDict[58] = "GL000191.1"
PysamToChrDict[59] = "GL000227.1"
PysamToChrDict[60] = "GL000228.1"
PysamToChrDict[61] = "GL000214.1"
PysamToChrDict[62] = "GL000221.1"
PysamToChrDict[63] = "GL000209.1"
PysamToChrDict[64] = "GL000218.1"
PysamToChrDict[65] = "GL000220.1"
PysamToChrDict[66] = "GL000213.1"
PysamToChrDict[67] = "GL000211.1"
PysamToChrDict[68] = "GL000199.1"
PysamToChrDict[69] = "GL000217.1"
PysamToChrDict[70] = "GL000216.1"
PysamToChrDict[71] = "GL000215.1"
PysamToChrDict[72] = "GL000205.1"
PysamToChrDict[73] = "GL000219.1"
PysamToChrDict[74] = "GL000224.1"
PysamToChrDict[75] = "GL000223.1"
PysamToChrDict[76] = "GL000195.1"
PysamToChrDict[77] = "GL000212.1"
PysamToChrDict[78] = "GL000222.1"
PysamToChrDict[79] = "GL000200.1"
PysamToChrDict[80] = "GL000193.1"

PysamToChrDict[82] = "GL000225.1"
PysamToChrDict[83] = "GL000192.1"

ChrToPysamDict = {}

for i in range(22):
    PysamToChrDict[str(i + 1)] = i
ChrToPysamDict["GL000207.1"] = "25"
ChrToPysamDict["GL000226.1"] = "26"
ChrToPysamDict["GL000229.1"] = "27"
ChrToPysamDict["GL000231.1"] = "28"
ChrToPysamDict["GL000210.1"] = "29"
ChrToPysamDict["GL000239.1"] = "30"
ChrToPysamDict["GL000235.1"] = "31"
ChrToPysamDict["GL000201.1"] = "32"
ChrToPysamDict["GL000247.1"] = "33"
ChrToPysamDict["GL000245.1"] = "34"
ChrToPysamDict["GL000197.1"] = "35"
ChrToPysamDict["GL000203.1"] = "36"
ChrToPysamDict["GL000246.1"] = "37"
ChrToPysamDict["GL000249.1"] = "38"
ChrToPysamDict["GL000196.1"] = "39"
ChrToPysamDict["GL000248.1"] = "40"
ChrToPysamDict["GL000244.1"] = "41"
ChrToPysamDict["GL000238.1"] = "42"
ChrToPysamDict["GL000202.1"] = "43"
ChrToPysamDict["GL000234.1"] = "44"
ChrToPysamDict["GL000232.1"] = "45"
ChrToPysamDict["GL000206.1"] = "46"
ChrToPysamDict["GL000240.1"] = "47"
ChrToPysamDict["GL000236.1"] = "48"
ChrToPysamDict["GL000241.1"] = "49"
ChrToPysamDict["GL000243.1"] = "50"
ChrToPysamDict["GL000242.1"] = "51"
ChrToPysamDict["GL000230.1"] = "52"
ChrToPysamDict["GL000237.1"] = "53"
ChrToPysamDict["GL000233.1"] = "54"
ChrToPysamDict["GL000204.1"] = "55"
ChrToPysamDict["GL000198.1"] = "56"
ChrToPysamDict["GL000208.1"] = "57"
ChrToPysamDict["GL000191.1"] = "58"
ChrToPysamDict["GL000227.1"] = "59"
ChrToPysamDict["GL000228.1"] = "60"
ChrToPysamDict["GL000214.1"] = "61"
ChrToPysamDict["GL000221.1"] = "62"
ChrToPysamDict["GL000209.1"] = "63"
ChrToPysamDict["GL000218.1"] = "64"
ChrToPysamDict["GL000220.1"] = "65"
ChrToPysamDict["GL000213.1"] = "66"
ChrToPysamDict["GL000211.1"] = "67"
ChrToPysamDict["GL000199.1"] = "68"
ChrToPysamDict["GL000217.1"] = "69"
ChrToPysamDict["GL000216.1"] = "70"
ChrToPysamDict["GL000215.1"] = "71"
ChrToPysamDict["GL000205.1"] = "72"
ChrToPysamDict["GL000219.1"] = "73"
ChrToPysamDict["GL000224.1"] = "74"
ChrToPysamDict["GL000223.1"] = "75"
ChrToPysamDict["GL000195.1"] = "76"
ChrToPysamDict["GL000212.1"] = "77"
ChrToPysamDict["GL000222.1"] = "78"
ChrToPysamDict["GL000200.1"] = "79"
ChrToPysamDict["GL000193.1"] = "80"
ChrToPysamDict["GL000194.1"] = "81"
ChrToPysamDict["GL000225.1"] = "82"
ChrToPysamDict["GL000192.1"] = "83"
