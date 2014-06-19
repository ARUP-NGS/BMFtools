#!/home/daniel/anaconda/bin/python
# ####MapSplitReads
# ####Map Split Reads is a program designed to allow for trimming fastq files, aligning, and splicing them back in for the purpose of catching spanner reads for detecting chimeras, translocations, and gene fusions.
# ####Author: Daniel Baker
# ####Contact: daniel.baker@aruplab.com, d.nephi.baker@gmail.com
# ####License: Explicit permission from the author
# ####See README.md for further details.
def main():
    
    #Includes/Imports
    import os
    import time
    from Bio import SeqIO
    import pysam
    import argparse
    import subprocess
    from datetime import datetime as dt
    startTime=time.time()
    from random import randint

    #Parse Args
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--fastq-files', help="Provide your paired fastq files", nargs = 2, metavar = ('read1','read2'),required=True)
    parser.add_argument('-a','--aligner', help="Provide your aligner. E.g., 'bwa', 'bowtie2', or 'snap'. Currently only 'bwa' is supported. Default: bwa", nargs='?', metavar='aligner', default='bwa')
    parser.add_argument('-b','--BAM', help="BAM file, if alignment has already run.")
    parser.add_argument('-n','--trim-to', help="Provide the integer to which you would like to trim. Default: 30", nargs = '?', metavar = 'trim-to-N', default=30, type=int) 
    parser.add_argument('-r','--ref', help="Location of your reference indices", nargs = 1, metavar = 'reference',required=True)
    parser.add_argument('-o','--opts', help="Please place additional arguments after this tag. E.g.: --opts '-L 0' ", nargs='?', default='')
    parser.add_argument('-s','--sort', help="Whether the SAM file is to be sorted. Nota bene: if you are multithreading the aligner, you must select this option for it to go through.",default=True)
    parser.add_argument('--sv',help="Automatically select alignment settings with minimal structural variation penalty",default=False)
    args = parser.parse_args()
    
    #Set variables from arguments
    rds1=args.fastq_files[0]
    sv=args.sv
    cwd=os.getcwd()+'/'
    copy_command = "cp {} {}".format(os.path.abspath(rds1),cwd+os.path.basename(rds1))
    print(copy_command)
    subprocess.call(copy_command.split(' '))
    reads1=cwd+os.path.basename(rds1)
    print("reads 1 file is {}".format(reads1))
    del rds1
    rds2=args.fastq_files[1]
    copy_command = "cp {} {}".format(os.path.abspath(rds2),cwd+os.path.basename(rds2))
    print(copy_command)
    subprocess.call(copy_command.split(' '))
    reads2=cwd+os.path.basename(rds2)
    print("reads 2 file is {}".format(reads2))
    del rds2
    rf = args.ref[0]
    reference=os.path.abspath(rf)
    del rf
    extra_opts=args.opts
    basename=reads1.split('.')[0]
    trim = args.trim_to
    
    #Start Logging
    bufsize=0
    log_name = str(basename)+'-'+str(randint(0,1000)).zfill(4)+ '-' +str(dt.now().strftime("%m-%d-%H-%M")) + '.log'
    logger=open(log_name,'w+',bufsize)
    logger.write("Logger Started.\n")
    logger.write("It's a secret mission in uncharted space!\n")
    
    if(args.BAM!=None):
        AMfile=args.BAM
        logger.write('\nSkipping alignment. A BAM or SAM file was already provided.\n')
        print("Still sorting fastq's for compatibility.")
        reads1_sort=reads1.split('.')[0]+'_sort.fastq'
        reads2_sort=reads2.split('.')[0]+'_sort.fastq'
        fastq_command_1=fastq_sort(reads1,reads1_sort)
        logger.write("For Read 1 sorting, command was: ")
        logger.write(fastq_command_1+'\n')
        fastq_command_2=fastq_sort(reads2,reads2_sort)
        logger.write("For Read 2 sorting, command was: ")
        logger.write(fastq_command_2+'\n')

    else:
        #Sort the Fastq's
        logger.write("\nSorting the fastq's in preparation for alignment.\n")
        reads1_sort=reads1.split('.')[0]+'_sort.fastq'
        reads2_sort=reads2.split('.')[0]+'_sort.fastq'
       
        fastq_command_1=fastq_sort(reads1,reads1_sort)
        logger.write("For Read 1 sorting, command was: ")
        logger.write(fastq_command_1+'\n')
    
        fastq_command_2=fastq_sort(reads2,reads2_sort)
        logger.write("For Read 2 sorting, command was: ")
        logger.write(fastq_command_2+'\n')

        #Trim the fastq's
        logger.write("\nTrimming the fastq's in preparation for alignment.\n")
        reads1_trim=reads1_sort.split('.')[0]+'_N{}.fastq'.format(trim)
        reads2_trim=reads2_sort.split('.')[0]+'_N{}.fastq'.format(trim)
        basename=reads1_trim.split('.')[0]+'-' + args.aligner
        fastx_command1=fastx_trim(reads1_sort,reads1_trim,trim)
        fastx_command2=fastx_trim(reads2_sort,reads2_trim,trim)
        logger.write("For Read 1 trimming, command was: ")
        logger.write(' '.join(fastx_command1))
        logger.write("For Read 2 trimming, command was: ")
        logger.write(' '.join(fastx_command2))
        
        #Alignment step
        logger.write("\nBeginning the alignment portion of the process.\n")
        if(args.aligner=='bwa'):
            if(sv==True):
                extra_opts+=" -U 9 -t 8"
            bwa_command=align_bwa(reads1_trim,reads2_trim,reference,extra_opts,basename+'.sam')
            AMfile=basename+'.sam'
            print(bwa_command)
            logger.write(' '.join(bwa_command))
        elif(args.aligner=='snap'):
            snap_command=align_snap(reads1_trim,reads2_trim,reference,extra_opts,basename+'.bam')
            AMfile=basename+'.bam'
            print(snap_command)
            logger.write(' '.join(snap_command))
        elif(args.aligner=='bowtie2'):
            bowtie2_command=align_bowtie2(reads1_trim,reads2_trim,reference,extra_opts,basename+'.sam')
            AMfile=basename+'.sam'
            print(bowtie2_command)
            logger.write(' '.join(bowtie2_command))
        else:
            print('No other aligner supported!\n')
            return
        logger.write("Alignment completed with {}!\n".format(args.aligner))
        outsam=open(basename+'.sam')
        if(AMfile[-3:]=='bam'):
            subprocess.call(['samtools','view','-h'],stdout=outsam,shell=True)
        outsam.close()
        AMfile=basename+'.sam' 
    if(args.sort==True):
        sortAMfile=basename+'_sort.sam'
        sam_command=sam_sort(AMfile,sortAMfile)
        logger.write("For SAM sort: {} was the command\n".format(sam_command))
        print("For SAM sort: {} was the command\n".format(sam_command))
    else:
        sortAMfile=AMfile
        print('Sorting the SAM file was opted out of.')
        logger.write('Sorting the SAM file was opted out of.')
    #Parse in *amfile
    SamReads=pysam.Samfile(sortAMfile,'r')
 
    #Split *amfile into two, one for each end of a pair
    SamReads1=pysam.Samfile(basename+'_reads1.sam','wh',template=SamReads)
    SamReads2=pysam.Samfile(basename+'_reads2.sam','wh',template=SamReads)
    for samread in SamReads:
        if(samread.is_read1):
            SamReads1.write(samread)
        elif(samread.is_read2):
            SamReads2.write(samread)
        else:
            raise NameError("This does not seem to be paired-end.")
    SamReads.close()
    SamReads1.close()
    SamReads2.close()
   
    #Sort the samreads1 and samreads2 files
    '''
    sort_reads1=sam_sort(basename+'_reads1.sam',basename+'_reads1_sort.sam')
    print(sort_reads1)
    logger.write(sort_reads1+'was the command for sorting reads 1')
    subprocess.call('cp {} {}_reads1.sam'.format(basename+'_reads1_sort.sam',basename),shell=True)
    sort_reads2=sam_sort(basename+'_reads2.sam',basename+'_reads2_sort.sam')
    subprocess.call('cp {} {}_reads2.sam'.format(basename+'_reads2_sort.sam',basename),shell=True)
    print(sort_reads2)
    logger.write(sort_reads2+'was the command for sorting reads 2')
    '''
    #Parse in fastq files
    fr1 = SeqIO.parse(reads1_sort,'fastq')
    fr2 = SeqIO.parse(reads2_sort,'fastq')
    #Loop through Reads1 filsame
    sam = pysam.Samfile(AMfile,'r')
    outsam = pysam.Samfile(basename+'-PE1-tf.sam','wh',template=sam)
    tempname = open(basename+'tempnames.txt','w+',0)
    
    #Fix pair 1's 
    SamReads1=pysam.Samfile(basename+'_reads1.sam','r',until_eof=True)
    for samline in SamReads1:
        WorkSam=samline
        WorkFq=fr1.next()
        if(WorkSam.qname not in WorkFq.name):
            tempname.write("Sam name = " +WorkSam.qname + " and fq name = " + WorkFq.name)
            raise NameError("Warning - SAM and fastq files are not in the same order")
        qual_holder=WorkFq.letter_annotations['phred_quality']
        qual_string=''
        for qscore in qual_holder:
            qual_string+=str(chr(qscore+33))
        fullseq_holder=str(WorkFq.seq)
        seq_len=len(fullseq_holder)
        ttrim=len(WorkSam.seq)
        WorkSam.seq=fullseq_holder
        #Replace qualities
        if(WorkSam.is_unmapped):
            logger.write('Unmapped read: {}\n'.format(WorkSam.qname))
            WorkSam.qual=qual_string
            outsam.write(WorkSam)
            continue
        if(WorkSam.cigar[-1][0] != 4):
            WorkSam.cigar+=[(4,seq_len-ttrim)]
        else:
            S_hold=seq_len-ttrim+WorkSam.cigar[-1][1]
            WorkSam.cigar=WorkSam.cigar[:-1]
            WorkSam.cigar+=[(4,S_hold)]
        WorkSam.qual=qual_string
        outsam.write(WorkSam)
    outsam.close()
    #Fix pair 2's
    outsam2 = pysam.Samfile(basename+'-PE2-tf.sam','wh',template=sam)
    tempname = open(basename+'tempnames.txt','w+',0)
    SamReads2=pysam.Samfile(basename+'_reads2.sam','r',until_eof=True)
    for samline in SamReads2:
        WorkSam=samline
        WorkFq=fr2.next()
        if(WorkSam.qname not in WorkFq.name):
            tempname.write("Sam name = " +WorkSam.qname + " and fq name = " + WorkFq.name+'\n')
            raise NameError("Warning - SAM and fastq files are not in the same order")
        qual_holder=WorkFq.letter_annotations['phred_quality']
        qual_string=''
        for qscore in qual_holder:
            qual_string+=str(chr(qscore+33))
        fullseq_holder=str(WorkFq.seq)
        seq_len=len(fullseq_holder)
        ttrim=len(WorkSam.seq)
        WorkSam.seq=fullseq_holder
        #Replace qualities
        if(WorkSam.is_unmapped):
            logger.write('Unmapped read: {}\n'.format(WorkSam.qname))
            WorkSam.qual=qual_string
            outsam.write(WorkSam)
            continue
        if(WorkSam.cigar[-1][0] != 4):
            WorkSam.cigar+=[(4,seq_len-ttrim)]
        else:
            S_hold=seq_len-ttrim+WorkSam.cigar[-1][1]
            WorkSam.cigar=WorkSam.cigar[:-1]
            WorkSam.cigar+=[(4,S_hold)]
        WorkSam.qual=qual_string
        outsam2.write(WorkSam)
    outsam2.close()
    tempname.close()
    #If tempname (the out of order pairs) is empty, delete it.
    if(os.stat(basename+'tempnames.txt')[6]==0):
        subprocess.call(['rm',basename+'tempnames.txt'])
    #Clean up
    #Comment this code if you want to save the trimmed fastq files or other working files
    subprocess.call(['rm',reads1_sort])
    subprocess.call(['rm',reads2_sort])
    subprocess.call(['rm',basename+'_reads1.sam'])
    subprocess.call(['rm',basename+'_reads2.sam'])
    #Convert each samfile into a BAM file for merging
    sam1 = basename+'-PE1-tf.sam'.format(trim,args.aligner)
    sam2 = basename+'-PE2-tf.sam'.format(trim,args.aligner)
    bam1 = basename+'.{}.{}.final-1.bam'.format(trim,args.aligner)
    bam2 = basename+'.{}.{}.final-2.bam'.format(trim,args.aligner)
    
    #Convert SAMs to BAM for merging, then remove the SAM files if the BAM files aren't empty
    bam1Command=Sam2Bam(sam1,bam1)
    print(bam1Command)
    #logger.write('\nSam2Bam command is: {}\n'.format(' '.join(bam1Command)))
    #subprocess.call(bam1Command)
    if(os.stat(bam1)[6]!=0):
        subprocess.call(['rm',sam1])

    bam2Command=Sam2Bam(sam2,bam2)
    print(bam2Command)
    #logger.write('\nSam2Bam command is: {}\n'.format(' '.join(bam2Command)))
    #subprocess.call(bam2Command)
    if(os.stat(bam2)[6]!=0):
        subprocess.call(['rm',sam2])
    
    #Merge the two BAM files
    final_bam=basename+'.{}.final.bam'.format(trim)
    final_str='samtools merge -n -h {} {} {} {}'.format(sam,final_bam,bam1,bam2)
    subprocess.call(final_str,shell=True)
    sam.close()
    #If things go wrong, you can comment this block for debugging
    subprocess.call(['rm',AMfile])
    subprocess.call(['rm',bam1])
    subprocess.call(['rm',bam2])
    subprocess.call(['rm',sam])
    
    logger.write("Realignment completed\n")
    print("Total time taken is {}s.\n".format(str(time.time()-startTime)))

    logger.write("Total time taken is {}s.\n".format(str(time.time()-startTime)))
    print("You are awesome. Congratulations.\n")

    logger.close()
    #Fin!

def sam_sort(insam,outsam):
    #Skip header and funnel all reads into a temp file
    import random
    output=open(outsam,'w+')
    tmpname='temp{}.txt'.format(random.randint(0,200))
    tmp=open(tmpname,'w',0)
    import subprocess
    command_str=str('grep -v "@SQ\|@PG\|VN:\|@HD" {}'.format(insam))
    print(command_str)
    subprocess.call(command_str,stdout=tmp,shell=True)
    tmp.close()
    #Save the header to the outsam
    command_header = 'grep "@SQ\|@PG\|@HD" {}'.format(insam)
    subprocess.call(command_header, stdout=output, shell=True)
    #sort the reads by query name
    tmp=open(tmpname,'r')
    command_str1=str('sort -k1,1 -t " " {}'.format(tmpname))
    print(command_str1)
    subprocess.call(command_str1,stdout=output,shell=True)
    output.close()
    tmp.close()
    subprocess.call('rm {}'.format(tmpname),shell=True)
    both_cmds=command_str+"\n"+command_str1
    return(both_cmds)

def Sam2Bam(insam,outbam):
    from subprocess import call
    output = open(outbam,'w',0)
    command_str='samtools view -Sb {}'.format(insam,shell=True)
    print(command_str)
    call(command_str,stdout=output,shell=True)
    return(command_str,outbam)

def Bam2Sam(inbam,outsam):
    from subprocess import call
    output = open(outsam,'w',0)
    command_str='samtools view -h {}'.format(inbam)
    print(command_str)
    call(command_str,stdout=output,shell=True)
    return(command_str,outsam)

def fastx_trim(infq, outfq, n):
    import subprocess
    command_str = ['fastx_trimmer','-l',str(n),'-i',infq,'-o',outfq]
    print(command_str)
    subprocess.call(command_str)
    return(command_str)

def align_bwa(R1,R2,ref,opts,outsam):
    import subprocess
    opt_concat = ""
    if(opts== ""):
        opts='-t 4 -T 0'
    output = open(outsam,'w',0)
    for i, opt_it in enumerate(opts.split()):
        opt_concat+=opt_it+" "
    command_str = 'bwa mem {} {} {} {}'.format(opt_concat,ref,R1,R2)
    #command_list = command_str.split(' ')
    print(command_str)
    subprocess.call(command_str, stdout=output,shell=True)
    output.close()
    return command_str;

def align_snap(R1,R2,ref,opts,outbam):
    import subprocess
    opt_concat = ""
    for i, opt_it in enumerate(opts.split()):
        opt_concat+=opt_it+" "
    command_str = "snap paired {} {} {} -o {} {}".format(ref,R1,R2,outbam,opt_concat)
    print(command_str)
    subprocess.call(command_str)
    return(command_str)

def align_bowtie2(R1,R2,ref,opts,outsam):
    import subprocess
    opt_concat=""
    output=open(outsam,'w',0) #Zero buffer size means it gets sent immediately straight to the file
    if(opts==""):
        opts='--threads 4 '
    if('--reorder' not in opts):
        opts+='--reorder '
    if('--mm' not in opts):
        opts+=' --mm '
    for i, opt_it in enumerate(opts.split()):
        opt_concat+=opt_it+" "
    command_str = 'bowtie2 {} --local --very-sensitive-local -x {} -1 {} -2 {}'.format(opt_concat,ref,R1,R2)
    print(command_str)
    #command_list=command_str.split(' ')
    subprocess.call(command_str, stdout=output,shell=True)
    output.close()
    return(command_str)

def fastq_sort(in_fastq,out_fastq):
    import subprocess
    outfile=open(out_fastq,'w')
    command_str='cat {} | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n"'.format(in_fastq)
    subprocess.call(command_str,stdout=outfile,shell=True)
    outfile.close()
    return(command_str)

if __name__ == "__main__":
    main()
