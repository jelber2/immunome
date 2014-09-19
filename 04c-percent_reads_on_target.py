#! /usr/bin/env python

# PBS cluster job submission in Python
# Calculates percent of reads on target
# By Jean P. Elbers
# jelber2@lsu.edu
# Last modified 15 Sep 2014
###############################################################################
Usage = """

04c-percent_reads_on_target.py

Command:
cd InDir = /work/jelber2/immunome/merged-bams
1.Bedtools intersect
        ~/bin/bedtools2-2.19.1/bin/bedtools intersect -abam Sample.bam \
        -b /work/jelber2/reference/immunome_baits_C_picta-3.0.3.bed > Sample.intersect.bam

2.Samtools flagstat to count mapped reads
        ~/bin/samtools-0.1.19/samtools flagstat \
        Sample.intersect.bam > Sample.intersect.bam.flagstat

Directory info:

File Info:
InDir = /work/jelber2/immunome/merged-bams
Input Files = Sample.bam
              Sample.bam.flagstat
OutDir = InDir
Output Files = Sample.intersect.bam
               Sample.intersect.bam.flagstat


Usage (execute following code in InDir):

find . -name '*.bam' -not -name '*-header.bam' -exec ~/scripts/immunome/04c-percent_reads_on_target.py {} \;

"""
###############################################################################
import os, sys, subprocess, re #imports os, sys, subprocess, re modules


if len(sys.argv)<2:
    print Usage
else:
    FileList = sys.argv[1:]
    RefDir = "/work/jelber2/reference"
    InDir = "/work/jelber2/immunome/merged-bams"
    for InFileName in FileList: # so samtools grabs only the file names (i.e., Samples)
        FileSuffix = ".bam" # string to remove from InFileName
        FilePrefix = "./"
        Samplepre = InFileName.replace(FileSuffix,'') # creates Samplepre string
        Sample = Samplepre.replace(FilePrefix,'') # creates Sample string
        # Customize your job options here
        Queue = "single"
        Allocation = "hpc_gopo02"
        Processors = "nodes=1:ppn=1"
        WallTime = "00:10:00"
        LogOut = InDir
        LogMerge = "oe"
        JobName = "percent_reads_on_target-%s" % (Sample)
        Command ="""
        ~/bin/bedtools2-2.19.1/bin/bedtools intersect -abam %s.bam -b /work/jelber2/reference/immunome_baits_C_picta-3.0.3.bed > %s.intersect.bam

        ~/bin/samtools-0.1.19/samtools flagstat %s.intersect.bam > %s.intersect.bam.flagstat""" % \
        (Sample, Sample,
        Sample, Sample)

        JobString = """
        #!/bin/bash
        #PBS -q %s
        #PBS -A %s
        #PBS -l %s
        #PBS -l walltime=%s
        #PBS -o %s
        #PBS -j %s
        #PBS -N %s

        cd %s
        %s\n""" % (Queue, Allocation, Processors, WallTime, LogOut, LogMerge, JobName, InDir, Command)

        #Create pipe to qsub
        proc = subprocess.Popen(['qsub'], shell=True,
          stdin=subprocess.PIPE, stdout=subprocess.PIPE, close_fds=True)
        (child_stdout, child_stdin) = (proc.stdout, proc.stdin)

        #Print JobString
        JobName = proc.communicate(JobString)[0]
        print JobString
        print JobName
