#! /usr/bin/env python

# PBS cluster job submission in Python
# BWA alignment to Painted Turtle Genome
# By Jean P. Elbers
# jelber2@lsu.edu
# Last modified 7 August 2014
###############################################################################
Usage = """

02-bwa.py - version 1.1
STEPS:
1.Uses BWA Mem to align paired reads to Painted Turtle Genome
    ~/bin/bwa-0.7.8/bwa mem \
    -M \
    -t16 \
    RefDir/C_picta-3.0.3 \
    Sample-R1-unmerged.trim.fastq.gz \
    Sample-R2-unmerged.trim.fastq.gz > ../bwa-alignment/Sample-paired.bwa.sam

2.Uses BWA Mem to align single and merged reads to Painted Turtle Genome
    ~/bin/bwa-0.7.8/bwa mem \
    -M \
    -t16 \
    RefDir/C_picta-3.0.3 \
    Sample-singlesANDmerged.trim.fastq.gz > ../bwa-alignment/Sample-singlesANDmerged.bwa.sam


Directory info:
InDir = /work/jelber2/immunome/trimmed-data
Input Files = Sample-R1-unmerged.trim.fastq.gz
              Sample-R2-unmerged.trim.fastq.gz
              Sample-singlesANDmerged.trim.fastq.gz
OutDir = /work/jelber2/immunome/bwa-alignment
Output Files = *-singlesANDmerged.bwa.sam
               *-paired.bwa.sam

Usage (execute following code in InDir):

~/scripts/immunome/02-bwa.py *.trim.fastq.gz

"""
###############################################################################
import os, sys, subprocess, re #imports os, sys, subprocess, re modules


if len(sys.argv)<2:
    print Usage
else:
    FileList = sys.argv[1:]
    RefDir = "/work/jelber2/reference"
    InDir = "/work/jelber2/immunome/trimmed-data/"
    OutDir = "/work/jelber2/immunome/bwa-alignment"
    NewFolderName = "bwa-alignment"
    os.chdir(InDir)
    os.chdir("..") # go up one directory
    if not os.path.exists(NewFolderName):
        os.mkdir(NewFolderName) # if NewFolderName does not exist, then make it
    os.chdir(InDir)
    R1 = re.compile(r"\w+\-R1\-unmerged\.trim\.fastq\.gz") # search string for obtaining only -R1-paired.fastq.gz files
    FileList = [f for f in FileList if R1.match(f)] # keeps only R1 files in FileList
    for InFileName in FileList: # so bwa only does command once for each R1 and R2 pair
        FileSuffix = "-R1-unmerged.trim.fastq.gz" # string to remove from InFileName
        Sample = InFileName.replace(FileSuffix,'') # creates Sample string
        # Customize your options here
        Queue = "workq"
        Allocation = "hpc_gopo01"
        Processors = "nodes=1:ppn=16"
        WallTime = "10:00:00"
        LogOut = OutDir
        LogMerge = "oe"
        JobName = "BWA-%s" % (Sample)
        Command ="""
        ~/bin/bwa-0.7.8/bwa mem \
        -M \
        -t16 \
        %s/C_picta-3.0.3 \
        %s-R1-unmerged.trim.fastq.gz \
        %s-R2-unmerged.trim.fastq.gz > ../bwa-alignment/%s-paired.bwa.sam

        ~/bin/bwa-0.7.8/bwa mem \
        -M \
        -t16 \
        %s/C_picta-3.0.3 \
        %s-singlesANDmerged.trim.fastq.gz > ../bwa-alignment/%s-singlesANDmerged.bwa.sam""" % \
        (RefDir, Sample, Sample, Sample,
        RefDir, Sample, Sample)

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
        jobname = proc.communicate(JobString)[0]
        print JobString
        print jobname
