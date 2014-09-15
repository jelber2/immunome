#! /usr/bin/env python

# PBS cluster job submission in Python
# To calculate percent of bases on target, etc using Picard
# By Jean P. Elbers
# jelber2@lsu.edu
# Last modified 15 Sep 2014
###############################################################################
Usage = """

12-seq_metrics.py - version 1.1: limited ram to 2GB per core
Command:
cd InDir = /work/jelber2/immunome/call-SNPs-recal03
1.CalculateHSMetrics:
    java -Xmx2g -jar ~/bin/picard-tools-1.118/CalculateHsMetrics.jar \
    BAIT_INTERVALS=/work/jelber2/reference/immunome_baits_C_picta-3.0.3.list \
    BAIT_SET_NAME=Immunome \
    TARGET_INTERVALS=/work/jelber2/reference/immunome_targetregion_C_picta-3.0.3.list \
    METRIC_ACCUMULATION_LEVEL=SAMPLE \
    R=/work/jelber2/reference/C_picta-3.0.3.fa \
    I=Sample.bam \
    O=Sample.hsmetrics.txt

InDir = /work/jelber2/immunome/call-SNPs-recal03
Input Files = *.bam


Usage (execute following code in InDir):

~/scripts/immunome/12-seq_metrics.py Sample.bam

"""
###############################################################################
import os, sys, subprocess #imports os, sys, subprocess modules


if len(sys.argv)<2:
    print Usage
else:
    FileList = sys.argv[1:]
    RefDir = "/work/jelber2/reference"
    InDir = os.getcwd()
    OutDir = InDir
    for InFileName in FileList: # do the following steps for each file in the inputstream
        FileSuffix = ".bam"
        Sample = InFileName.replace(FileSuffix,'') # create file prefix string
        # Customize your job options here
        Queue = "single"
        Allocation = "hpc_gopo02"
        Processors = "nodes=1:ppn=1"
        WallTime = "01:00:00"
        LogOut = OutDir
        LogMerge = "oe"
        JobName = "seq-metrics-%s" % (Sample)
        Command ="""
        java -Xmx2g -jar ~/bin/picard-tools-1.118/CalculateHsMetrics.jar \
        BAIT_INTERVALS=/work/jelber2/reference/immunome_baits_C_picta-3.0.3.list \
        BAIT_SET_NAME=Immunome \
        TARGET_INTERVALS=/work/jelber2/reference/immunome_targetregion_C_picta-3.0.3.list \
        METRIC_ACCUMULATION_LEVEL=SAMPLE \
        R=/work/jelber2/reference/C_picta-3.0.3.fa \
        I=%s.bam \
        O=%s.hsmetrics.txt""" % (Sample, Sample)

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
