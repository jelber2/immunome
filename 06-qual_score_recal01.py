#! /usr/bin/env python

# PBS cluster job submission in Python
# Uses GATK-3.2.2 BaseRecalibrator to recalibrate quality scores
# By Jean P. Elbers
# jelber2@lsu.edu
# Last modified 12 Aug 2014
###############################################################################
Usage = """

06-qual_score_recal01.py - version 1.0
Command:
cd InDir = /work/jelber2/immunome/realign-around-indels
1.Analyze patterns of covariation in the sequence dataset
    java -jar ~/bin/GATK-3.2.2/GenomeAnalysisTK.jar \
    -T BaseRecalibrator \
    -R RefDir/C_picta-3.0.3.fa \
    -I Sample-realigned.bam \
    -L /work/jelber2/reference/immunome_C_picta-3.0.3.list \
    -knownSites ../call-SNPs-initial/ALL-samples-initial-Q30-SNPs.vcf \
    -o ../call-SNPs-recal01/Sample-recal-data.table

2.Do a second pass to analyze covariation remaining after recalibration
    java -jar ~/bin/GATK-3.2.2/GenomeAnalysisTK.jar \
    -T BaseRecalibrator \
    -R RefDir/C_picta-3.0.3.fa \
    -I Sample-realigned.bam \
    -L /work/jelber2/reference/immunome_C_picta-3.0.3.list \
    -knownSites ../call-SNPs-initial/ALL-samples-initial-Q30-SNPs.vcf \
    -BQSR ../call-SNPs-recal01/Sample-recal-data.table \
    -o ../call-SNPs-recal01/Sample-post-recal-data.table

3.Generate before/after plots
    java -jar ~/bin/GATK-3.2.2/GenomeAnalysisTK.jar \
    -T AnalyzeCovariates \
    -R RefDir/C_picta-3.0.3.fa \
    -before ../call-SNPs-recal01/Sample-recal-data.table \
    -after ../call-SNPs-recal01/Sample-post-recal-data.table \
    -plots ../call-SNPs-recal01/Sample-recalibration_plots.pdf

4.Apply the recalibration to your sequence data
    java -jar ~/bin/GATK-3.2.2/GenomeAnalysisTK.jar \
    -T PrintReads \
    -R RefDir/C_picta-3.0.3.fa \
    -I Sample-realigned.bam \
    -L /work/jelber2/reference/immunome_C_picta-3.0.3.list \
    -BQSR ../call-SNPs-recal01/Sample-recal-data.table \
    -o ../call-SNPs-recal01/Sample-recal01.bam


File Info
InDir = /work/jelber2/immunome/realign-around-indels
Input Files =
       Sample-realigned.bam
       ../call-SNPs-initial/ALL-samples-initial-Q30-SNPs.vcf
OutDir = /work/jelber2/immunome/call-SNPs-recal01
Output Files = Sample-recal01.bam


Usage (execute following code in InDir):

~/scripts/immunome/06-qual_score_recal01.py *-realigned.bam

"""
###############################################################################
import os, sys, subprocess #imports os, sys, subprocess modules


if len(sys.argv)<2:
    print Usage
else:
    FileList = sys.argv[1:]
    RefDir = "/work/jelber2/reference"
    InDir = "/work/jelber2/immunome/realign-around-indels"
    OutDir1 = "call-SNPs-recal01"
    os.chdir(InDir)
    os.chdir("..") # go up one directory
    if not os.path.exists(OutDir1):
        os.mkdir(OutDir1) # if OutDir1 does not exist, then make it
    os.chdir(InDir)
    for InFileName in FileList: # do the following steps for each file in the inputstream
        FileSuffix = "-realigned.bam"
        Sample = InFileName.replace(FileSuffix,'') # create file prefix string
        # Customize your options here
        Queue = "single"
        Allocation = "hpc_gopo01"
        Processors = "nodes=1:ppn=1"
        WallTime = "12:00:00"
        LogOut = "/work/jelber2/immunome/call-SNPs-recal01"
        LogMerge = "oe"
        JobName = "qual_score_recal-%s" % (Sample)
        Command ="""
        java -jar ~/bin/GATK-3.2.2/GenomeAnalysisTK.jar \
        -T BaseRecalibrator \
        -R %s/C_picta-3.0.3.fa \
        -I %s-realigned.bam \
        -L /work/jelber2/reference/immunome_C_picta-3.0.3.list \
        -knownSites ../call-SNPs-initial/ALL-samples-initial-Q30-SNPs.vcf \
        -o ../call-SNPs-recal01/%s-recal-data.table

        java -jar ~/bin/GATK-3.2.2/GenomeAnalysisTK.jar \
        -T BaseRecalibrator \
        -R %s/C_picta-3.0.3.fa \
        -I %s-realigned.bam \
        -L /work/jelber2/reference/immunome_C_picta-3.0.3.list \
        -knownSites ../call-SNPs-initial/ALL-samples-initial-Q30-SNPs.vcf \
        -BQSR ../call-SNPs-recal01/%s-recal-data.table \
        -o ../call-SNPs-recal01/%s-post-recal-data.table

        java -jar ~/bin/GATK-3.2.2/GenomeAnalysisTK.jar \
        -T AnalyzeCovariates \
        -R %s/C_picta-3.0.3.fa \
        -before ../call-SNPs-recal01/%s-recal-data.table \
        -after ../call-SNPs-recal01/%s-post-recal-data.table \
        -plots ../call-SNPs-recal01/%s-recalibration_plots.pdf

        java -jar ~/bin/GATK-3.2.2/GenomeAnalysisTK.jar \
        -T PrintReads \
        -R %s/C_picta-3.0.3.fa \
        -I %s-realigned.bam \
        -L /work/jelber2/reference/immunome_C_picta-3.0.3.list \
        -BQSR ../call-SNPs-recal01/%s-recal-data.table \
        -o ../call-SNPs-recal01/%s-recal01.bam""" % \
        (RefDir, Sample, Sample,
        RefDir, Sample, Sample, Sample,
        RefDir, Sample, Sample, Sample,
        RefDir, Sample, Sample, Sample)

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
