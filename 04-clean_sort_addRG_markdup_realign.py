#! /usr/bin/env python

# PBS cluster job submission in Python
# Clean, sort, add Read groups, and Mark duplicates, realign around indels
# By Jean P. Elbers
# jelber2@lsu.edu
# Last modified 11 Aug 2014
###############################################################################
Usage = """

04-clean_sort_addRG_markdup_realign.py - version 1.1 : added GPU= barcode
                                                       to add Read groups
                                                     fixed interval outdir
Command:
cd InDir = /work/jelber2/immunome/stampy-alignment
1.Uses samtools merge to combine stampy bam files
    ~/bin/samtools-0.1.19/samtools merge \
    Sample.stampy.bam Sample-paired.stampy.bam Sample-singlesANDmerged.stampy.bam

2.Uses samtools flagstat to get alignment metrics on stampy aligned bam file
    ~/bin/samtools-0.1.19/samtools flagstat \
    Sample.stampy.bam > ../stampy-alignment/Sample.stampy.bam.flagstat

3.Clean the initial stampy BAM file:
    java -Xmx2g -jar ~/bin/picard-tools-1.118/CleanSam.jar \
    I=Sample.stampy.bam \
    O=../clean-sort-addRG-markdup/Sample-CL.bam

4. Add read groups and sort:
    java -Xmx2g -jar ~/bin/picard-tools-1.118/AddOrReplaceReadGroups.jar \
    I=../clean-sort-addRG-markdup/Sample-CL.bam \
    O=../clean-sort-addRG-markdup/Sample-CL-RG.bam \
    SORT_ORDER=coordinate \
    RGPL=illumina \
    RGPU=barcode \
    RGLB=Lib1 \
    RGID=Sample \
    RGSM=Sample \
    VALIDATION_STRINGENCY=LENIENT

5.Mark PCR duplicates and optical duplicates:
    java -Xmx2g -jar ~/bin/picard-tools-1.118/MarkDuplicates.jar \
    I=../clean-sort-addRG-markdup/Sample-CL-RG.bam \
    O=../clean-sort-addRG-markdup/Sample-CL-RG-MD.bam \
    METRICS_FILE=../clean-sort-addRG-markdup/Sample-CL-RG-MD.metrics \
    MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=250 \
    CREATE_INDEX=true \
    ASSUME_SORTED=true \
    REMOVE_DUPLICATES=false

6.Find INDEL regions within individual BAM files
    java -Xmx2g -jar ~/bin/GATK-3.2.2/GenomeAnalysisTK.jar \
    -T RealignerTargetCreator \
    -R RefDir/C_picta-3.0.3.fa \
    -fixMisencodedQuals \ ##needed? - removed from code
    -I ../clean-sort-addRG-markdup/Sample-CL-RG-MD.bam \
    --minReadsAtLocus 4 \
    -o ../realign-around-indels/Sample.merged.intervals

7.Realign the BAM based on indel intervals:
    java -Xmx2g -jar ~/bin/GATK-3.2.2/GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -R RefDir/C_picta-3.0.3.fa \
    -fixMisencodedQuals \ ##needed? - removed from code
    -I ../clean-sort-addRG-markdup/Sample-CL-RG-MD.bam \
    -targetIntervals ../realign-around-indels/Sample.merged.intervals \
    -LOD 3.0 \
    -o ../realign-around-indels/Sample-realigned.bam

Directory info:

(1)/work/jelber2/immunome/stampy-alignment
(2)                      /clean-sort-addRG-markdup
(3)                      /realign-around-indels

InDir = /work/jelber2/immunome/stampy-alignment
Input Files = *-singlesANDmerged.stampy.bam
              *-paired.stampy.bam

Usage (execute following code in InDir):

~/scripts/immunome/04-clean_sort_addRG_markdup_realign *.stampy.bam

"""
###############################################################################
import os, sys, subprocess, re #imports os, sys, subprocess, re modules


if len(sys.argv)<2:
    print Usage
else:
    FileList = sys.argv[1:]
    RefDir = "/work/jelber2/reference"
    InDir = "/work/jelber2/immunome/stampy-alignment"
    OutDir1 = "clean-sort-addRG-markdup"
    OutDir2 = "realign-around-indels"
    os.chdir(InDir)
    os.chdir("..") # go up one directory
    if not os.path.exists(OutDir1):
        os.mkdir(OutDir1) # if OutDir1 does not exist, make it
    if not os.path.exists(OutDir2):
        os.mkdir(OutDir2) # if OutDir2 does not exist, make it
    os.chdir(InDir)
    Paired = re.compile(r"\w+\-paired\.stampy\.bam") # search string for obtaining paired sam files
    FileList = [f for f in FileList if Paired.match(f)] # keeps only paired files in FileList
    for InFileName in FileList: # so samtools grabs only the file names (i.e., Samples)
        FileSuffix = "-paired.stampy.bam" # string to remove from InFileName
        Sample = InFileName.replace(FileSuffix,'') # creates Sample string
        # Customize your job options here
        Queue = "single"
        Allocation = "hpc_gopo01"
        Processors = "nodes=1:ppn=1"
        WallTime = "10:00:00"
        LogOut = "/work/jelber2/immunome/clean-sort-addRG-markdup"
        LogMerge = "oe"
        JobName = "clean-sort-addRG-markdup-realign-%s" % (Sample)
        Command ="""
        ~/bin/samtools-0.1.19/samtools merge \
        %s.stampy.bam %s-paired.stampy.bam %s-singlesANDmerged.stampy.bam

        ~/bin/samtools-0.1.19/samtools flagstat \
        %s.stampy.bam > %s.stampy.bam.flagstat

        java -Xmx2g -jar ~/bin/picard-tools-1.118/CleanSam.jar \
        I=%s.stampy.bam \
        O=../clean-sort-addRG-markdup/%s-CL.bam

        java -Xmx2g -jar ~/bin/picard-tools-1.118/AddOrReplaceReadGroups.jar \
        I=../clean-sort-addRG-markdup/%s-CL.bam \
        O=../clean-sort-addRG-markdup/%s-CL-RG.bam \
        SORT_ORDER=coordinate \
        RGPL=illumina \
        RGPU=barcode \
        RGLB=Lib1 \
        RGID=%s \
        RGSM=%s \
        VALIDATION_STRINGENCY=LENIENT

        java -Xmx2g -jar ~/bin/picard-tools-1.118/MarkDuplicates.jar \
        I=../clean-sort-addRG-markdup/%s-CL-RG.bam \
        O=../clean-sort-addRG-markdup/%s-CL-RG-MD.bam \
        METRICS_FILE=../clean-sort-addRG-markdup/%s-CL-RG-MD.metrics \
        MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=250 \
        CREATE_INDEX=true \
        ASSUME_SORTED=true \
        REMOVE_DUPLICATES=false

        java -Xmx2g -jar ~/bin/GATK-3.2.2/GenomeAnalysisTK.jar \
        -T RealignerTargetCreator \
        -R %s/C_picta-3.0.3.fa \
        -I ../clean-sort-addRG-markdup/%s-CL-RG-MD.bam \
        --minReadsAtLocus 4 \
        -o ../realign-around-indels/%s.merged.intervals

        java -Xmx2g -jar ~/bin/GATK-3.2.2/GenomeAnalysisTK.jar \
        -T IndelRealigner \
        -R %s/C_picta-3.0.3.fa \
        -I ../clean-sort-addRG-markdup/%s-CL-RG-MD.bam \
        -targetIntervals ../realign-around-indels/%s.merged.intervals \
        -LOD 3.0 \
        -o ../realign-around-indels/%s-realigned.bam""" % \
        (Sample, Sample, Sample,
        Sample, Sample,
        Sample, Sample,
        Sample, Sample, Sample, Sample,
        Sample, Sample, Sample,
        RefDir, Sample, Sample,
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

