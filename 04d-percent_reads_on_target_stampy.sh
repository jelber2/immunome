#!/bin/sh

######
# Program Steps
# 1. Creates tab-delimited output file: percent_reads_on_target_stampy.txt
#    with 4 columns: sample_id, total_mapped_reads, mapped_reads_on_target,
#    and pct_mapped_reads_on_target
# 2. Creates samplelist
# 3. Reads one sample from samplelist at a time doing a set of commands
#    filling in the four columns of the percent_reads_on_target.txt
# 4. Repeats for all samples.
######

cd /work/jelber2/immunome/merged-bams
touch percent_reads_on_target_stampy.txt #1 creates output file percent_reads_on_target.txt
echo -e sample_id"\t"mapped_reads_on_target"\t"total_mapped_reads"\t"pct_mapped_reads_on_target >> percent_reads_on_target_stampy.txt #creates column headers

ls *.bam | grep -Po '^\w+' | sort -u > samplelist #2 makes samplelist

while read i #syntax to read one line at a time from samplelist
do #do the following commands
export imappedreadsontarget="$(awk NR==3 $i.intersect.bam.flagstat | grep -Po '^\d+')" #3a. create the bash variable imappedreadsontarget
export itotalmappedreads="$(awk NR==3 $i.bam.flagstat | grep -Po '^\d+')" #3b. creates the bash variable itotalmappedreads
export ipctmappedreadsontarget="$(bc -l <<< "($imappedreadsontarget/$itotalmappedreads)*100")" #3c. calculates ipct_reads_on_target
echo -e $i"\t"$imappedreadsontarget"\t"$itotalmappedreads"\t"$ipctmappedreadsontarget >> percent_reads_on_target_stampy.txt #3b. inserts info for current sample into percent_reads_on_target.txt
done < samplelist #syntax to read one line at a time from samplelist
