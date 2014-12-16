#immunome
========
###Python scripts for running PBS job submissions on LSU's SuperMikeII cluster.
###The jobs in this repository are for analyzing Illumina NGS reads from a target enrichment experiment.
###Probes target painted turtle immune response genes and thus capture the "immunome" of other chelonians.

========
#STEPS FOR QUALITY CONTROL, MAPPING, & SNP CALLING
##Download fastq.gz.zip files for the two MiSeq runs from BaseSpace
###Copy data to supermikeII
    #immunome = miseq data from 9Sep2014
    rsync --archive --stats --progress /work/jelber2/immunome/analysis_13944931_fastq.zip jelber2@mike.hpc.lsu.edu:/work/jelber2/immunome/
    #immunome2 = miseq data from 15Sep2014
    rsync --archive --stats --progress /work/jelber2/immunome2/analysis_14120117_fastq.zip jelber2@mike.hpc.lsu.edu:/work/jelber2/immunome2/
###Unzip data (on SuperMikeII)
    unzip /work/jelber2/immunome/analysis_13944931_fastq.zip
    cd Data/Intensities/BaseCalls
    mkdir /work/jelber2/immunome/fastq
    mv * /work/jelber2/immunome/fastq
    cd /work/jelber2
    cd immunome2
    unzip /work/jelber2/immunome2/analysis_14120117_fastq.zip
    cd Data/Intensities/BaseCalls
    mkdir /work/jelber2/immunome2/fastq
    mv * /work/jelber2/immunome2/fastq
###Rename fastq files using mv command and regular expressions:
    cd /work/jelber2/immunome/fastq
#Quality Control
###Run 01-trimmomatic.py on fastq.gz in:
    /work/jelber2/immunome/fastq
    /work/jelber2/immumome2/fastq (change directory in script to immunome2)
#Mapping
###Run 02-bwa.py on trim.fastq.gz on trim.fastq.gz in:
    /work/jelber2/immunome/trimmed-data
    /work/jelber2/immunome2/trimmed-data (change directory in script to immunome2)
###Run 03-stampy.py on bwa.sam in:
    /work/jelber2/immunome/bwa-alignment
    /work/jelber2/immunome2/bwa-alignment (change directory in script to immunome2)
#SNP Calling
###Run 04a-clean_sort_addRG.py stampy.bam in:
    /work/jelber2/immunome/stampy-alignment (make RGID=%s_9Sep2014, and directory to immunome)
    /work/jelber2/immunome2/stampy-alignment (make RGID=%s_15Sep2014, and directory to immunome2)
####Run 04b-clean_sort_addRG_markdup_realign.py
####Run 04c-percent_reads_on_target.py
####Run 04d-percent_reads_on_target_stampy.sh
####Run 05-mergeBAM_callSNPs_initial.py
####Run 06-qual_score_recal01.py
####Run 07-mergeBAM_callSNPs_recal01.py
####Run 08-qual_score_recal02.py
####Run 09-mergeBAM_callSNPs_recal02.py
####Run 10-qual_score_recal03.py
####Run 11-mergeBAM_callSNPs_recal03.py
####Run 12-seq_metrics.py
####Need to use beagle to improve SNPs called by Unified Genotyper
#####Downloaded beagle
    cd ~/bin
    wget http://faculty.washington.edu/browning/beagle/beagle.r1398.jar
#####Ran beagle
    cd /work/jelber2/immunome
    mkdir beagle
    java -Xmx4000m -jar ~/bin/beagle.r1398.jar \
    gtgl=/work/jelber2/immunome/call-SNPs-recal03/ALL-samples-recal03-Q30-SNPs.vcf \
    nthreads=2 \
    out=/work/jelber2/immunome/beagle/ALL-samples-recal03-Q30-SNPs-beagle
========
#STEPS FOR VARIANT PREDICTION
##Download Tools First
###Downloaded snpEff
    #ideally want to know if variants will affect protein structure and possibly immune gene function
    cd /work/jelber2/reference
    wget http://iweb.dl.sourceforge.net/project/snpeff/snpEff_latest_core.zip
    unzip snpEff_latest_core.zip
###Added Chrysemys_picta_bellii-3.0.3 to snpEff.config using nano
    cd /work/jelber2/reference/snpEff
    nano snpEff.config # added the following four lines after the Capsella_rubella_v1.0 entry (remove 4 spaces on left if cut and pasting)
    # Chrysemys_picta_bellii-3.0.3
    Chrysemys_picta_bellii-3.0.3.genome : western painted turtle
    	Chrysemys_picta_bellii-3.0.3.reference : ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3/
    	Chrysemys_picta_bellii-3.0.3.M.codonTable : Standard
###Created data directory for Chrysemys_picta_bellii-3.0.3 genome
    cd /work/jelber2/reference/snpEff
    mkdir data
    cd data
    mkdir Chrysemys_picta_bellii-3.0.3
    cd Chrysemys_picta_bellii-3.0.3
    # downloaded FASTA file
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna.gz
    # snpEff requires genome.fa file to be called "sequences.fa"
    mv GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna.gz sequences.fa.gz
    # have to unzip sequences.fa.gz
    gunzip sequences.fa.gz
    # downloaded gff3 file (i.e., gene annotation file)
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.gff.gz
    # snpEff requires gene annotation file be called "genes.gff"
    mv GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.gff.gz genes.gff.gz
    # unzipped genes.gff.gz
    gunzip genes.gff.gz
    # download protein sequences
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_protein.faa.gz
    mv GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_protein.faa.gz protein.fa.gz
    gunzip protein.fa.gz
###Built snpEff database for Chrysemys_picta_bellii-3.0.3
    cd /work/jelber2/reference/snpEff
    # used snpEff_build.py script to implement command below, which took < 30 minutes
    java -jar -Xmx8g /work/jelber2/reference/snpEff/snpEff.jar build -gff3 -v Chrysemys_picta_bellii-3.0.3 2>&1 | tee Chrysemys_picta_bellii-3.0.3.build
###Downloaded vcftools
    cd /home/jelber2/bin/
    wget http://downloads.sourceforge.net/project/vcftools/vcftools_0.1.12b.tar.gz?r=http%3A%2F%2Fsourceforge.net%2Fprojects%2Fvcftools%2Ffiles%2F&ts=1411515317&use_mirror=superb-dca2
    tar -xzf vcftools_0.1.12b.tar.gz 
    mv vcftools_0.1.12b.tar.gz vcftools_0.1.12b
    cd vcftools_0.1.12b/
    nano ~/.soft #add the following two lines to using nano.soft file
    PATH+=/home/jelber2/bin/tabix-0.2.6
    PERL5LIB = /home/jelber2/bin/vcftools_0.1.12b/perl
    resoft #to refresh soft file
    cd /home/jelber2/bin/vcftools_0.1.12b/
    make #compile vcftools
    # Path to vcftools executable
    /home/jelber2/bin/vcftools_0.1.12b/bin/vcftools
###Downloaded bcftools
    cd /home/jelber2/bin
    git clone --branch=develop git://github.com/samtools/htslib.git
    git clone --branch=develop git://github.com/samtools/bcftools.git
    cd bcftools; make
#Need to look for protein altering variants shared by samples in the same population
###Split vcf file from GATK for snpEff
    #snpEff needs ALL-samples*.vcf file split by sample (i.e., into Sample1.vcf, Sample2.vcf)
    cd /work/jelber2/immunome/call-SNPs-recal03/
    ls *-recal03.bam | grep -Po '^\w+'| sort -u | grep -v 'ALL' > samplelist
    mkdir ../split-vcfs
    cd ../split-vcfs
    zcat ../beagle/ALL-samples-recal03-Q30-SNPs-beagle.vcf.gz > ALL-samples-recal03-Q30-SNPs-beagle.vcf
    cp ../call-SNPs-recal03/samplelist .
    ~/bin/samtools-1.1/htslib-1.1/bgzip ALL-samples-recal03-Q30-SNPs-beagle.vcf
    ~/bin/samtools-1.1/htslib-1.1/tabix -p vcf ALL-samples-recal03-Q30-SNPs-beagle.vcf.gz
    #code to split each vcf file
    while read i
    do
    ~/bin/bcftools/bcftools view -s $i ALL-samples-recal03-Q30-SNPs-beagle.vcf.gz -O v -o $i.vcf
    done < samplelist
###Ran snpEff on each split vcf file
    cd /work/jelber2/immunome/split-vcfs/
    # command below to run snpEff on all samples in samplelist
    # not implemented on SuperMikeII b/c process was < 15 min
    while read i
    do
    java -Xmx4g -jar /work/jelber2/reference/snpEff/snpEff.jar \
    -v -i vcf -o gatk \
    Chrysemys_picta_bellii-3.0.3 \
    /work/jelber2/immunome/split-vcfs/$i.vcf \
    > /work/jelber2/immunome/split-vcfs/$i-snpeff.vcf
    mv snpEff_genes.txt $i-snpeff-genes.txt
    mv snpEff_summary.html $i-snpeff-summary.html
    done < samplelist
###Ran VariantAnnotator on each file
    while read i
    do
    java -Xmx4g -jar ~/bin/GATK-3.2.2/GenomeAnalysisTK.jar \
    -T VariantAnnotator \
    -R /work/jelber2/reference/C_picta-3.0.3.fa \
    -A SnpEff \
    --variant $i.vcf \
    --snpEffFile $i-snpeff.vcf \
    -L $i.vcf \
    -o $i-annotated.vcf
    done < samplelist
###Merge split, annotated vcfs
    #compress then index split files
    while read i
    do
    ~/bin/samtools-1.1/htslib-1.1/bgzip -f $i-annotated.vcf
    ~/bin/samtools-1.1/htslib-1.1/tabix -p vcf $i-annotated.vcf.gz
    done < samplelist
    #finally merge them and index
    ~/bin/bcftools/bcftools merge -f PASS -o ALL-samples-annotated.vcf.gz -O z -m both ../split-vcfs/*-annotated.vcf.gz
    ~/bin/samtools-1.1/htslib-1.1/tabix -p vcf ALL-samples-annotated.vcf.gz
###Get only high quality non-synonymous snp alleles
    cd /work/jelber2/immunome
    mkdir /venny-data
    cd venny-data
    cp ../split-vcfs/samplelist .
    #script to automate unzipping bgzipped files
    #and get only  with non-synonymous snps
    #NOTE COMMENT CODE FOR GENOTYPE AND UNCOMMENT
    #CODE FOR HAPLOTYPE IF YOU WANT HAPLOTYPES
    #PRECEEDING STEPS ASSUME YOU WANT GENOTYPES
    while read i
    do
    zcat ../split-vcfs/$i-annotated.vcf.gz | \
    grep -v '#' | grep 'PASS' | grep 'NON_SYNONYMOUS_CODING' | \
    cut -f 1-2,10 | \
    #code for same haplotypes (i.e., 0/1 =! 1/0)
    #perl -pe "s/(\w+\.\d)\t(\d+)\t(\d\|\d).+\n/\1\t\2\t\3\n/" > $i-haplotype.txt
    #code for same genotypes (i.e., 0/1 = 1/0)
    perl -pe "s/(\w+\.\d)\t(\d+)\t(\d)\|(\d).+\n/\1\t\2\t\3\t\4\n/" | awk -v OFS='\t' '{a=$3+$4;print $1,$2,a;}' - > $i-genotype.txt
    done < samplelist
####For getting snps alleles shared amongst populations
####1.Opened each set of text files (4 files per set) separately with gedit
    #e.g.,
    gedit AL*-genotype.txt
    #opens AL102-genotype.txt AL106-genotype.txt AL103-genotype.txt AL108-genotype.txt
####2.Copy and pasted the contents of each file into Venny
    #Oliveros, J.C. (2007) VENNY. An interactive tool for comparing lists with Venn Diagrams.
    #http://bioinfogp.cnb.csic.es/tools/venny/index.html
####3.Saved the intersection of all 4 samples as pop-shared-genotypes.txt (e.g., AL-shared-genotypes.txt)
####4.Saved venn diagram image as immunome_venny_shared_pop_nonsynonymous_snps_genotypes.png
    #where 'pop' is either LA, AL, GG, or FL
####5.Saved unique parts (i.e., complement) of each Sample as Sample-unique-genotypes.txt
####6.Repeated steps 1-4 using the four populations
    gedit LA-shared-genotypes.txt AL-shared-genotypes.txt GG-shared-genotypes.txt FL-shared-genotypes.txt
####7.Saved venn diagram as immunome_venny_shared_all_nonsynonymous_snps.png
####8.Saved snp alleles shared by individuals in each pop as:
    LA-only-shared-genotypes.txt
    AL-only-shared-genotypes.txt
    GG-only-shared-genotypes.txt
    FL-only-shared-genotypes.txt
####9.Get only snps unique to a population
    cat LA*-unique-genotypes.txt > LA-unique-genotypes.txt
    cat AL*-unique-genotypes.txt > AL-unique-genotypes.txt
    cat GG*-unique-genotypes.txt > GG-unique-genotypes.txt
    cat FL*-unique-genotypes.txt > FL-unique-genotypes.txt
####10.Copy and pasted the contents of each file into Venny
    gedit LA-unique-genotypes.txt AL-unique-genotypes.txt GG-unique-genotypes.txt FL-unique-genotypes.txt
####11.Saved snp alleles unique to each pop as:
    LA-only-unique-genotypes.txt
    AL-only-unique-genotypes.txt
    GG-only-unique-genotypes.txt
    FL-only-unique-genotypes.txt
###Used IGV to visualize shared snps
####Get igv
    cd ~/bin
    wget http://www.broadinstitute.org/igv/projects/downloads/IGV_2.3.40.zip
    unzip IGV_2.3.40.zip 
    mv IGV_2.3.40.zip IGV_2.3.40
    cd IGV_2.3.40/
####Make igv .genome file from C_picta genome
    cd /work/jelber2/reference
    #get FASTA(FNA) file for C_picta-3.0.3
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna.gz
    #get gff gene annotation for C_picta-3.0.3
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3/GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.gff.gz
    #unzip the two files
    gunzip GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna.gz
    gunzip GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.gff.gz
    #modify igv.sh and change 2000m to 4000m, so you have at least 4GB ram
    nano ~/bin/IGV_2.3.40/igv.sh
    #run igv
    ~/bin/IGV_2.3.40/igv.sh
    #Click on the 'Genomes' tab on the main menu then select, 'Create .genome File...'
    #Under 'Unique identifier' type 'C_picta-3.0.3'
    #Under 'Descriptive name' type 'C_picta-3.0.3'
    #Under 'FASTA file' select 'GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.fna'
    #Under 'Gene file' select 'GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.gff'
    #Click 'Ok'
    #Save as GCF_000241765.3_Chrysemys_picta_bellii-3.0.3_genomic.genome in the 
    # /work/jelber2/reference/ directory
####Make directory for images
    cd /work/jelber2/immunome
    mkdir igv-snps
    cd igv-snps
####Make target files to filter vcfs
    cut -f 1-2 ../venny-data/LA-only-unique-genotypes.txt | grep -v 'Elements' | perl -pe "s/(\w+\.\d)\t(\d+)\n/goto \1:\2\n/" > LA-only_nonysn_unique_snps.txt
    cut -f 1-2 ../venny-data/LA-only-shared-genotypes.txt | grep -v 'Elements' | perl -pe "s/(\w+\.\d)\t(\d+)\n/goto \1:\2\n/" > LA-only_nonysn_snps.txt
    cut -f 1-2 ../venny-data/AL-only-shared-genotypes.txt | grep -v 'Elements' | perl -pe "s/(\w+\.\d)\t(\d+)\n/goto \1:\2\n/" > AL-only_nonysn_snps.txt
    cut -f 1-2 ../venny-data/GG-only-shared-genotypes.txt | grep -v 'Elements' | perl -pe "s/(\w+\.\d)\t(\d+)\n/goto \1:\2\n/" > GG-only_nonysn_snps.txt
    cut -f 1-2 ../venny-data/FL-only-shared-genotypes.txt | grep -v 'Elements' | perl -pe "s/(\w+\.\d)\t(\d+)\n/goto \1:\2\n/" > FL-only_nonysn_snps.txt
    #split files into 50 sites
    split -l 100 -d LA-only_nonysn_unique_snps.txt LA-only_nonysn_unique_snps.txt-
    split -l 50 -d LA-only_nonysn_snps.txt LA-only_nonysn_snps.txt-
    split -l 50 -d AL-only_nonysn_snps.txt AL-only_nonysn_snps.txt-
    split -l 50 -d GG-only_nonysn_snps.txt GG-only_nonysn_snps.txt-
    split -l 50 -d FL-only_nonysn_snps.txt FL-only_nonysn_snps.txt-
####Run IGV-
    #loading the ALLsamples-annotated.vcf file
    ~/bin/IGV_2.3.40/igv.sh /work/jelber2/immunome/split-vcfs/ALL-samples-annotated.vcf.gz
    #Once loaded, Click on 'Tools' tab in the main menu
    #Select the batch file LA-only_nonysn_snps_target.txt-00 in the 
    #/work/jelber2/immunome/igv-snps/ directory
    #let IGV do it's thing (it will take a while)
    #repeat for each number (i.e., 01,02,03) and population (AL,GG,FL)

========
#STEPS FOR LOOKING FOR SNPs UNDER SELECTION
##Download tools
###Download BayeScan
    cd /home/jelber2/bin
    wget http://cmpg.unibe.ch/software/BayeScan/files/BayeScan2.1.zip
    unzip BayeScan2.1.zip
    mv BayeScan2.1.zip BayeScan2.1
    cd BayeScan2.1/
    cd binaries/
    chmod u+x BayeScan2.1_linux64bits # makes the file executable
    # Path to BayeScan
    /home/jelber2/bin/BayeScan2.1/binaries/BayeScan2.1_linux64bits
###Download Simple Fool's Guide to RNA-seq scripts to convert vcf file to BayeScan input format
    #runs make_bayescan_input.py using:
    #30 = min genotype quality
    #4 = min number of good quality genotypes required from each population in order for a given SNP to be included in the analysis
    #1 = min number of copies of the minor allele that are necc. for a locus to be considered trustworthy enough to be used in BayeScan
    #1 = make outfile file (used_snp_genos.txt) showing what snp genotypes were used
    #output = bayes_input.tx, snpkey.txt, low_freq_snps.txt, used_snp_genos.txt
    python /Users/jelbers/Documents/Documents/LSU/Dissertation/Immunome/fromSFG/Scripts_for_SFG/make_bayescan_input.py /Users/jelbers/Desktop/sandbox/bayescan/ALL-samples-recal03-Q30-SNPs.vcf /Users/jelbers/Desktop/sandbox/bayescan_with_sfg_scripts/populations.txt 30 4 1 1
    #copy files to SuperMikeII
    rsync --stats --progress --archive /Users/jelbers/Desktop/sandbox/bayescan_with_sfg_scripts/ jelber2@mike.hpc.lsu.edu:/work/jelber2/immunome/bayescan_with_sfg_scripts/ -n
###Run BayeScan with all loci
    ~/bin/BayeScan2.1/binaries/BayeScan2.1_linux64bits \
    /work/jelber2/immunome/bayescan_with_sfg_scripts/bayescan_input.txt \
    -snp \
    -od . \
    -o bayescan_all_loci \
    -threads 16
###Run BayeScan excluding loci with low frequency minor allele snps
    ~/bin/BayeScan2.1/binaries/BayeScan2.1_linux64bits \
    /work/jelber2/immunome/bayescan_with_sfg_scripts/bayescan_input.txt \
    -d low_freq_snps.txt
    -snp \
    -od . \
    -o bayescan_no_loci_with_low_freq_minor_alleles \
    -threads 16
