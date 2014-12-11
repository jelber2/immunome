immunome
========
#Python scripts for running PBS job submissions on LSU's SuperMikeII cluster.

##The jobs in this repository are for analyzing Illumina NGS reads from a target enrichment experiment.
##Probes target painted turtle immune response genes and thus capture the "immunome" of other chelonians.


STEPS FOR QUALITY CONTROL, MAPPING, & SNP CALLING
#Download fastq.gz.zip files for the two MiSeq runs from BaseSpace
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
###Run 04b-clean_sort_addRG_markdup_realign.py
###Run 04c-percent_reads_on_target.py
###Run 04d-percent_reads_on_target_stampy.sh
###Run 05-mergeBAM_callSNPs_initial.py
###Run 06-qual_score_recal01.py
###Run 07-mergeBAM_callSNPs_recal01.py
###Run 08-qual_score_recal02.py
###Run 09-mergeBAM_callSNPs_recal02.py
###Run 10-qual_score_recal03.py
###Run 11-mergeBAM_callSNPs_recal03.py
###Run 12-seq_metrics.py

========
STEPS FOR VARIANT PREDICTION
#Download Tools First
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
    cp ../call-SNPs-recal03/ALL-samples-recal03-Q30-SNPs.vcf .
    cp ../call-SNPs-recal03/samplelist .
    ~/bin/samtools-1.1/htslib-1.1/bgzip ALL-samples-recal03-Q30-SNPs.vcf
    ~/bin/samtools-1.1/htslib-1.1/tabix -p vcf ALL-samples-recal03-Q30-SNPs.vcf.gz
    #code to split each vcf file (option '-e' to exclude rows not containing variants.)
    while read i
    do
    ~/bin/vcftools_0.1.12b/bin/vcf-subset -c $i -e ALL-samples-recal03-Q30-SNPs.vcf.gz > $i.vcf
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
###Use vcftools/bcftools on annotated snpEff files
    #code to compress then index each split-vcf, annotated by snpEff
    while read i
    do
    ~/bin/samtools-1.1/htslib-1.1/bgzip $i-annotated.vcf
    ~/bin/samtools-1.1/htslib-1.1/tabix -p vcf $i-annotated.vcf.gz
    done < samplelist
###Get only high quality snp alleles shared among:
    #Florida torts
    ~/bin/bcftools/bcftools isec -f PASS --include 'SNPEFF_EFFECT="NON_SYNONYMOUS_CODING"' -p FLsamples_shared -n=4 -w1 FL846-annotated.vcf.gz FL855-annotated.vcf.gz FL857-annotated.vcf.gz FL880-annotated.vcf.gz
    #Alabama torts
    ~/bin/bcftools/bcftools isec -f PASS --include 'SNPEFF_EFFECT="NON_SYNONYMOUS_CODING"' -p ALsamples_shared -n=4 -w1 AL102-annotated.vcf.gz AL103-annotated.vcf.gz AL106-annotated.vcf.gz AL108-annotated.vcf.gz
    #Georgia torts
    ~/bin/bcftools/bcftools isec -f PASS --include 'SNPEFF_EFFECT="NON_SYNONYMOUS_CODING"' -p GGsamples_shared -n=4 -w1 GG1044-annotated.vcf.gz GG1435-annotated.vcf.gz GG1835-annotated.vcf.gz GG462-annotated.vcf.gz
    #Louisiana torts
    ~/bin/bcftools/bcftools isec -f PASS --include 'SNPEFF_EFFECT="NON_SYNONYMOUS_CODING"' -p LAsamples_shared -n=4 -w1 LA62-annotated.vcf.gz LA66-annotated.vcf.gz LA77-annotated.vcf.gz LA78-annotated.vcf.gz
    #'ALL' torts
    #~/bin/bcftools/bcftools isec -f PASS --include 'SNPEFF_AMINO_ACID_CHANGE' -p LAsamples_shared -n=16 -w1 LA62-annotated.vcf.gz LA66-annotated.vcf.gz LA77-annotated.vcf.gz LA78-annotated.vcf.gz GG1044-annotated.vcf.gz GG1435-annotated.vcf.gz GG1835-annotated.vcf.gz GG462-annotated.vcf.gz AL102-annotated.vcf.gz AL103-annotated.vcf.gz AL106-annotated.vcf.gz AL108-annotated.vcf.gz FL846-annotated.vcf.gz FL855-annotated.vcf.gz FL857-annotated.vcf.gz FL880-annotated.vcf.gz
###Rename vcf files for each state's torts, compressed them, then index them
    mv ALsamples_shared/0000.vcf ALsamples_shared.vcf
    mv GGsamples_shared/0000.vcf GGsamples_shared.vcf
    mv FLsamples_shared/0000.vcf FLsamples_shared.vcf
    mv LAsamples_shared/0000.vcf LAsamples_shared.vcf
    ~/bin/samtools-1.1/htslib-1.1/bgzip ALsamples_shared.vcf
    ~/bin/samtools-1.1/htslib-1.1/bgzip GGsamples_shared.vcf
    ~/bin/samtools-1.1/htslib-1.1/bgzip FLsamples_shared.vcf
    ~/bin/samtools-1.1/htslib-1.1/bgzip LAsamples_shared.vcf
    ~/bin/samtools-1.1/htslib-1.1/tabix -p vcf ALsamples_shared.vcf.gz
    ~/bin/samtools-1.1/htslib-1.1/tabix -p vcf GGsamples_shared.vcf.gz
    ~/bin/samtools-1.1/htslib-1.1/tabix -p vcf FLsamples_shared.vcf.gz
    ~/bin/samtools-1.1/htslib-1.1/tabix -p vcf LAsamples_shared.vcf.gz
###Used Venny to create VennDiagram
    #http://bioinfogp.cnb.csic.es/tools/venny/index.html
    #Oliveros, J.C. (2007) VENNY. An interactive tool for comparing lists with Venn Diagrams.
    #http://bioinfogp.cnb.csic.es/tools/venny/index.html.
####1.Made three column file (chr    position    allele)
    cut -f 1-2,4 ALsamples_shared/sites.txt > ALsamples_snps.txt
    cut -f 1-2,4 GGsamples_shared/sites.txt > GGsamples_snps.txt
    cut -f 1-2,4 FLsamples_shared/sites.txt > FLsamples_snps.txt
    cut -f 1-2,4 LAsamples_shared/sites.txt > LAsamples_snps.txt
####2.Open files with gedit
    gedit LAsamples_snps.txt
    gedit ALsamples_snps.txt
    gedit GGsamples_snps.txt
    gedit FLsamples_snps.txt
####3.Copy and pasted each file into venny in order from west to east
    List1= LA torts
    List2= AL torts
    List3= GG torts
    List4= FL torts
####3.Created Venn diagram and then saved output as
    /work/jelber2/immunome/split-vcfs/immunome_venny_shared_nonsynonymous_snps.png
    /work/jelber2/immunome/split-vcfs/immunome_venny_shared_all_snps.png

========
STEPS FOR LOOKING FOR SNPs UNDER SELECTION
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



========
OLD CODE
###Download and Run PGDSpider to convert vcf file to bayescan input format
    wget http://www.cmpg.unibe.ch/software/PGDSpider/PGDSpider_2.0.7.1.zip
    unzip PGDSpider_2.0.7.1.zip 
    mv PGDSpider_2.0.7.1.zip PGDSpider_2.0.7.1

    cd ~/bin/PGDSpider_2.0.7.1/
    #use following command to generate spider.conf.xml
    java -Xmx1024m -Xms512m -jar ~/bin/PGDSpider_2.0.7.1/PGDSpider2-cli.jar \
    -inputfile /work/jelber2/immunome/bayescan/AL.vcf \
    -inputformat VCF \
    -outputfile /work/jelber2/immunome/bayescan/AL.bayescan \
    -outputformat GESTE_BAYE_SCAN
    
    #edit spider.conf.xml
    nano spider.conf.xml #to add path to samtools
    #change
    <entry key="PathBcftools"></entry>
    #to
    <entry key="PathBcftools">/home/jelber2/bin/samtools-0.1.19/bcftools/bcftools</entry>
    #change
    <entry key="PathSamtools"></entry>
    #to
    <entry key="PathSamtools">/home/jelber2/bin/samtools-0.1.19/samtools</entry>
    #save and exit
    
    #run following command to generate spid file
    java -Xmx1024m -Xms512m -jar ~/bin/PGDSpider_2.0.7.1/PGDSpider2-cli.jar \
    -inputfile /work/jelber2/immunome/bayescan/AL.vcf \
    -inputformat VCF \
    -outputfile /work/jelber2/immunome/bayescan/AL.bayescan \
    -outputformat GESTE_BAYE_SCAN
    
    #should see error:
    No SPID file specified containing preanswered conversion questions.
    A template SPID file was saved under: /work/jelber2/immunome/bayescan/template_VCF_GESTE_BAYE_SCAN.spid

    #edit the spid file
    nano /work/jelber2/immunome/bayescan/template_VCF_GESTE_BAYE_SCAN.spid
        example spid file:
        saved as:/work/jelber2/immunome/bayescan/VCF_GESTE_BAYE_SCAN.spid
            # spid-file generated: Thu Sep 25 21:32:59 CDT 2014

            # VCF Parser questions
            PARSER_FORMAT=VCF

            # Do you want to include a file with population definitions?
            VCF_PARSER_POP_QUESTION=TRUE
            # Only input following regions (refSeqName:start:end, multiple regions: whitespace sepa$
            VCF_PARSER_REGION_QUESTION=
            # What is the ploidy of the data?
            VCF_PARSER_PLOIDY_QUESTION=DIPLOID
            # Only output following individuals (ind1, ind2, ind4, ...):
            VCF_PARSER_IND_QUESTION=
            # Output genotypes as missing if the read depth of a position for the sample is below:
            VCF_PARSER_READ_QUESTION=20
            # Take most likely genotype if "PL" or "GL" is given in the genotype field?
            VCF_PARSER_PL_QUESTION=TRUE
            # Do you want to exclude loci with only missing data?
            VCF_PARSER_EXC_MISSING_LOCI_QUESTION=TRUE
            # Select population definition file:
            VCF_PARSER_POP_FILE_QUESTION=/work/jelber2/immunome/bayescan/population_definitions.txt
            # Only output SNPs with a phred-scaled quality of at least:
            VCF_PARSER_QUAL_QUESTION=30
            # Do you want to include non-polymorphic SNPs?
            VCF_PARSER_MONOMORPHIC_QUESTION=Yes
            # Output genotypes as missing if the phred-scale genotype quality is below:
            VCF_PARSER_GTQUAL_QUESTION=30

            # GESTE / BayeScan Writer questions
            WRITER_FORMAT=GESTE_BAYE_SCAN

            # Specify which data type should be included in the GESTE / BayeScan file  (GESTE / Bay$
            GESTE_BAYE_SCAN_WRITER_DATA_TYPE_QUESTION=SNP

        /work/jelber2/immunome/bayescan/population_definitions.txt
            AL102	pop_1
            AL103	pop_1
            AL106	pop_1
            AL108	pop_1
            FL846	pop_2
            FL855	pop_2
            FL857	pop_2
            FL880	pop_2
            GG1044	pop_3
            GG1435	pop_3
            GG1835	pop_3
            GG462	pop_3
            LA62	pop_4
            LA66	pop_4
            LA77	pop_4
            LA78	pop_4

    #PGDSpider_run.py
        java -Xmx1024m -Xms512m -jar ~/bin/PGDSpider_2.0.7.1/PGDSpider2-cli.jar \
        -inputfile /work/jelber2/immunome/call-SNPs-recal03/ALL-samples-recal03-Q30-SNPs.vcf \
        -inputformat VCF \
        -outputfile /work/jelber2/immunome/bayescan/ALL-samples-recal03-Q30-SNPs.bayescan.input \
        -outputformat GESTE_BAYE_SCAN \
        -spid /work/jelber2/immunome/bayescan/VCF_GESTE_BAYE_SCAN.spid


-Run BayeScan

        ~/bin/BayeScan2.1/binaries/BayeScan2.1_linux64bits \
        /work/jelber2/immunome/bayescan/ALL-samples-recal03-Q30-SNPs.bayescan.input \
        -snp \
        -od . \
        -threads 16

OLD
#Do all possible comparision
    ~/bin/bcftools/bcftools isec -n -4 -f PASS -p comparisions ALsamples_shared.vcf.gz GGsamples_shared.vcf.gz FLsamples_shared.vcf.gz LAsamples_shared.vcf.gz
#does not work    ~/bin/vcftools_0.1.12b/bin/vcf-isec -a -p comparisions ALsamples_shared.vcf.gz GGsamples_shared.vcf.gz FLsamples_shared.vcf.gz LAsamples_shared.vcf.gz


#Florida vs. Alabama
    ~/bin/bcftools/bcftools isec -f PASS -p GGsamples_shared -n=2 -w1 GG1044.vcf.gz GG1435.vcf.gz GG1835.vcf.gz GG462.vcf.gz



-Use bcftools stats then plot
    ~/bin/bcftools/bcftools stats -c some -f PASS -S FLsamples ALL-samples-recal03-Q30-SNPs.vcf.gz > out
    ~/bin/bcftools/plot-vcfstats -p FLsample_plots -s out


~/bin/bcftools/bcftools isec -f PASS -p FLsamples_shared --n=4 -w1 FL846.vcf.gz FL855.vcf.gz FL857.vcf.gz FL880.vcf.gz


    # command to split vcf file by sample (not implemented on SuperMikeII b/c process was quick, < 5 min)
    cd /work/jelber2/immunome/call-SNPs-recal03/ 
    while read i
    do
    java -Xmx4g -jar ~/bin/GATK-3.2.2/GenomeAnalysisTK.jar \
    -T SelectVariants \
    -R /work/jelber2/reference/C_picta-3.0.3.fa \
    -V ../call-SNPs-recal03/ALL-samples-recal03-Q30-SNPs.vcf \
    -sn $i \
    -o ../split-vcfs/$i.vcf
    done < samplelist
