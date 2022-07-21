#!bin/bash

#------------------------------------------------------------------
# Script to convert CEL to VCF and filter based on CONF (snparray and missingness)
# Wrote by Thais Crippa de Oliveira
# Based on https://github.com/freeseek/gtc2vcf#convert-affymetrix-cel-files-to-chp-files
# Last update: Oct, 2021
#------------------------------------------------------------------

######################################################################
# REQUIREMENTS

# percentil_m.R

## Install bioconda
# curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
# sh Miniconda3-latest-Linux-x86_64.sh
# conda config --add channels defaults
# conda config --add channels bioconda
# conda config --add channels conda-forge

### Intall bcftools
# conda install bcftools==1.12

### Install Plugin
    # mkdir -p ~/bin/bcftools/plugins
    # cd ~/bin/bcftools/plugins
    # wget http://software.broadinstitute.org/software/gtc2vcf/gtc2vcf_1.11-20210315.zip
    # unzip gtc2vcf_1.11-20210315.zip
    # cd
    
    export BCFTOOLS_PLUGINS="$HOME/bin/bcftools_m/plugins"

#REF B37: /home/bioinfo/ref/b37/human_g1k_v37_decoy.fasta
#REF hg38: /home/bioinf/ref/Broad_hg38/Homo_sapiens_assembly38.fasta
#REF hg19: /home/bioinf/ref/UCSC_hg19/ucsc.hg19.fasta
#REF h19Broad: /home/bioinf/ref/Broad_hg19/Homo_sapiens_assembly19.with_chr.fasta (IN USE)
#CHP files (Estela): /home/bioinf/res/mtle-snparrays-estela/20210608_185242_CHP
#n35_annotation:

#Affimetriz APT (Array Power Tools)
wget https://www.thermofisher.com/br/en/home/life-science/microarray-analysis/microarray-analysis-partners-programs/affymetrix-developers-network/affymetrix-power-tools.html
unzip unzip apt_2.11.4_linux_64_bit_x86_binaries.zip

# Dependencies
wget http://www.affymetrix.com/Auth/support/downloads/library_files/genomewidesnp6_libraryfile.zip
wget http://www.affymetrix.com/Auth/analysis/downloads/lf/genotyping/GenomeWideSNP_6/SNP6_supplemental_axiom_analysis_files.zip
wget http://www.affymetrix.com/Auth/analysis/downloads/na35/genotyping/GenomeWideSNP_6.na35.annot.csv.zip
wget http://www.affymetrix.com/Auth/support/downloads/library_files/genomewidesnp6_libraryfile.zip
unzip genomewidesnp6_libraryfile.zip CD_GenomeWideSNP_6_rev3/Full/GenomeWideSNP_6/LibFiles/GenomeWideSNP_6.{cdf,chrXprobes,chrYprobes,specialSNPs}
unzip -o SNP6_supplemental_axiom_analysis_files.zip GenomeWideSNP_6.{generic_prior.txt,apt-probeset-genotype.AxiomGT1.xml,AxiomGT1.sketch}
unzip -o GenomeWideSNP_6.na35.annot.csv.zip GenomeWideSNP_6.na35.annot.csv

# Original CEL Directory
CELFILES="/home/nfs/bipmed/raw/array/"
ls -1 "$PWD/"* ${CELFILES} > CElfiles.txt #to generate a full path of each of the sample files - generate inside directory that contain CEL files

# Variables
OUTDIR="home/thais/bipmed/chp"
NAME="bipmed"
LISTCEL="/home/thais/bipmed/chp/CElfiles.txt"
REF_FASTA="/home/bioinf/ref/Broad_hg19/Homo_sapiens_assembly19.with_chr.fasta"
CSV="/home/nfs/ref/GenomeWideSNP_6/GenomeWideSNP_6.na35.annot.csv"

##################################################
## STEP 1: CEL 2 CHP
##################################################

echo "##### STEP 1: CEL 2 CHP"

#### To axion SNParrays
#apt-probeset-genotype \
  --analysis-files-path . \
  --xml-file GenomeWideSNP_6.apt-probeset-genotype.AxiomGT1.xml \
  --out-dir $path_to_output_folder \
  --cel-files $cel_list_file \
  --special-snps GenomeWideSNP_6.specialSNPs \
  --chip-type GenomeWideEx_6 \
  --chip-type GenomeWideSNP_6 \
  --table-output false \
  --cc-chp-output \
  --write-models \
  --read-models-brlmmp GenomeWideSNP_6.generic_prior.txt

#### To genome-wide 6.0 snps
apt-probeset-genotype \
  --analysis-files-path . \
  -out-dir $OUTDIR \
  --cel-files $LISTCEL \
  -c GenomeWideSNP_6.cdf \
  --set-gender-method cn-probe-chrXY-ratio \
  --chrX-probes GenomeWideSNP_6.chrXprobes \
  --chrY-probes GenomeWideSNP_6.chrYprobes \
  --special-snps GenomeWideSNP_6.specialSNPs \
  --cc-chp-output \
  --write-models \
  --read-models-birdseed GenomeWideSNP_6.birdseed-v2.models \
  -a birdseed-v2



##################################################
## STEP 2: CHP 2 VCF
##################################################

echo "##### STEP 2: CHP 2 VCF"

### ------- Generate VCF multi-sample ---------------------------------------------------------
### Managing CHP files

#Inside cc-chp folder
ls -1 "$PWD/"*  > out.txt
grep -v "/home/thais/bipmed/chp/cc-chp/out.txt" out.txt > out1.txt
chp=$(tr '\n' ' ' < out1.txt)

### Run gtc2vcf, plugin affy2vcf

# Ex. bcftools +affy2vcf --csv <AnnotFile.csv> --fasta-ref <ref.fasta> <A.chp> ... <An.chp> --output <out.vcf>
# PS: Ref need to be the same used on the snparray
# PS1: Can run with one or more VCFs
## SNParray case of BIPMED: use hg19 USCS as reference

bcftools +affy2vcf --csv ${CSV} --fasta-ref ${REF_FASTA} $chp --output ${NAME}.vcf

# Remove Mitochondrial genoma - code as chrM
vcftools --vcf ${NAME}.vcf --not-chr chrM --recode --recode-INFO-all --out ${NAME}_m

# Sort VCF
bcftools sort -Ov ${NAME}_m.recode.vcf -o ${NAME}_m_sorted.vcf

##################################################
## STEP 3: LiftOver (if necessary)
##################################################

echo "##### STEP 3: LiftOver"

# http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/ # Chains

# Download chain --> For a single file, e.g. hg38ToHg19.over.chain.gz
# The file names reflect the assembly conversion data contained within in the format <db1>To<Db2>.over.chain.gz.
    rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg18.over.chain.gz .
        
# hg19/b37 to hg38
# if 100g of RAM was not sufficient, change to a greater value
# for a vcf with 400 samples and over 10gb, need at least 80G RAM

### Uncomment this part!!!!!!
#picard -Xms100g LiftoverVcf I=dnabr.hg38.700.vcf.gz O=dnabra_hg38ToHg19.vcf CHAIN=/home/bioinf/ref/Liftover_Chain_Files/hg38ToHg19.over.chain REJECT=rejected_variants.vcf R=/home/bioinf/ref/Broad_hg19/Homo_sapiens_assembly19.with_chr.fasta

#If there are any CHR like that: chr_*_alt, remove this snps
#grep "^chr" lifted_overb37toHg38.vcf | cut -f1 | sort -V | uniq -c
#grep -v "<patter>" lifted_overb37toHg38.vcf > lifted_overb37toHg38_m.vcf

##################################################
## STEP 4: Filter VCF
##################################################

echo "##### STEP 4: Filter VCF"

# VCF Filter based on CONF column - Similar to GQ -genotype quality-, is the phred-score of the analogous probability in a sequencing experiment - (probability that the GT is the real genotype) - set as 90% percentil.

#### Change path to input and name of output
R CMD BATCH percentil.Rscript

## ------- Back to bash ----------------------------------------------------------

# Make a test on VCFtools:
# Common Error: Expected at least 2 parts in FILTER definition: ll filters passed
# Change that:##FILTER=All filters passed
# For that: ##FILTER=<ID=PASS,Description="All filters passed">

sed 's/FILTER=All filters passed/FILTER=<ID=PASS,Description="All filters passed">/g' ${NAME}_m_sorted_Perc.vcf > ${NAME}_m_sorted_Perc_m.vcf


##################################################
## STEP 5: Merge new VCF with other VCFs (Optional)
##################################################

echo "##### STEP 5: Merge new VCF with other VCFs"

# Bcftools Merge
# bcftools: bgzip file.vcf
# bcftools index for all the gz VCF files: bcftools index file.vcf.gz
bcftools merge A.vcf.gz B.vcf.gz -o Total_EstelaMaira.vcf







