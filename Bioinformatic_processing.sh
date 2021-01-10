#####
# Supplementary code A
# Bioinformatic processing pipeline
# by Antton Alberdi
# antton.alberdi@sund.ku.dk
# January 2021
#####

###############
# A1) DATA INPUT
###############

# Each individual dataset is processed based on 5 pieces of information.
## 1) Datafile name (as stored in ERDA)
## 2) Forward primer location at 16S rRNA gene nucleotide sequence (e.g. 515 for primer 515f)
## 3) Reverse primer location at 16S rRNA gene nucleotide sequence (e.g. 806 for primer 806r)
## 4) Forward primer sequence
## 5) Reverse primer sequence
## e.g. D00368,515,806,GTGCCAGCMGCCGCGGTAA,GGACTACHVGGGTWTCTAAT
# This information is included as a shell string object called "line"
# line="D00368,515,806,GTGCCAGCMGCCGCGGTAA,GGACTACHVGGGTWTCTAAT"
#
# The working directory also needs to be specified
# workdir="/home/thisstudy"

#
# A1.1) Get data from input line
#

# Information necessary for processing the file is parsed from the data input file.

sample=$(echo $line | cut -d',' -f1)
start=$(echo $line | cut -d',' -f2)
end=$(echo $line | cut -d',' -f3)
forward=$(echo $line | cut -d',' -f4)
reverse=$(echo $line | cut -d',' -f5)

#Get complements
forward_rev=$(echo $forward | tr "[ATGCYRSWKMBDHVatgcyrswkmbdhv]" "[TACGRYSWMKVHDBtacgryswmkvhdb]" | rev)
reverse_rev=$(echo $reverse | tr "[ATGCYRSWKMBDHVatgcyrswkmbdhv]" "[TACGRYSWMKVHDBtacgryswmkvhdb]" | rev)

#Create directory
mkdir ${workdir}/${sample}

#Print to log file
echo "" > ${workdir}/${sample}/${sample}.log
echo "####################################" >> ${workdir}/${sample}/${sample}.log
echo "##### Processing sample ${sample} #####" >> ${workdir}/${sample}/${sample}.log
echo "####################################" >> ${workdir}/${sample}/${sample}.log
echo "" >> ${workdir}/${sample}/${sample}.log
echo "Start position: ${start}" >> ${workdir}/${sample}/${sample}.log
echo "End position: ${end}" >> ${workdir}/${sample}/${sample}.log
echo "Forward primer: ${forward}" >> ${workdir}/${sample}/${sample}.log
echo "Reverse primer: ${reverse}" >> ${workdir}/${sample}/${sample}.log
echo "" >> ${workdir}/${sample}/${sample}.log

###############
# A2) TRANSFER FILES FROM REPOSITORY TO SERVER
###############

# Raw data files are stored at the Electronic Research Data Archive of the University of Copenhagen.
# The data are retrieved through a sftp connection for their processing.

#Print to log
now=$(date +"%Y-%m-%d %H:%M:%S")
echo "$now | Getting data from ERDA" >> ${workdir}/${sample}/${sample}.log

#Create directory
mkdir ${workdir}/${sample}/raw
cd ${workdir}/${sample}/raw

#Generate transfer file
echo "cd microbiome_database" > transfer.sh
echo "get -r ${sample}* " >> transfer.sh
echo "exit" >> transfer.sh

#Perform the transfer
sftp -b transfer.sh [USERNAME@DOMAIN]
gunzip ${sample}*

#Test whether the transfer was accomplished succesfully
file=${sample}.fastq
filea=${sample}.fasta
file1=${sample}_1.fastq
file2=${sample}_2.fastq

if test -f "$file"; then
    echo "  File ${file} succesfully transferred from ERDA" >> ${workdir}/${sample}/${sample}.log
    format="FASTQ"
elif test -f "$filea"; then
    echo "  File ${filea} succesfully transferred from ERDA" >> ${workdir}/${sample}/${sample}.log
    format="FASTA"
elif test -f "$file1"; then
   if test -f "$file2"; then
      echo "  Files ${file1} and ${file2} succesfully transferred from ERDA" >> ${workdir}/${sample}/${sample}.log
      format="FASTQ"
   else
      echo "  There was an error when transferring file(s) from ERDA" >> ${workdir}/${sample}/${sample}.log
      continue
   fi
else
    echo "  There was an error when transferring file(s) from ERDA" >> ${workdir}/${sample}/${sample}.log
    continue
fi

#Remove transfer file
rm transfer.sh

cd ${workdir}

###############
# A3) GET DATA PROPERTIES
###############

# General data properties such as average read length, number of reads, phred score format and
# whether it is raw SR, PE or pre-processed sequence (merged paired reads) are inferred from the files.
# Based on this information, the pipeline conducts different procedures for processing the data.

#Identify sequencing mode (SR=1, PE=2)
seqmode=$(ls ${workdir}/${sample}/raw/${sample}*.fast[a,q] | wc -l)

#Identify read length, sequencing type and phred type > print to stats
#PE
if [ "$seqmode" -eq 2 ] && [ "$format" == "FASTQ" ];then
readlength=$(awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count}' ${workdir}/${sample}/raw/${sample}_1.fastq)
rawlinenum=$(cat ${workdir}/${sample}/raw/${sample}_1.fastq | wc -l)
rawseqnum=$(( $rawlinenum / 4 ))
type=PE
phred=$(head -n 100 ${workdir}/${sample}/raw/${sample}_1.fastq | awk '{if(NR%4==0) printf("%s",$0);}' |  od -A n -t u1 | awk 'BEGIN{min=100;max=0;}{for(i=1;i<=NF;i++) {if($i>max) max=$i; if($i<min) min=$i;}}END{if(max<=74 && min<59) print "Phred+33"; else if(max>73 && min>=64) print "Phred+64"; else if(min>=59 && min<64 && max>73) print "Solexa+64"; else print "Unknown score encoding";}')
#SR
elif [ "$seqmode" -eq 1 ] && [ "$format" == "FASTQ" ];then
readlength=$(awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count}' ${workdir}/${sample}/raw/${sample}.fastq)
rawlinenum=$(cat ${workdir}/${sample}/raw/${sample}.fastq | wc -l)
rawseqnum=$(( $rawlinenum / 4 ))
type=SR
phred=$(head -n 100 ${workdir}/${sample}/raw/${sample}.fastq | awk '{if(NR%4==0) printf("%s",$0);}' |  od -A n -t u1 | awk 'BEGIN{min=100;max=0;}{for(i=1;i<=NF;i++) {if($i>max) max=$i; if($i<min) min=$i;}}END{if(max<=74 && min<59) print "Phred+33"; else if(max>73 && min>=64) print "Phred+64"; else if(min>=59 && min<64 && max>73) print "Solexa+64"; else print "Unknown score encoding";}')
#SR
elif [ "$seqmode" -eq 1 ] && [ "$format" == "FASTA" ];then
readlength=$(awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count}' ${workdir}/${sample}/raw/${sample}.fasta)
rawlinenum=$(cat ${workdir}/${sample}/raw/${sample}.fasta | wc -l)
rawseqnum=$(( $rawlinenum / 2 ))
type=SR
phred=NA
else
echo "The file format is not correct"
fi

# If no reverse primer is expected, change to SRF.
# This is the data type in which only the forward reads are provided. This is a common data type, as it was
# standardised in the Earth Microbiome Project for generating 16S V4 seequences.
if [ -z "$end" ]; then
type=SRF
fi

#If no forward primer is expected, change to SRR
# This type of data is far less common, yet some studies only provide reverse reads due to errors or poor quality of forward reads.
if [ -z "$start" ]; then
type=SRR
fi

#Print to stats file
printf "File format\t$format\n" > ${workdir}/${sample}/${sample}.stats
printf "Phred format\t$phred\n" > ${workdir}/${sample}/${sample}.stats
printf "Read length\t$readlength\n" > ${workdir}/${sample}/${sample}.stats
printf "Raw\t$rawseqnum\n" >> ${workdir}/${sample}/${sample}.stats

#Print to log file
echo "  Average read length of the files is ${readlength}" >> ${workdir}/${sample}/${sample}.log
echo "  The data format is ${format}" >> ${workdir}/${sample}/${sample}.log
echo "  The phred format is ${phred}" >> ${workdir}/${sample}/${sample}.log
echo "  The number of ${type} raw sequences is ${rawseqnum}" >> ${workdir}/${sample}/${sample}.log

###############
# A4) CHECK WHETHER THE FILE CONTAINS PRIMERS
###############

# Some data files, mostly pre-processed single read sequences have primers already trimmed.
# This script detects whether the reads contain the declared primers.
# It is also useful to identify whether the primer information provided in the publication is correct.

#Print to log
now=$(date +"%Y-%m-%d %H:%M:%S")
echo "$now | Checking whether sequences contain primers" >> ${workdir}/${sample}/${sample}.log

module load tools anaconda3/4.4.0

if [ "$type" == "SR" ] && [ "$format" == "FASTQ" ];then
seqswithprimers2=$(cutadapt -e 0.15 -a ${forward}...${reverse_rev} -a ${forward_rev}...${reverse} --trimmed-only ${workdir}/${sample}/raw/${sample}.fastq | wc -l)
seqswithprimers=$(($seqswithprimers2 / 4))
primerper1=$(($seqswithprimers * 100))
   if [ $primerper1 -ne 0 ];then
   primerper=$(($primerper1 / $rawseqnum))
   else
   primerper=0
   fi
elif [ "$type" == "SR" ] && [ "$format" == "FASTA" ];then
seqswithprimers2=$(cutadapt -e 0.15 -a ${forward}...${reverse_rev} -a ${forward_rev}...${reverse} --trimmed-only ${workdir}/${sample}/raw/${sample}.fasta | wc -l)
seqswithprimers=$(($seqswithprimers2 / 2))
primerper1=$(($seqswithprimers * 100))
   if [ $primerper1 -ne 0 ];then
   primerper=$(($primerper1 / $rawseqnum))
   else
   primerper=0
   fi
elif [ "$type" == "SRF" ] && [ "$format" == "FASTQ" ];then
seqswithprimers2=$(cutadapt -e 0.15 -g ${forward} --trimmed-only ${workdir}/${sample}/raw/${sample}.fastq | wc -l)
seqswithprimers=$(($seqswithprimers2 / 4))
primerper1=$(($seqswithprimers * 100))
   if [ $primerper1 -ne 0 ];then
   primerper=$(($primerper1 / $rawseqnum))
   else
   primerper=0
   fi
elif [ "$type" == "SRF" ] && [ "$format" == "FASTA" ];then
seqswithprimers2=$(cutadapt -e 0.15 -g ${forward} --trimmed-only ${workdir}/${sample}/raw/${sample}.fasta | wc -l)
seqswithprimers=$(($seqswithprimers2 / 2))
primerper1=$(($seqswithprimers * 100))
   if [ $primerper1 -ne 0 ];then
   primerper=$(($primerper1 / $rawseqnum))
   else
   primerper=0
   fi
elif [ "$type" == "SRR" ] && [ "$format" == "FASTQ" ];then
seqswithprimers2=$(cutadapt -e 0.15 -g ${reverse} --trimmed-only ${workdir}/${sample}/raw/${sample}.fastq | wc -l)
seqswithprimers=$(($seqswithprimers2 / 4))
primerper1=$(($seqswithprimers * 100))
   if [ $primerper1 -ne 0 ];then
   primerper=$(($primerper1 / $rawseqnum))
   else
   primerper=0
   fi
elif [ "$type" == "SRR" ] && [ "$format" == "FASTA" ];then
seqswithprimers2=$(cutadapt -e 0.15 -g ${reverse} --trimmed-only ${workdir}/${sample}/raw/${sample}.fasta | wc -l)
seqswithprimers=$(($seqswithprimers2 / 2))
primerper1=$(($seqswithprimers * 100))
   if [ $primerper1 -ne 0 ];then
   primerper=$(($primerper1 / $rawseqnum))
   else
   primerper=0
   fi
elif [ "$type" == "PE" ];then
seqswithprimers2=$(cutadapt -e 0.15 -a ${forward} -a ${reverse} -a ${forward_rev} -a ${reverse_rev} --trimmed-only ${workdir}/${sample}/raw/${sample}_1.fastq | wc -l)
seqswithprimers=$(($seqswithprimers2 / 4))
primerper1=$(($seqswithprimers * 100))
   if [ $primerper1 -ne 0 ];then
   primerper=$(($primerper1 / $rawseqnum))
   else
   primerper=0
   fi
else
echo "The file format is not correct"
fi

#Are primers present?
if [ "$primerper" -gt 5 ];then
primersin=TRUE
echo "  The sequencing reads contain primer sequences" >> ${workdir}/${sample}/${sample}.log
else
primersin=FALSE
echo "  The sequencing reads do not contain primer sequences" >> ${workdir}/${sample}/${sample}.log
fi

printf "Contains primers\t$primersin\n" >> ${workdir}/${sample}/${sample}.stats

###############
# A5) QUALITY-FILTER AND COLLAPSE
###############

# Reads with average phred scores below 25 are filtered out, and PE reads merged based on their overlapping fragments

if [ "$format" == "FASTQ" ];then

#Print to log
now=$(date +"%Y-%m-%d %H:%M:%S")
echo "$now | Quality filtering and merging" >> ${workdir}/${sample}/${sample}.log

mkdir ${workdir}/${sample}/qualfilt
module load tools gcc AdapterRemoval/2.2.2

#Get marker length data
if [ "$type" == "SRF" ] || [ "$type" == "SRR" ];then
   markerlength=${readlength%.*}
else
   if [ "$primersin" == "TRUE" ];then
      #Marker length with primers
      markerlength=$(($end - $start))
   else
      #Marker length without primers
      forwardlength=$(echo $forward | wc -c)
      reverselength=$(echo $reverse | wc -c)
      primerlength=$(($forwardlength + $reverselength))
      markerlength=$((($end - $start)-$primerlength))
   fi
fi
echo "  The estimated marker length is ${markerlength}" >> ${workdir}/${sample}/${sample}.log

#Security range
minlength=$(($markerlength - 75))
maxlength=$(($markerlength + 75))

#PE (quality-filter & collapse)
if [ "$type" == "PE" ];then
   if [ "$phred" == "Phred+64" ];then
   #Phred type Phred+64
   AdapterRemoval --file1 ${workdir}/${sample}/raw/${sample}_1.fastq --file2 ${workdir}/${sample}/raw/${sample}_2.fastq --basename ${workdir}/${sample}/qualfilt/${sample} --threads 40 --collapse --minalignmentlength 5 --minquality 20 --qualitybase 64 --minlength ${minlength} --maxlength ${maxlength} --qualitymax 62
   else
   #Phred type Phred+33 (most common)
   AdapterRemoval --file1 ${workdir}/${sample}/raw/${sample}_1.fastq --file2 ${workdir}/${sample}/raw/${sample}_2.fastq --basename ${workdir}/${sample}/qualfilt/${sample} --threads 40 --collapse --minalignmentlength 5 --minquality 20 --minlength ${minlength} --maxlength ${maxlength} --qualitymax 93
   fi
cat ${workdir}/${sample}/qualfilt/${sample}.collapsed ${workdir}/${sample}/qualfilt/${sample}.collapsed.truncated > ${workdir}/${sample}/qualfilt/${sample}
rm ${workdir}/${sample}/qualfilt/${sample}.*
mv ${workdir}/${sample}/qualfilt/${sample} ${workdir}/${sample}/qualfilt/${sample}.fastq
#SR (only quality-filter)
elif [ "$type" == "SR" ] || [ "$type" == "SRF" ] || [ "$type" == "SRR" ];then
   if [ "$phred" == "Phred+64" ];then
   #Phred type Phred+64
   AdapterRemoval --file1 ${workdir}/${sample}/raw/${sample}.fastq --basename ${workdir}/${sample}/qualfilt/${sample} --threads 40 --minquality 25 --qualitybase 64 --minlength ${minlength} --maxlength ${maxlength}
   else
   #Phred type Phred+33 (most common)
   AdapterRemoval --file1 ${workdir}/${sample}/raw/${sample}.fastq --basename ${workdir}/${sample}/qualfilt/${sample} --threads 40 --minquality 25 --minlength ${minlength} --maxlength ${maxlength}
   fi
mv ${workdir}/${sample}/qualfilt/${sample}.truncated ${workdir}/${sample}/qualfilt/${sample}.fastq

#Use BBDuk to filter quality if AdapterRemoval fails (sometimes happens with IonTorrent data) (added in December 4th 2020)
if [ ! -s ${workdir}/${sample}/qualfilt/${sample}.fastq ];then
module load tools bbmap/38.35 java/1.8.0
bbduk.sh in=${workdir}/${sample}/raw/${sample}.fastq out=${workdir}/${sample}/qualfilt/${sample}.fastq maq=18
fi

else
echo "The file format is not correct"
fi

#Print to stats file
qflinenum=$(cat ${workdir}/${sample}/qualfilt/${sample}.fastq | wc -l)
qfseqnum=$(( $qflinenum /4 ))
printf "Quality filtered\t$qfseqnum\n" >> ${workdir}/${sample}/${sample}.stats
echo "  The number of quality-filtered sequences is ${qfseqnum}" >> ${workdir}/${sample}/${sample}.log

else
#Print to log
now=$(date +"%Y-%m-%d %H:%M:%S")
echo "$now | Skipping quality filtering and merging because data are already in fasta format." >> ${workdir}/${sample}/${sample}.log
printf "Quality filtered\t$rawseqnum\n" >> ${workdir}/${sample}/${sample}.stats
fi

#Identify sequence length and print to stats
if [ "$format" == "FASTQ" ];then
qfseqlength=$(awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count}' ${workdir}/${sample}/qualfilt/${sample}.fastq)
elif [ "$format" == "FASTA" ];then
qfseqlength=$readlength
else
echo "  The format is incorrect" >> ${workdir}/${sample}/${sample}.log
fi
printf "Sequence length (before trimming)\t$qfseqlength\n" >> ${workdir}/${sample}/${sample}.stats
echo "  The sequence length after quality filtering is ${qfseqlength}" >> ${workdir}/${sample}/${sample}.log

###############
# A6) REMOVE PRIMERS
###############

# Primer sequences are trimmed if these have been detected previously.
# Otherwise, files are directly copied from qualfilt (in case of fastq reads) or raw (in case of fasta sequences) directories

#Print to log
now=$(date +"%Y-%m-%d %H:%M:%S")
echo "$now | Trimming primers" >> ${workdir}/${sample}/${sample}.log

mkdir ${workdir}/${sample}/primertrimmed

if [ "$primersin" == "TRUE" ];then
module load tools anaconda3/4.4.0 libgtextutils/0.7 fastx_toolkit/0.0.14
#Trim primers
   if [ "$format" == "FASTQ" ];then
        if [ "$type" == "SRF" ];then
         cutadapt -e 0.15 -g ${forward} -o ${workdir}/${sample}/primertrimmed/${sample}.fastq ${workdir}/${sample}/qualfilt/${sample}.fastq
        elif [ "$type" == "SRR" ];then
         cutadapt -e 0.15 -g ${reverse} -o ${workdir}/${sample}/primertrimmed/${sample}.fastq ${workdir}/${sample}/qualfilt/${sample}.fastq
        else
         cutadapt -e 0.15 -a ${forward}...${reverse_rev} -a ${forward_rev}...${reverse} -o ${workdir}/${sample}/primertrimmed/${sample}.fastq ${workdir}/${sample}/qualfilt/${sample}.fastq
        fi
   else
        if [ "$type" == "SRF" ];then
         cutadapt -e 0.15 -g ${forward} -o ${workdir}/${sample}/primertrimmed/${sample}.fasta ${workdir}/${sample}/raw/${sample}.fasta
        elif [ "$type" == "SRR" ];then
         cutadapt -e 0.15 -g ${reverse} -o ${workdir}/${sample}/primertrimmed/${sample}.fasta ${workdir}/${sample}/raw/${sample}.fasta
        else
         cutadapt -e 0.15 -a ${forward}...${reverse_rev} -a ${forward_rev}...${reverse} -o ${workdir}/${sample}/primertrimmed/${sample}.fasta ${workdir}/${sample}/raw/${sample}.fasta
        fi
   fi
else
   if [ "$format" == "FASTQ" ];then
   cp ${workdir}/${sample}/qualfilt/${sample}.fastq ${workdir}/${sample}/primertrimmed/${sample}.fastq
   else
   cp ${workdir}/${sample}/raw/${sample}.fasta ${workdir}/${sample}/primertrimmed/${sample}.fasta
   fi
fi

#Print to stats file
if [ "$format" == "FASTQ" ];then
pritrilinenum=$(cat ${workdir}/${sample}/primertrimmed/${sample}.fastq | wc -l)
pritriseqnum=$(( $pritrilinenum / 4 ))
pritriseqlength=$(awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count}' ${workdir}/${sample}/primertrimmed/${sample}.fastq)
else
pritrilinenum=$(cat ${workdir}/${sample}/primertrimmed/${sample}.fasta | wc -l)
pritriseqnum=$(( $pritrilinenum / 2 ))
pritriseqlength=$(awk '{if(NR%4==2) {count++; bases += length} } END{print bases/count}' ${workdir}/${sample}/primertrimmed/${sample}.fasta)
fi
printf "Primer trimmed\t$pritriseqnum\n" >> ${workdir}/${sample}/${sample}.stats

#Identify sequence length and print to stats
printf "Sequence length (after trimming)\t$pritriseqlength\n" >> ${workdir}/${sample}/${sample}.stats
echo "  The sequence length after primer trimming is ${pritriseqlength}" >> ${workdir}/${sample}/${sample}.log

###############
# A7) CONVERT TO FASTA
###############

# Fastq file is converted into fasta.

if [ "$format" == "FASTQ" ];then
sed -n '1~4s/^@/>/p;2~4p' ${workdir}/${sample}/primertrimmed/${sample}.fastq > ${workdir}/${sample}/primertrimmed/${sample}.fasta
fi

###############
# A8) ADD SAMPLE NAME
###############

# The sample name is added to the headers of the fasta file.
# This is necessary to avoid issues in the downstream OTU table creation with usearch otutab.

sed -i "s/^>/>sample=${sample};/g" ${workdir}/${sample}/primertrimmed/${sample}.fasta

###############
# A9) CALCULATE MINIMUM CUTOFF VALUE
###############

# The minimum copy number threshold based on sequencing depth is calculated.
# The threshold is set at 0.01% of the total sequencing depth.
# If the sequencing depth is below 10000 reads, the copy number threshold is set at 2,
# which implies removing singletones
# This filtering is applied to avoid singletons and other rare reads (mostly artefactual) to exponentially increase
# the computational resources needed to process the data.

cutoff=$((pritriseqnum / 10000))
if [ "$cutoff" -le 1 ];then
cutoff=2
fi
printf "Cutoff value\t$cutoff\n" >> ${workdir}/${sample}/${sample}.stats

###############
# A10) DEREPLICATE
###############

# Duplicated reads are removed and the number of copies of each sequence is added to the header.
# The minimum copy number threshold defined above is implemented here, as the sequences with less copies
# than the specified threshold are removed.

#Print to log
now=$(date +"%Y-%m-%d %H:%M:%S")
echo "$now | Dereplicating" >> ${workdir}/${sample}/${sample}.log
echo "  Using a minimum copy number threshold of ${cutoff}" >> ${workdir}/${sample}/${sample}.log

mkdir ${workdir}/${sample}/derep
module load tools usearch/11.0.667
usearch -fastx_uniques ${workdir}/${sample}/primertrimmed/${sample}.fasta -minuniquesize ${cutoff} -sizeout -relabel ${sample}_ -fastaout ${workdir}/${sample}/derep/${sample}.fasta

#Print to stats file
uniquenum=$(grep -c ">" ${workdir}/${sample}/derep/${sample}.fasta)
printf "Dereplicated\t$uniquenum\n" >> ${workdir}/${sample}/${sample}.stats
echo "  The number of dereplicated (unique) sequences is ${uniquenum}" >> ${workdir}/${sample}/${sample}.log

###############
# A11) DENOISE
###############

# The unoise algorithm is used to create zero-ratio OTUs (zOTUs)
# DADA2 is not used because we don't have control over the experimental procedures and how samples were pooled
# for sequencing, which is needed for accurate error pattern learning.

#Print to log
now=$(date +"%Y-%m-%d %H:%M:%S")
echo "$now | Clustering OTUs" >> ${workdir}/${sample}/${sample}.log

mkdir ${workdir}/${sample}/denoised
usearch -unoise3 ${workdir}/${sample}/derep/${sample}.fasta -minsize ${cutoff} -zotus ${workdir}/${sample}/denoised/${sample}.fasta

#Convert to oneliner
perl /home/projects/ku-cbd/people/antalb/scripts/FastaToOnliner.pl ${workdir}/${sample}/denoised/${sample}.fasta > ${workdir}/${sample}/${sample}.fasta
rm -rf ${workdir}/${sample}/derep/deblur_working_dir ${workdir}/${sample}/derep/split

#Print to stats file
linenum=$(cat ${workdir}/${sample}/${sample}.fasta | wc -l)
zotunum=$(( $linenum /2 ))
printf "zOTUs\t$zotunum\n" >> ${workdir}/${sample}/${sample}.stats
echo "  The number of OTUs is ${zotunum}" >> ${workdir}/${sample}/${sample}.log

###############
# A12) CREATE ZOTU TABLE
###############

#Print to log
now=$(date +"%Y-%m-%d %H:%M:%S")
echo "$now | Generating OTU table" >> ${workdir}/${sample}/${sample}.log

usearch -otutab ${workdir}/${sample}/primertrimmed/${sample}.fasta -otus ${workdir}/${sample}/${sample}.fasta -otutabout ${workdir}/${sample}/${sample}.tab
mappednum=$(sed '1d;$d' ${workdir}/${sample}/${sample}.tab | cut -f2 | awk '{SUM+=$1}END{print SUM}')
printf "Mapped\t$mappednum\n" >> ${workdir}/${sample}/${sample}.stats
echo "  The number of reads mapped to OTU reference sequences is ${mappednum}" >> ${workdir}/${sample}/${sample}.log

###############
# A13) ASSIGN TAXONOMY
###############

#Print to log
now=$(date +"%Y-%m-%d %H:%M:%S")
echo "$now | Assigning taxonomy" >> ${workdir}/${sample}/${sample}.log

#Load required modules
module load tools hs-blastn/20160816 parallel/20190522

#Create taxonomy directory
mkdir ${workdir}/${sample}/taxonomy

#Split derep file in chuncks of 200 reads to avoid memory issues
split -a1 -l400 --additional-suffix=.fasta ${workdir}/${sample}/${sample}.fasta ${workdir}/${sample}/${sample}
samplelist=$(ls ${workdir}/${sample}/${sample}[a-z].fasta | sed 's/.*\///' | sed 's/\.fasta//')

#Run alignments
for subsample in ${samplelist}; do
    echo "Running ${subsample}"
    hs-blastn align -db /home/projects/ku-cbd/people/antalb/databases/silva/SILVA_132_SSURef_Nr99_tax_silva_modified.fasta -query ${workdir}/${sample}/${subsample}.fasta -num_threads 40 -max_target_seqs 1 -outfmt 6 > ${workdir}/${sample}/taxonomy/${subsample}.txt
done
#Merge results
cat ${workdir}/${sample}/taxonomy/${sample}[a-z].txt > ${workdir}/${sample}/taxonomy/blast.txt
rm ${workdir}/${sample}/taxonomy/${sample}[a-z].txt
rm ${workdir}/${sample}/${sample}[a-z].fasta

#Run alignment
#hs-blastn align -db /home/projects/ku-cbd/people/antalb/databases/silva/SILVA_132_SSURef_Nr99_tax_silva_modified.fasta -query ${workdir}/${sample}/${sample}.fasta -num_threads 40 -max_target_seqs 1 -outfmt 6 > ${workdir}/${sample}/taxonomy/blast.txt

if [ ! -s ${workdir}/${sample}/taxonomy/blast.txt ]; then
    echo " There was an error when assigning taxonomy (a probable memory issue)" >> ${workdir}/${sample}/${sample}.log
    continue
fi

rm ${workdir}/${sample}/${sample}.tax
#Obtain taxonomy (in parallel)
function taxononmy_parallel() {
line=${1}
workdir=${2}
sample=${3}
OTU=$(echo $line | cut -d' ' -f1)
TARGET=$(echo $line | cut -d' ' -f2)
COUNTS=$(grep -w $OTU ${workdir}/${sample}/${sample}.tab | cut -f2)
ID=$(echo $line | cut -d' ' -f3)
TAX=$(grep -w $TARGET /home/projects/ku-cbd/people/antalb/databases/silva/SILVA_132_SSURef_Nr99_tax_silva_modified.txt | sed 's/[^ ]* //' | tr ';' '\t' | awk -F '\t' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}')
printf "${OTU}\t${COUNTS}\t${ID}\t${TAX}\n" >> ${workdir}/${sample}/${sample}.tax
}
export -f taxononmy_parallel
parallel -j 40 --tmpdir ${workdir}/${sample}/taxonomy/ -k taxononmy_parallel {} ${workdir} ${sample} < ${workdir}/${sample}/taxonomy/blast.txt

#Add header
sed -i '1i OTU\tCounts\tIdentity\tDomain\tPhylum\tClass\tOrder\tFamily\tGenus' ${workdir}/${sample}/${sample}.tax

###############
# A14) AGGREGATE TAXONOMY
###############

#Print to log
now=$(date +"%Y-%m-%d %H:%M:%S")
echo "$now | Aggregating taxonomy" >> ${workdir}/${sample}/${sample}.log

module load tools gcc intel/perflibs R/4.0.0

export WORKDIR=${workdir}
export SAMPLE=${sample}
Rscript ${workdir}/bin/aggregate_taxonomy.R

#Get taxonomy numbers
domain_n=$(cat ${workdir}/${sample}/${sample}.tax_domain.txt | wc -l)
phylum_n=$(cat ${workdir}/${sample}/${sample}.tax_phylum.txt | wc -l)
class_n=$(cat ${workdir}/${sample}/${sample}.tax_class.txt | wc -l)
order_n=$(cat ${workdir}/${sample}/${sample}.tax_order.txt | wc -l)
family_n=$(cat ${workdir}/${sample}/${sample}.tax_family.txt | wc -l)
genus_n=$(cat ${workdir}/${sample}/${sample}.tax_genus.txt | wc -l)

#Add taxonomic stats to log file
echo "  Number of domains: ${domain_n}" >> ${workdir}/${sample}/${sample}.log
echo "  Number of phyla: ${phylum_n}" >> ${workdir}/${sample}/${sample}.log
echo "  Number of classes: ${class_n}" >> ${workdir}/${sample}/${sample}.log
echo "  Number of orders: ${order_n}" >> ${workdir}/${sample}/${sample}.log
echo "  Number of families: ${family_n}" >> ${workdir}/${sample}/${sample}.log
echo "  Number of genera: ${genus_n}" >> ${workdir}/${sample}/${sample}.log
