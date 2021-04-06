# Epigenomics tasks

## Task 5: Distal regulatory activity

Run the docker:
```bash
$ sudo docker run -v $PWD:$PWD -w $PWD --rm -it dgarrimar/epigenomics_course
```
Create folders for adequate storing of exercise (Task1):

```bash
cd epigenomics_uvic
mkdir regulatory_elements
cd regulatory elements
```

Select from ENCODE the experiments regarding H3K27ac and H3K4me1 and download their metadata (Task2):

```bash
../bin/download.metadata.sh "https://www.encodeproject.org/metadata/?type=Experiment&replicates.library.biosample.donor.uuid=d370683e-81e7-473f-8475-7716d027849b&status=released&status=submitted&status=in+progress&biosample_ontology.term_name=stomach&biosample_ontology.term_name=sigmoid+colon&target.label=H3K4me1&target.label=H3K27ac"
```

Download bigBed files associated to metadata:

```bash
mkdir data
mkdir analyses
mkdir data/bigBed.files

#bigBed files
grep -F H3K27ac metadata.tsv |\     
grep -F "bigBed_narrowPeak" |\
grep -F "pseudoreplicated_peaks" |\
grep -F "GRCh38" |\
awk 'BEGIN{FS=OFS="\t"}{print $1, $10, $22}' |\
sort -k2,2 -k1,1r |\
sort -k2,2 -u > analyses/bigBed.H3K27ac.peaks.ids.txt

grep -F H3K4me1 metadata.tsv |\
grep -F "bigBed_narrowPeak" |\
grep -F "pseudoreplicated_peaks" |\
grep -F "GRCh38" |\
awk 'BEGIN{FS=OFS="\t"}{print $1, $10, $22}' |\
sort -k2,2 -k1,1r |\
sort -k2,2 -u > analyses/bigBed.H3K4me1.peaks.ids.txt

for histone_modif in H3K4me1 H3K27ac; do
    cut -f1 analyses/bigBed."$histone_modif".peaks.ids.txt |\
    while read filename; do
        wget -P data/bigBed.files "https://www.encodeproject.org/files/$filename/@@download/$filename.bigBed"
        done
done
```

Integrity of downloaded files is checked:

```bash
for file_type in bigBed bigWig; do

  # retrieve original MD5 hash from the metadata
  ../bin/selectRows.sh <(cut -f1 analyses/"$file_type".*.ids.txt) metadata.tsv | cut -f1,45 > data/"$file_type".files/md5sum.txt

  # compute MD5 hash on the downloaded files 
  cat data/"$file_type".files/md5sum.txt |\
  while read filename original_md5sum; do 
    md5sum data/"$file_type".files/"$filename"."$file_type" |\
    awk -v filename="$filename" -v original_md5sum="$original_md5sum" 'BEGIN{FS=" "; OFS="\t"}{print filename, original_md5sum, $1}' 
  done > tmp 
  mv tmp data/"$file_type".files/md5sum.txt

  # make sure there are no files for which original and computed MD5 hashes differ
  awk '$2!=$3' data/"$file_type".files/md5sum.txt

done
```

Convert bigBed files to Bed files:

```bash
mkdir data/bed.files

for histone in H3K4me1 H3K27ac; do
    cut -f1 analyses/bigBed."$histone".peaks.ids.txt |\
    while read filename; do
        bigBedToBed data/bigBed.files/"$filename".bigBed data/bed.files/"$filename".bed
    done
done
```

Retrieve peaks outside gene body from previous exercise for each tissue and convert them into Bed files:

```bash
mkdir annotation

for tissue in stomach sigmoid_colon; do
    cat ../ATAC-seq/analyses/peaks.analysis/peaks.outside.gene.body."$tissue".ATAC-seq.txt > annotation/peaks.outside.gene.body."$tissue".bed
done
```

Once all files are ready, intersection of peaks outside gene body for both tissues with H3K4me1 and H3K27ac is performed. First peaks outside gene body are intersected by one of the hystone modification peaks and then the obtained peaks are intersected by the second hystone modification peaks (first intersection output file has to be a bed file as it will be input for the second intersect):

```bash
mkdir analyses/peaks.analysis

#H3K4me1
cut -f-2 analyses/bigBed.H3K4me1.peaks.ids.txt |\
while read filename tissue; do 
    bedtools intersect -a annotation/peaks.outside.gene.body."$tissue".bed -b data/bed.files/"$filename".bed -u |\
    sort -u > analyses/peaks.analysis/open.regions.K4."$tissue".bed
done

#second intersection with H3K27ac
cut -f-2 analyses/bigBed.H3K27ac.peaks.ids.txt |\
while read filename tissue; do 
    bedtools intersect -a analyses/peaks.analysis/open.regions.K4."$tissue".bed -b data/bed.files/"$filename".bed -u |\
    sort -u > analyses/peaks.analysis/open.regions.candidate."$tissue".txt
done

wc -l open.regions.candidate.sigmoid_colon.txt
# 14215

wc -l open.regions.candidate.stomach.txt
# 8022
```

After this, we find 14,215 distal regulatory items for sigmoid colon and 8,022 for stomach

We proceed to filter by chromosome 1 regulatory elements and generate a file containing the name of the regulatory region and the start coordinate of the region (Task3):

```bash
for tissue in stomach sigmoid_colon; do
    awk '$1=="chr1"' analyses/peaks.analysis/open.regions.candidate."$tissue".txt |
    awk 'BEGIN{FS=OFS="\t"}{print $4,$2}' > analyses/peaks.analysis/regulatory.elements.starts."$tissue".tsv
done
```

After this, a tab-separated file called gene.starts.tsv is created. It  will store the name of the gene in the first column, and the start coordinate of the gene on the second column using as starting point a BED file generated in one of the previous sessions (Task4):

```bash
awk '$1=="chr1"' ../ATAC-seq/annotation/gencode.v24.protein.coding.gene.body.bed |
awk 'BEGIN{FS=OFS="\t"}{if ($6=="+"){start=$2} else {start=$3}; print $4, start}' > annotation/gene.starts.tsv
```

Download script `get.distance.py` and complete it so it returns the closest gene, the start of the gene and the distance of the regulatory element given a start coordinate (Task5):

```bash
cd ../bin
wget https://public-docs.crg.es/rguigo/Data/bborsari/UVIC/epigenomics_course/get.distance.py

nano get.distance.py

#!/usr/bin/env
#************
# LIBRARIES *
#************

import sys
from optparse import OptionParser


#*****************
# OPTION PARSING *
#*****************

parser = OptionParser()
parser.add_option("-i", "--input", dest="input")
parser.add_option("-s", "--start", dest="start")
options, args = parser.parse_args()

open_input = open(options.input)
enhancer_start = int(options.start)


#********
# BEGIN *
#********

x=1000000 # set maximum distance to 1 Mb
selectedGene="" # initialize the gene as empty
selectedGeneStart=0 # initialize the start coordinate of the gene as empty

for line in open_input.readlines(): # for each line in the input file
    gene, y = line.strip().split('\t') # split the line into two columns based on a tab 
    # define a variable called position that correspond to the integer of the start of the gene
    pos= int(y)
    # compute the absolute value of the difference between position and enhancer_start
    diff = abs(pos - enhancer_start)
    if diff < x:    # if this absolute value is lower than x
        x = diff        # this value will now be your current x
        selectedGene = gene        # save gene as selectedGene
        selectedGeneStart = pos        # save position as selectedGeneStart

print "\t".join([selectedGene, str(selectedGeneStart), str(x)])
```

```bash
python ../bin/get.distance.py --input gene.starts.tsv --start 980000
```
When running the previous command the results are the ones indicated in the distal regulatory activity.

Now we proceed to retrieve the closest gene and the distance to the closest gene using the python script for each egulatory element contained in the file `regulatory.elements.starts.tsv` (Task6):

Stomach:
```bash
cat analyses/peaks.analysis/regulatory.elements.starts.stomach.tsv | while read element start; do 
    python ../bin/get.distance.py --input annotation/gene.starts.tsv --start $start       
done > regulatoryElements.genes.distances.stomach.tsv
```
Sigmoidal colon:

```bash
cat analyses/peaks.analysis/regulatory.elements.starts.sigmoid_colon.tsv |
while read element start; do
    python ../bin/get.distance.py --input annotation/gene.starts.tsv --start $start 
done > regulatoryElements.genes.distances.sigmoid_colon.tsv
```

We get the files to be used to calculate the mean and median in R:

```bash
R
data <- read.table("regulatoryElements.genes.distances.stomach.tsv")
values <- rbind(c("stomach", round(mean(data$V3),2), median(data$V3)))

data <- read.table("regulatoryElements.genes.distances.sigmoid_colon.tsv")
values <- rbind(values, c("sigmoid colon", round(mean(data$V3),2), median(data$V3)))

colnames(values) <- c("Tissue","Mean", "Median")
values
#     Tissue          Mean       Median 
#[1,] "stomach"       "45227.05" "27735"
#[2,] "sigmoid colon" "73026.44" "35768"
```

Mean distance for stomach is 45,227.05 and median is 27,735. Sigmoid colon mean value is 73,026 and median 35,768






