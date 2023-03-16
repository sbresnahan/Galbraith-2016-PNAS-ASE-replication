# Testing the Kinship Theory of Intragenomic Conflict in Honey Bees
Replicating pipeline and results described in Galbraith et al. (2016). *PNAS*.   
https://www.pnas.org/content/113/4/1020

## Overview of pipeline
- [Requirements](#requirements)
- [Process F0 DNA-seq libraries](#process-f0-dna-seq-libraries)
  * Sample metadata
  * Define directory variables
  * Retrieve F0 DNA-seq libraries
  * Trim adapters from F0 DNA-seq libraries
  * Generate BWA alignment index for the Amel_HAv3.1 reference genome
  * Align F0 DNA-seq libraries to Amel_HAv3.1
  * Call F0 variants on Amel_HAv3.1
  * Generate F0 reference genomes
  * Create STAR indices for each F0 genome
- [Process F1 RNA-seq libraries](#process-f1-rna-seq-libraries)
  * Sample metadata
  * Retrieve F1 RNA-seq libraries
  * Trim adapters from F1 RNA-seq libraries
  * Align F1 RNA-seq libraries to F0 genomes
- [Count F1 reads over F0 variants](#count-f1-reads-over-f0-variants)
  * Define directory variables
  * Filter variants for analysis
  * Generate BED file of each SNP:gene
  * Compute strand-wise read coverage at each SNP:gene


# Requirements

##### miniconda
https://docs.conda.io/en/latest/miniconda.html

##### sra-tools
```
conda create --name sra-tools
conda install -c bioconda -n sra-tools sra-tools
```

##### fastp
```
conda create --name fastp
conda install -c bioconda -n fastp fastp
```

##### bwa
```
conda create --name bwa
conda install -c bioconda -n bwa bwa
```

##### tophat2
```
conda create --name tophat2
conda install -c bioconda -n tophat2 tophat
```

##### freebayes
```
conda create --name freebayers
conda install -c bioconda -n freebayes freebayes
```

##### samtools
```
conda create --name samtools
conda install -c bioconda -n samtools samtools
```

##### bcftools
```
conda create --name bcftools
conda install -c bioconda -n bcftools bcftools
```

##### gatk
```
conda create --name gatk
conda install -c bioconda -n gatk gatk
```

##### bedtools
```
conda create --name bedtools
conda install -c bioconda -n bedtools bedtools
```

# Process F0 DNA-seq libraries

### Sample metadata

|    SRA     | Cross | Parent | Lineage |
| ---------- | ----- | ------ | ------- |
| SRR3037350 |  875  |    Q   |   EHB   |
| SRR3037351 |  875  |    D   |   AHB   |
| SRR3037352 |  888  |    Q   |   AHB   |
| SRR3037353 |  888  |    D   |   EHB   |
| SRR3037354 |  882  |    Q   |   EHB   |
| SRR3037355 |  882  |    D   |   AHB   |
| SRR3037356 |  894  |    Q   |   AHB   |
| SRR3037357 |  894  |    D   |   EHB   |

## Define directory variables

```
DIR_SRA="/storage/home/stb5321/scratch/galbraith/raw"
DIR_TRIM="/storage/home/stb5321/scratch/galbraith/trimmed"
DIR_INDEX="/storage/home/stb5321/scratch/galbraith/index"
DIR_ALIGN="/storage/home/stb5321/scratch/galbraith/aligned"
DIR_VARIANTS="/storage/home/stb5321/scratch/galbraith/variants"
DIR_ARG="/storage/home/stb5321/scratch/galbraith/parent_genomes"
DIR_SORT="/storage/home/stb5321/scratch/galbraith/tophat2_nomm_sort"
DIR_COUNTS="/storage/home/stb5321/scratch/galbraith/tophat2_counts"
```

## Retrieve F0 DNA-seq libraries

```
conda activate sra-tools

cd ${DIR_SRA}

SRA=("SRR3037350" "SRR3037351" "SRR3037352" "SRR3037353" \
     "SRR3037354" "SRR3037355" "SRR3037356" "SRR3037357")

for i in "${SRA[@]}"
do
  prefetch -O ${DIR_SRA} ${i}
  fasterq-dump -O ${DIR_SRA} ${DIR_SRA}/${i}.sra
  rm {i}.sra
done

conda deactivate
```

## Trim adapters from F0 DNA-seq libraries

```
conda activate fastp

for i in "${SRA[@]}"
do
  fastp -w 8 -i ${i}_1.fastq -I ${i}_2.fastq \
  -o ${DIR_TRIM}/${i}_1.fastq -O ${DIR_TRIM}/${i}_2.fastq
done

conda deactivate
```

## Generate BWA alignment index for the Amel_HAv3.1 reference genome

```
cd ${DIR_INDEX}

wget -O Amel_HAv3.1.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/254/395/GCF_003254395.2_Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz
gunzip Amel_HAv3.1.fna.gz
mv Amel_HAv3.1.fna Amel_HAv3.1.fasta

wget -O Amel_HAv3.1.gff.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/254/395/GCF_003254395.2_Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_genomic.gff.gz
gunzip Amel_HAv3.1.gff.gz

conda activate bwa

bwa index Amel_HAv3.1.fasta

conda deactivate
```

## Align F0 DNA-seq libraries to Amel_HAv3.1

1) Align [`bwa mem`].
2) Convert SAM to sorted BAM [`samtools`].

```
cd ${DIR_TRIM}

conda activate bwa

for i in "${SRA[@]}"
do
  bwa mem -t 8 ${DIR_INDEX}/Amel_HAv3.1.fasta ${i}_1.fastq ${i}_2.fastq > ${DIR_ALIGN}/${i}.sam
done

conda deactivate

conda activate samtools

cd ${DIR_ALIGN}

for i in "${SRA[@]}"
do
  samtools view -@ 8 -O BAM ${i}.sam | samtools sort -@ 8 -O BAM > ${i}.bam
done

conda deactivate
```

## Call F0 variants on Amel_HAv3.1

Using freebayes to account for differences in ploidy between `DIPLOID` and `HAPLOID` libraries.

```
DIPLOID=("SRR3037350" "SRR3037352" "SRR3037354" "SRR3037356")
HAPLOID=("SRR3037351" "SRR3037353" "SRR3037355" "SRR3037357")

cd ${DIR_INDEX}

conda activate samtools

samtools faidx Amel_HAv3.1.fasta

conda deactivate

cd ${DIR_ALIGN}

conda activate freebayes

for i in "${DIPLOID[@]}"
do
  freebayes -f ${DIR_INDEX}/Amel_HAv3.1.fasta ${i}.bam > ${DIR_VARIANTS}/${i}.vcf
done

for i in "${HAPLOID[@]}"
do
  freebayes -f ${DIR_INDEX}/Amel_HAv3.1.fasta ${i}.bam -p 1 > ${DIR_VARIANTS}/${i}.vcf
done

conda deactivate
```

## Generate F0 reference genomes

1) Remove heterozygous variants and indels [`bcftools filter`] from VCFs. NOTE: variants, MNPs, and complex variants are retained.
2) Compress resultant VCFs and index [`bgzip`, `tabix`].  
3) Integrate homozygous variants into Amel_HAv3.1 for each F0 library, separately [`gatk FastaAlternateReferenceMaker`].

```
INDEX="${DIR_INDEX}/Amel_HAv3.1.fasta"

cd ${DIR_INDEX}

conda activate gatk

gatk CreateSequenceDictionary -R Amel_HAv3.1.fasta -O Amel_HAv3.1.dict

cd ${DIR_VARIANTS}

for i in "${SRA[@]}"
do
   conda activate bcftools
   bcftools filter -e 'GT="het"' ${i}.vcf \
   | bcftools filter -i 'TYPE="snp"' - > ${i}_het_snp.vcf
   bgzip -c ${i}_het_snp.vcf > ${i}_het_snp.vcf.gz
   tabix -p vcf ${i}_het_snp.vcf.gz
   conda deactivate
   conda activate gatk
   gatk FastaAlternateReferenceMaker -R ${INDEX} -O ${DIR_ARG}/${i}.fasta -V ${i}_het_snp.vcf.gz
   conda deactivate
done

conda deactivate

cd ${DIR_ARG}

grep ">" SRR3037350.fasta | sed 's/>//g' > bad_headers.txt
grep ">" ${DIR_INDEX}/Amel_HAv3.1.fasta | sed 's/\s.*$//' | sed 's/>//g' > good_headers.txt
paste -d"\t" bad_headers.txt good_headers.txt > replace_headers.tsv

for i in "${SRA[@]}"
do
awk 'FNR==NR{  a[">"$1]=$2;next}$1 in a{  sub(/>/,">"a[$1]"|",$1)}1' \
replace_headers.tsv ${i}.fasta | sed 's/:.*//' > ${i}.fa
done
```

## Create tophat2 indices for each F0 genome

```
conda activate tophat2

for i in "${SRA[@]}"
do
bowtie2-build ${i}.fa ${i}
done

conda deactivate
```

# Process F1 RNA-seq libraries

### Sample metadata

|    SRA     | Cross |   Phenotype  |
| ---------- | ----- | ------------ |
| SRR3033262 |  875  | Reproductive |
| SRR3033263 |  875  | Reproductive |
| SRR3033264 |  875  | Reproductive |
| SRR3033249 |  875  | Sterile      |
| SRR3033250 |  875  | Sterile      |
| SRR3033251 |  875  | Sterile      |
| SRR3033268 |  888  | Reproductive |
| SRR3033269 |  888  | Reproductive |
| SRR3033270 |  888  | Reproductive |
| SRR3033256 |  888  | Sterile      |
| SRR3033257 |  888  | Sterile      |
| SRR3033258 |  888  | Sterile      |
| SRR3033265 |  882  | Reproductive |
| SRR3033266 |  882  | Reproductive |
| SRR3033267 |  882  | Reproductive |
| SRR3033252 |  882  | Sterile      |
| SRR3033253 |  882  | Sterile      |
| SRR3033254 |  882  | Sterile      |
| SRR3033255 |  882  | Sterile      |
| SRR3033271 |  894  | Reproductive |
| SRR3033272 |  894  | Reproductive |
| SRR3033273 |  894  | Reproductive |
| SRR3033259 |  894  | Sterile      |
| SRR3033260 |  894  | Sterile      |
| SRR3033261 |  894  | Sterile      |

## Retrieve F1 RNA-seq libraries

```
conda activate sra-tools

cd ${DIR_SRA}

RNA_SRA=("SRR3033249" "SRR3033250" "SRR3033251" "SRR3033252" \
         "SRR3033253" "SRR3033254" "SRR3033255" "SRR3033256" \
         "SRR3033257" "SRR3033258" "SRR3033259" "SRR3033260" \
         "SRR3033261" "SRR3033262" "SRR3033263" "SRR3033264" \
         "SRR3033265" "SRR3033266" "SRR3033267" "SRR3033268" \
         "SRR3033269" "SRR3033270" "SRR3033271" "SRR3033272" "SRR3033273")

for i in "${RNA_SRA[@]}"
do
  prefetch ${i}
  fasterq-dump -O ${DIR_SRA} ${i}
done

conda deactivate
```

## Trim adapters from F1 RNA-seq libraries

```
conda activate fastp

for i in "${RNA_SRA[@]}"
do
  fastp -w 8 -i ${i}_1.fastq -I ${i}_2.fastq \
  -o ${DIR_TRIM}/${i}_1.fastq -O ${DIR_TRIM}/${i}_2.fastq
done

conda deactivate
```

## Align F1 RNA-seq libraries to F0 genomes

Align [`tophat2`] F1 libraries (SRA accessions of lists `l875Q`, `l888Q`, `l882Q`, and `l894Q`) to respective F0 genomes and output coordinate-sorted BAM.

```
cd ${DIR_TRIM}

## SRR3037350 (875Q|E) vs SRR3037351 (875D|A)
l875Q=("SRR3033262" "SRR3033263" "SRR3033264" \
       "SRR3033249" "SRR3033250" "SRR3033251")

for i in "${l875Q[@]}"
do
tophat2 --library-type fr-firststrand --b2-very-sensitive -p 8 -o ${DIR_ALIGN}/875Q_${i} \
${DIR_ARG}/SRR3037350 ${i}_1.fastq ${i}_2.fastq

tophat2 --library-type fr-firststrand --b2-very-sensitive -p 8 -o ${DIR_ALIGN}/875D_${i} \
${DIR_ARG}/SRR3037351 ${i}_1.fastq ${i}_2.fastq
done


## SRR3037352 (888Q|A) vs SRR3037353 (888D|E)
l888Q=("SRR3033268" "SRR3033269" "SRR3033270" \
       "SRR3033256" "SRR3033257" "SRR3033258")

for i in "${l888Q[@]}"
do
tophat2 --library-type fr-firststrand --b2-very-sensitive -p 8 -o ${DIR_ALIGN}/888Q_${i} \
${DIR_ARG}/SRR3037352 ${i}_1.fastq ${i}_2.fastq

tophat2 --library-type fr-firststrand --b2-very-sensitive -p 8 -o ${DIR_ALIGN}/888D_${i} \
${DIR_ARG}/SRR3037353 ${i}_1.fastq ${i}_2.fastq
done


## SRR3037354 (882Q|E) vs SRR3037355 (882D|A)
l882Q=("SRR3033265" "SRR3033266" "SRR3033267" \
       "SRR3033252" "SRR3033253" "SRR3033254" "SRR3033255")

for i in "${l882Q[@]}"
do
tophat2 --library-type fr-firststrand --b2-very-sensitive -p 8 -o ${DIR_ALIGN}/882Q_${i} \
${DIR_ARG}/SRR3037354 ${i}_1.fastq ${i}_2.fastq

tophat2 --library-type fr-firststrand --b2-very-sensitive -p 8 -o ${DIR_ALIGN}/882D_${i} \
${DIR_ARG}/SRR3037355 ${i}_1.fastq ${i}_2.fastq
done


## SRR3037356 (894Q|A) vs SRR3037357 (894D|E)
l894Q=("SRR3033271" "SRR3033272" "SRR3033273" \
       "SRR3033259" "SRR3033260" "SRR3033261")

for i in "${l894Q[@]}"
do
tophat2 --library-type fr-firststrand --b2-very-sensitive -p 8 -o ${DIR_ALIGN}/894Q_${i} \
${DIR_ARG}/SRR3037356 ${i}_1.fastq ${i}_2.fastq

tophat2 --library-type fr-firststrand --b2-very-sensitive -p 8 -o ${DIR_ALIGN}/894D_${i} \
${DIR_ARG}/SRR3037357 ${i}_1.fastq ${i}_2.fastq
done

conda deactivate
```

# Count F1 reads over F0 variants

### Reciprocal cross colonies

| Block | Cross 1 | Cross 2 |
| ----- | ------- | ------- |
|   A   |   875   |   888   |
|   B   |   882   |   894   |

## Filter variants for analysis

### Remove concordant variants between parents

Find intersect between homozygous variants in F0 parents and report those that do not overlap [`bedtools intersect`].

```
cd ${DIR_VARIANTS}

conda activate bedtools

bedtools intersect -header -v -a SRR3037350_het_snp.vcf.gz \
-b SRR3037351_het_snp.vcf.gz > SRR3037350_outer.vcf

bedtools intersect -header -v -a SRR3037351_het_snp.vcf.gz \
-b SRR3037350_het_snp.vcf.gz > SRR3037351_outer.vcf

bedtools intersect -header -v -a SRR3037352_het_snp.vcf.gz \
-b SRR3037353_het_snp.vcf.gz > SRR3037352_outer.vcf

bedtools intersect -header -v -a SRR3037353_het_snp.vcf.gz \
-b SRR3037352_het_snp.vcf.gz > SRR3037353_outer.vcf


bedtools intersect -header -v -a SRR3037354_het_snp.vcf.gz \
-b SRR3037355_het_snp.vcf.gz > SRR3037354_outer.vcf

bedtools intersect -header -v -a SRR3037355_het_snp.vcf.gz \
-b SRR3037354_het_snp.vcf.gz > SRR3037355_outer.vcf

bedtools intersect -header -v -a SRR3037356_het_snp.vcf.gz \
-b SRR3037357_het_snp.vcf.gz > SRR3037356_outer.vcf

bedtools intersect -header -v -a SRR3037357_het_snp.vcf.gz \
-b SRR3037356_het_snp.vcf.gz > SRR3037357_outer.vcf


conda deactivate
```

1) Remove lines beginning with "#" from each VCF [`grep`].  
2) Merge discordant variants to one VCF per cross [`cat`].  
3) Compress VCFs [`bgzip`].

```
conda activate bcftools

grep -v '^#' SRR3037350_outer.vcf | cat SRR3037351_outer.vcf - > 875.vcf
bgzip -c 875.vcf > 875.vcf.gz

grep -v '^#' SRR3037352_outer.vcf | cat SRR3037353_outer.vcf - > 888.vcf
bgzip -c 888.vcf > 888.vcf.gz

grep -v '^#' SRR3037354_outer.vcf | cat SRR3037355_outer.vcf - > 882.vcf
bgzip -c 882.vcf > 882.vcf.gz

grep -v '^#' SRR3037356_outer.vcf | cat SRR3037357_outer.vcf - > 894.vcf
bgzip -c 894.vcf > 894.vcf.gz

conda deactivate
```

### Keep variants that are concordant between crosses (consensus variants)

Find intersect between homozygous variants between crosses for each block and report those that overlap [`bedtools intersect`].

```
conda activate bedtools

bedtools intersect -header -u -a 875.vcf.gz -b 888.vcf.gz > 875_888_aSet.vcf
bedtools intersect -header -u -a 882.vcf.gz -b 894.vcf.gz > 882_894_aSet.vcf
```

1) Remove lines beginning with "#" from each VCF [`grep`].  
2) Keep only columns 1 (chromosome) and 2 (start coordinate of SNP), duplicate column 2 as column 3 (stop coordinate of SNP) [`awk`].  
3) Add column 4 (name) as sequential list of numbers [`awk`].  
4) Add "snp_" to beginning of each line at column 4 [`awk`].  
5) Find variants shared between blocks and report those that overlap [`bedtools intersect`].

```
grep -v '^#' 875_888_aSet.vcf | awk -v OFS="\t" '{print $1, $2, $2}' \
 | awk -v OFS="\t" '$4=(FNR FS $4)' \
 | awk -v OFS="\t" '{print $1, $2, $3, "snp_"$4}' > 875_888_aSet.bed

grep -v '^#' 882_894_aSet.vcf | awk -v OFS="\t" '{print $1, $2, $2}' \
  | awk -v OFS="\t" '$4=(FNR FS $4)' \
  | awk -v OFS="\t" '{print $1, $2, $3, "snp_"$4}' > 882_894_aSet.bed

bedtools intersect -u -a 875_888_aSet.bed -b 882_894_aSet.bed > consensus_aSet.bed
```

## Generate BED file of variants intersecting genes

```
awk '$3 == "gene" { print $0 }' Amel_HAv3.1.gff | awk -v OFS="\t" '{print $1, $2, $3, $4, $5, $6, $7, $8}' > Amel_HAv3.1_genes.txt
awk '$3 == "gene" { print $0 }' Amel_HAv3.1.gff | awk '{print $9}' | grep -Po 'GeneID[^\s]*' \
| cut -d':' -f2 | cut -d';' -f1 | cut -d',' -f1 | sed -e 's/^/LOC/' > Amel_HAv3.1_geneIDs.txt
paste -d'\t' Amel_HAv3.1_genes.txt Amel_HAv3.1_geneIDs.txt > Amel_HAv3.1_genes.gff3

bedtools intersect -wb -a Amel_HAv3.1_genes.gff3 -b ${DIR_VARIANTS}/consensus_aSet.bed \
| awk -v OFS="\t" '{print $10, $11, $12, $13 ":" $9, $6, $7}' \
| grep -v '^NC_001566.1' | sort -k1,1V -k2,2n > ${DIR_VARIANTS}/SNPs_for_analysis.bed

conda deactivate
```

## Compute read coverage at each SNP:gene

Count [`bedtools intersect`] strand-wise reads of F1 libraries (SRA accessions of lists `l875Q`, `l888Q`, `l882Q`, and `l894Q`) aligned to respective F0 genomes at each SNP-gene, requiring 0 mismatches.

```
sort --parallel=8 -k1,1 -k2,2n variants_for_analysis.bed > snps_for_analysis_sorted.bed


for i in "${l875Q[@]}"
do
conda activate bamtools
bamtools filter -tag XM:0 -in ${DIR_ALIGN}/875Q_${i}/accepted_hits.bam \
-out ${DIR_ALIGN}/875Q_${i}/nomm_hits.bam
bamtools filter -tag XM:0 -in ${DIR_ALIGN}/875D_${i}/accepted_hits.bam \
-out ${DIR_ALIGN}/875D_${i}/nomm_hits.bam

conda deactivate
conda activate bedtools

bedtools bamtobed -i ${DIR_ALIGN}/875Q_${i}/nomm_hits.bam \
> ${DIR_ALIGN}/875Q_${i}/nomm_hits.bed

sort --parallel=10 -k1,1 -k2,2n ${DIR_ALIGN}/875Q_${i}/nomm_hits.bed \
> ${DIR__SORT}/875Q_${i}_nomm_hits_sorted.bed

bedtools intersect -S -sorted -c -a ${DIR_VARIANTS}/snps_for_analysis_sorted.bed \
-b ${DIR_SORT}/875Q_${i}_nomm_hits_sorted.bed \
> ${DIR_COUNTS}/875Q_${i}.txt

bedtools bamtobed -i ${DIR_ALIGN}/875D_${i}/nomm_hits.bam \
> ${DIR_ALIGN}/875D_${i}/nomm_hits.bed

sort --parallel=10 -k1,1 -k2,2n ${DIR_ALIGN}/875D_${i}/nomm_hits.bed \
> ${DIR_SORT}/875D_${i}_nomm_hits_sorted.bed

bedtools intersect -S -sorted -c -a ${DIR_VARIANTS}/snps_for_analysis_sorted.bed \
-b ${DIR_SORT}/875D_${i}_nomm_hits_sorted.bed \
> ${DIR_COUNTS}/875D_${i}.txt

conda deactivate
done


for i in "${l888Q[@]}"
do
conda activate bamtools
bamtools filter -tag XM:0 -in ${DIR_ALIGN}/888Q_${i}/accepted_hits.bam \
-out ${DIR_ALIGN}/888Q_${i}/nomm_hits.bam
bamtools filter -tag XM:0 -in ${DIR_ALIGN}/888D_${i}/accepted_hits.bam \
-out ${DIR_ALIGN}/888D_${i}/nomm_hits.bam

conda deactivate
conda activate bedtools

bedtools bamtobed -i ${DIR_ALIGN}/888Q_${i}/nomm_hits.bam \
> ${DIR_ALIGN}/888Q_${i}/nomm_hits.bed

sort --parallel=10 -k1,1 -k2,2n ${DIR_ALIGN}/888Q_${i}/nomm_hits.bed \
> ${DIR_SORT}/888Q_${i}_nomm_hits_sorted.bed

bedtools intersect -S -sorted -c -a ${DIR_VARIANTS}/snps_for_analysis_sorted.bed \
-b ${DIR_SORT}/888Q_${i}_nomm_hits_sorted.bed \
> ${DIR_COUNTS}/888Q_${i}.txt

bedtools bamtobed -i ${DIR_ALIGN}/888D_${i}/nomm_hits.bam \
> ${DIR_ALIGN}/888D_${i}/nomm_hits.bed

sort --parallel=10 -k1,1 -k2,2n ${DIR_ALIGN}/888D_${i}/nomm_hits.bed \
> ${DIR_SORT}/888D_${i}_nomm_hits_sorted.bed

bedtools intersect -S -sorted -c -a ${DIR_VARIANTS}/snps_for_analysis_sorted.bed \
-b ${DIR_SORT}/888D_${i}_nomm_hits_sorted.bed \
> ${DIR_COUNTS}/888D_${i}.txt

conda deactivate
done


for i in "${l882Q[@]}"
do
conda activate bamtools
bamtools filter -tag XM:0 -in ${DIR_ALIGN}/882Q_${i}/accepted_hits.bam \
-out ${DIR_ALIGN}/882Q_${i}/nomm_hits.bam
bamtools filter -tag XM:0 -in ${DIR_ALIGN}/882D_${i}/accepted_hits.bam \
-out ${DIR_ALIGN}/882D_${i}/nomm_hits.bam

conda deactivate
conda activate bedtools

bedtools bamtobed -i ${DIR_ALIGN}/882Q_${i}/nomm_hits.bam \
> ${DIR_ALIGN}/882Q_${i}/nomm_hits.bed

sort --parallel=10 -k1,1 -k2,2n ${DIR_ALIGN}/882Q_${i}/nomm_hits.bed \
> ${DIR_SORT}/882Q_${i}_nomm_hits_sorted.bed

bedtools intersect -S -sorted -c -a ${DIR_VARIANTS}/snps_for_analysis_sorted.bed \
-b ${DIR_SORT}/882Q_${i}_nomm_hits_sorted.bed \
> ${DIR_COUNTS}/882Q_${i}.txt

bedtools bamtobed -i ${DIR_ALIGN}/882D_${i}/nomm_hits.bam \
> ${DIR_ALIGN}/882D_${i}/nomm_hits.bed

sort --parallel=10 -k1,1 -k2,2n ${DIR_ALIGN}/882D_${i}/nomm_hits.bed \
> ${DIR_SORT}/882D_${i}_nomm_hits_sorted.bed

bedtools intersect -S -sorted -c -a ${DIR_VARIANTS}/snps_for_analysis_sorted.bed \
-b ${DIR_SORT}/882D_${i}_nomm_hits_sorted.bed \
> ${DIR_COUNTS}/882D_${i}.txt

conda deactivate
done


for i in "${l894Q[@]}"
do
conda activate bamtools
bamtools filter -tag XM:0 -in ${DIR_ALIGN}/894Q_${i}/accepted_hits.bam \
-out ${DIR_ALIGN}/894Q_${i}/nomm_hits.bam
bamtools filter -tag XM:0 -in ${DIR_ALIGN}/894D_${i}/accepted_hits.bam \
-out ${DIR_ALIGN}/894D_${i}/nomm_hits.bam

conda deactivate
conda activate bedtools

bedtools bamtobed -i ${DIR_ALIGN}/894Q_${i}/nomm_hits.bam \
> ${DIR_ALIGN}/894Q_${i}/nomm_hits.bed

sort --parallel=10 -k1,1 -k2,2n ${DIR_ALIGN}/894Q_${i}/nomm_hits.bed \
> ${DIR_SORT}/894Q_${i}_nomm_hits_sorted.bed

bedtools intersect -S -sorted -c -a ${DIR_VARIANTS}/snps_for_analysis_sorted.bed \
-b ${DIR_SORT}/894Q_${i}_nomm_hits_sorted.bed \
> ${DIR_COUNTS}/894Q_${i}.txt

bedtools bamtobed -i ${DIR_ALIGN}/894D_${i}/nomm_hits.bam \
> ${DIR_ALIGN}/894D_${i}/nomm_hits.bed

sort --parallel=10 -k1,1 -k2,2n ${DIR_ALIGN}/894D_${i}/nomm_hits.bed \
> ${DIR_SORT}/894D_${i}_nomm_hits_sorted.bed

bedtools intersect -S -sorted -c -a ${DIR_VARIANTS}/snps_for_analysis_sorted.bed \
-b ${DIR_SORT}/894D_${i}_nomm_hits_sorted.bed \
> ${DIR_COUNTS}/894D_${i}.txt

conda deactivate
done
```
