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

##### star
```
conda create --name star
conda install -c bioconda -n star star
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
INDEX_GTF="/storage/home/stb5321/scratch/galbraith/index/Amel_HAv3.1.gtf"
DIR_RNA="/storage/home/stb5321/scratch/galbraith/STAR_variants"
DIR_COUNTS="/storage/home/stb5321/scratch/galbraith/counts"
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

wget -O Amel_HAv3.1.gff.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/254/395/GCF_003254395.2_Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_genomic.gff.gz
gunzip Amel_HAv3.1.gff.gz

wget -O Amel_HAv3.1.gtf.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/254/395/GCF_003254395.2_Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_genomic.gtf.gz
gunzip Amel_HAv3.1.gtf.gz

conda activate bwa

bwa index Amel_HAv3.1.fna

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
  bwa mem -t 8 ${DIR_INDEX}/Amel_HAv3.1.fna ${i}_1.fastq ${i}_2.fastq > ${DIR_ALIGN}/${i}.sam
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

samtools faidx Amel_HAv3.1.fna

conda deactivate

cd ${DIR_ALIGN}

conda activate freebayes

for i in "${DIPLOID[@]}"
do
  freebayes -f ${DIR_INDEX}/Amel_HAv3.1.fna ${i}.bam > ${DIR_VARIANTS}/${i}.vcf
done

for i in "${HAPLOID[@]}"
do
  freebayes -f ${DIR_INDEX}/Amel_HAv3.1.fna ${i}.bam -p 1 > ${DIR_VARIANTS}/${i}.vcf
done

conda deactivate
```

## Generate F0 reference genomes

1) Remove heterozygous variants and indels [`bcftools filter`] from VCFs. NOTE: variants, MNPs, and complex variants are retained.
2) Compress resultant VCFs and index [`bgzip`, `tabix`].  
3) Integrate homozygous variants into Amel_HAv3.1 for each F0 library, separately [`bcftools consensus`].

```
conda activate bcftools

cd ${DIR_VARIANTS}

for i in "${SRA[@]}"
do
  bcftools filter -e 'GT="het" | 'TYPE="indel"' - > ${i}_typefilter.vcf
  bgzip -c ${i}_typefilter.vcf > ${i}_typefilter.vcf.gz
  tabix -p vcf ${i}_typefilter.vcf.gz
  bcftools consensus -f ${INDEX} -o ${DIR_ARG}/${i}.fna ${i}_typefilter.vcf.gz
done

cd ${DIR_ARG}

for i in "${SRA[@]}"
do
   sed '/^>/ s/ .*//' ${i}.fa > ${i}.fasta
   rm ${i}.fa
done

conda deactivate
```

## Create STAR indices for each F0 genome

```
cd ${DIR_ARG}

conda activate star

for i in "${SRA[@]}"
do
  mkdir ${i}_STAR
  STAR --runThreadN 8 \
  --runMode genomeGenerate \
  --genomeDir ${i}_STAR \
  --genomeFastaFiles ${i}.fna \
  --sjdbGTFfile ${INDEX_GTF} \
  --sjdbOverhang 99 \
  --genomeSAindexNbases 12
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

Align [`STAR`] F1 libraries (SRA accessions of lists `l875Q`, `l888Q`, `l882Q`, and `l894Q`) to respective F0 genomes, allowing for 0 mismatches and output coordinate-sorted BAM.

```
cd ${DIR_TRIM}

conda activate star


l875Q=("SRR3033262" "SRR3033263" "SRR3033264" \
      "SRR3033249" "SRR3033250" "SRR3033251")

for i in "${l875Q[@]}"
do
STAR --genomeDir ${DIR_ARG}/SRR3037350_STAR \
     --runThreadN 8 \
     --outFilterMismatchNmax 0 \
     --readFilesIn ${i}_1.fastq ${i}_2.fastq \
     --outFileNamePrefix ${DIR_ALIGN}/875Q_${i} \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMunmapped None \
     --outSAMattributes Standard
     --alignEndsType EndToEnd

STAR --genomeDir ${DIR_ARG}/SRR3037351_STAR \
     --runThreadN 8 \
     --outFilterMismatchNmax 0 \
     --readFilesIn ${i}_1.fastq ${i}_2.fastq \
     --outFileNamePrefix ${DIR_ALIGN}/875D_${i} \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMunmapped None \
     --outSAMattributes Standard
     --alignEndsType EndToEnd
done


l888Q=("SRR3033268" "SRR3033269" "SRR3033270" \
      "SRR3033256" "SRR3033257" "SRR3033258")

for i in "${l888Q[@]}"
do
STAR --genomeDir ${DIR_ARG}/SRR3037352_STAR \
     --runThreadN 8 \
     --outFilterMismatchNmax 0 \
     --readFilesIn ${i}_1.fastq ${i}_2.fastq \
     --outFileNamePrefix ${DIR_ALIGN}/888Q_${i} \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMunmapped None \
     --outSAMattributes Standard
     --alignEndsType EndToEnd

STAR --genomeDir ${DIR_ARG}/SRR3037353_STAR \
     --runThreadN 8 \
     --outFilterMismatchNmax 0 \
     --readFilesIn ${i}_1.fastq ${i}_2.fastq \
     --outFileNamePrefix ${DIR_ALIGN}/888D_${i} \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMunmapped None \
     --outSAMattributes Standard
     --alignEndsType EndToEnd
done


l882Q=("SRR3033265" "SRR3033266" "SRR3033267" \
      "SRR3033252" "SRR3033253" "SRR3033254" "SRR3033255")

for i in "${l882Q[@]}"
do
STAR --genomeDir ${DIR_ARG}/SRR3037354_STAR \
     --runThreadN 8 \
     --outFilterMismatchNmax 0 \
     --readFilesIn ${i}_1.fastq ${i}_2.fastq \
     --outFileNamePrefix ${DIR_ALIGN}/882Q_${i} \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMunmapped None \
     --outSAMattributes Standard
     --alignEndsType EndToEnd

STAR --genomeDir ${DIR_ARG}/SRR3037355_STAR \
     --runThreadN 8 \
     --outFilterMismatchNmax 0 \
     --readFilesIn ${i}_1.fastq ${i}_2.fastq \
     --outFileNamePrefix ${DIR_ALIGN}/882D_${i} \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMunmapped None \
     --outSAMattributes Standard
     --alignEndsType EndToEnd
done


l894Q=("SRR3033271" "SRR3033272" "SRR3033273" \
      "SRR3033259" "SRR3033260" "SRR3033261")

for i in "${l894Q[@]}"
do
STAR --genomeDir ${DIR_ARG}/SRR3037356_STAR \
     --runThreadN 8 \
     --outFilterMismatchNmax 0 \
     --readFilesIn ${i}_1.fastq ${i}_2.fastq \
     --outFileNamePrefix ${DIR_ALIGN}/894Q_${i} \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMunmapped None \
     --outSAMattributes Standard
     --alignEndsType EndToEnd

STAR --genomeDir ${DIR_ARG}/SRR3037357_STAR \
     --runThreadN 8 \
     --outFilterMismatchNmax 0 \
     --readFilesIn ${i}_1.fastq ${i}_2.fastq \
     --outFileNamePrefix ${DIR_ALIGN}/894D_${i} \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMunmapped None \
     --outSAMattributes Standard
     --alignEndsType EndToEnd
done

conda deactivate
```

# Count F1 reads over F0 variants

### Reciprocal cross colonies

| Block | Cross 1 | Cross 2 |
| ----- | ------- | ------- |
|   A   |   875   |   888   |
|   B   |   882   |   894   |

## Define directory variables

```
DIR_RNA="/storage/home/stb5321/scratch/galbraith/STAR_variants"
DIR_COUNTS="/storage/home/stb5321/scratch/galbraith/counts"
```

## Filter variants for analysis

### Remove concordant variants between parents

Find intersect between homozygous variants in F0 parents and report those that do not overlap [`bedtools intersect`].

```
cd ${DIR_VARIANTS}

conda activate bedtools


bedtools intersect -header -v -a SRR3037350_typefilter.vcf.gz \
-b SRR3037351_typefilter.vcf.gz > SRR3037350_outer.vcf

bedtools intersect -header -v -a SRR3037351_typefilter.vcf.gz \
-b SRR3037350_typefilter.vcf.gz > SRR3037351_outer.vcf

bedtools intersect -header -v -a SRR3037352_typefilter.vcf.gz \
-b SRR3037353_typefilter.vcf.gz > SRR3037352_outer.vcf

bedtools intersect -header -v -a SRR3037353_typefilter.vcf.gz \
-b SRR3037352_typefilter.vcf.gz > SRR3037353_outer.vcf


bedtools intersect -header -v -a SRR3037354_typefilter.vcf.gz \
-b SRR3037355_typefilter.vcf.gz > SRR3037354_outer.vcf

bedtools intersect -header -v -a SRR3037355_typefilter.vcf.gz \
-b SRR3037354_typefilter.vcf.gz > SRR3037355_outer.vcf

bedtools intersect -header -v -a SRR3037356_typefilter.vcf.gz \
-b SRR3037357_typefilter.vcf.gz > SRR3037356_outer.vcf

bedtools intersect -header -v -a SRR3037357_typefilter.vcf.gz \
-b SRR3037356_typefilter.vcf.gz > SRR3037357_outer.vcf


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
awk '$3 == "gene" { print $0 }' Amel_HAv3.1_OGSv3.2.gff3 | sed 's/;//g' > Amel_HAv3.1_OGSv3.2_genes.gff3

bedtools intersect -wb -a Amel_HAv3.1_OGSv3.2_genes.gff3 -b ${DIR_VARIANTS}/consensus_aSet.bed \
| awk -v OFS="\t" '{print $10, $11, $12, $13 ":" $9, $6, $7}' \
| grep -v '^NC_001566.1' | sort -k1,1V -k2,2n > ${DIR_VARIANTS}/SNPs_for_analysis.bed
```

## Compute strand-wise read coverage at each SNP:gene

Count [`bedtools intersect`] reads of F1 libraries (SRA accessions of lists `l875Q`, `l888Q`, `l882Q`, and `l894Q`) aligned to respective F0 genomes at each SNP-gene, accounting for strandedness [`-S` since gene transcripts are antisense.]

```
sort --parallel=8 -k1,1 -k2,2n variants_for_analysis.bed > variants_for_analysis_sorted.bed


for i in "${l875Q[@]}"
do
bedtools bamtobed -i ${DIR_RNA}/875Q_${i}Aligned.sortedByCoord.out.bam \
> ${DIR_RNA}/875Q_${i}Aligned.sortedByCoord.out.bed

sort --parallel=8 -k1,1 -k2,2n ${DIR_RNA}/875Q_${i}Aligned.sortedByCoord.out.bed \
> ${DIR_RNA}/875Q_${i}Aligned.sortedByCoord.out.sorted.bed

bedtools intersect -sorted -S -c -a variants_for_analysis_sorted.bed \
-b ${DIR_RNA}/875Q_${i}Aligned.sortedByCoord.out.sorted.bed \
> ${DIR_COUNTS}/875Q_${i}.txt

bedtools bamtobed -i ${DIR_RNA}/875D_${i}Aligned.sortedByCoord.out.bam \
> ${DIR_RNA}/875D_${i}Aligned.sortedByCoord.out.bed

sort --parallel=8 -k1,1 -k2,2n ${DIR_RNA}/875D_${i}Aligned.sortedByCoord.out.bed \
> ${DIR_RNA}/875D_${i}Aligned.sortedByCoord.out.sorted.bed

bedtools intersect -sorted -S -c -a variants_for_analysis_sorted.bed \
-b ${DIR_RNA}/875D_${i}Aligned.sortedByCoord.out.sorted.bed \
> ${DIR_COUNTS}/875D_${i}.txt
done


for i in "${l888Q[@]}"
do
bedtools bamtobed -i ${DIR_RNA}/888Q_${i}Aligned.sortedByCoord.out.bam \
> ${DIR_RNA}/888Q_${i}Aligned.sortedByCoord.out.bed

sort --parallel=8 -k1,1 -k2,2n ${DIR_RNA}/888Q_${i}Aligned.sortedByCoord.out.bed \
> ${DIR_RNA}/888Q_${i}Aligned.sortedByCoord.out.sorted.bed

bedtools intersect -sorted -S -c -a variants_for_analysis_sorted.bed \
-b ${DIR_RNA}/888Q_${i}Aligned.sortedByCoord.out.sorted.bed \
> ${DIR_COUNTS}/888Q_${i}.txt

bedtools bamtobed -i ${DIR_RNA}/888D_${i}Aligned.sortedByCoord.out.bam \
> ${DIR_RNA}/888D_${i}Aligned.sortedByCoord.out.bed

sort --parallel=8 -k1,1 -k2,2n ${DIR_RNA}/888D_${i}Aligned.sortedByCoord.out.bed \
> ${DIR_RNA}/888D_${i}Aligned.sortedByCoord.out.sorted.bed

bedtools intersect -sorted -S -c -a variants_for_analysis_sorted.bed \
-b ${DIR_RNA}/888D_${i}Aligned.sortedByCoord.out.sorted.bed \
> ${DIR_COUNTS}/888D_${i}.txt
done


for i in "${l882Q[@]}"
do
bedtools bamtobed -i ${DIR_RNA}/882Q_${i}Aligned.sortedByCoord.out.bam \
> ${DIR_RNA}/882Q_${i}Aligned.sortedByCoord.out.bed

sort --parallel=8 -k1,1 -k2,2n ${DIR_RNA}/882Q_${i}Aligned.sortedByCoord.out.bed \
> ${DIR_RNA}/882Q_${i}Aligned.sortedByCoord.out.sorted.bed

bedtools intersect -sorted -S -c -a variants_for_analysis_sorted.bed \
-b ${DIR_RNA}/882Q_${i}Aligned.sortedByCoord.out.sorted.bed \
> ${DIR_COUNTS}/882Q_${i}.txt

bedtools bamtobed -i ${DIR_RNA}/882D_${i}Aligned.sortedByCoord.out.bam \
> ${DIR_RNA}/882D_${i}Aligned.sortedByCoord.out.bed

sort --parallel=8 -k1,1 -k2,2n ${DIR_RNA}/882D_${i}Aligned.sortedByCoord.out.bed \
> ${DIR_RNA}/882D_${i}Aligned.sortedByCoord.out.sorted.bed

bedtools intersect -sorted -S -c -a variants_for_analysis_sorted.bed \
-b ${DIR_RNA}/882D_${i}Aligned.sortedByCoord.out.sorted.bed \
> ${DIR_COUNTS}/882D_${i}.txt
done


for i in "${l894Q[@]}"
do
bedtools bamtobed -i ${DIR_RNA}/894Q_${i}Aligned.sortedByCoord.out.bam \
> ${DIR_RNA}/894Q_${i}Aligned.sortedByCoord.out.bed

sort --parallel=8 -k1,1 -k2,2n ${DIR_RNA}/894Q_${i}Aligned.sortedByCoord.out.bed \
> ${DIR_RNA}/894Q_${i}Aligned.sortedByCoord.out.sorted.bed

bedtools intersect -sorted -S -c -a variants_for_analysis_sorted.bed \
-b ${DIR_RNA}/894Q_${i}Aligned.sortedByCoord.out.sorted.bed \
> ${DIR_COUNTS}/894Q_${i}.txt

bedtools bamtobed -i ${DIR_RNA}/894D_${i}Aligned.sortedByCoord.out.bam \
> ${DIR_RNA}/894D_${i}Aligned.sortedByCoord.out.bed

sort --parallel=8 -k1,1 -k2,2n ${DIR_RNA}/894D_${i}Aligned.sortedByCoord.out.bed \
> ${DIR_RNA}/894D_${i}Aligned.sortedByCoord.out.sorted.bed

bedtools intersect -sorted -S -c -a variants_for_analysis_sorted.bed \
-b ${DIR_RNA}/894D_${i}Aligned.sortedByCoord.out.sorted.bed \
> ${DIR_COUNTS}/894D_${i}.txt
done


conda deactivate
```
