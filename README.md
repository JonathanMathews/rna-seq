# RNA-Seq Pipeline

**RNA-Seq pipeline: Learning Nextflow and introduction to Bioinformatics**

## Step 0: Download data

`docker run -t --rm -v $PWD:/output:rw -w /output ncbi/sra-tools fasterq-dump -e 2 -p SRR6875349`

`docker run -t --rm -v $PWD:/output:rw -w /output ncbi/sra-tools fasterq-dump -e 2 -p SRR6875339`

## Step 1: Trim data

**Trim Galore**

- Wrapper around Cutadapt and FastQC
- User friendly and most widely used

**fastp**

- Faster version of trimming
- All-in-one preprocessing and QC for FastQ short reads
- *For long reads used fastplong*

`docker run -v $(pwd):/data biocontainers/fastp:v0.20.1_cv1 fastp -i SRR6875349_1.fq.gz -I SRR6875349_2.fq.gz -o SRR6875349.trim.R1.fq.gz -O SRR6875349.trim.R2.fq.gz -w 4 -h SRR6875349.trim_report.html -j SRR6875349.trim_report.json -q 25 -l 30 --length_limit 0 -g -x`

`docker run -v $(pwd):/data biocontainers/fastp:v0.20.1_cv1 fastp -i SRR6875339_1.fq.gz -I SRR6875339_2.fq.gz -o SRR6875339.trim.R1.fq.gz -O SRR6875339.trim.R2.fq.gz -w 4 -h SRR6875339.trim_report.html -j SRR6875339.trim_report.json -q 25 -l 30 --length_limit 0 -g -x`

- i: Read 1 (positive strand) input
- I: Read 2 (negative strand) input
- o: Read 1 trimmed output
- O: Read 2 trimmed output
- w: Number of processing units used for job (multithreading)
- h: QC report of trimming in HTML format
- j: QC report of trimming in JSON format
- q: Minimum QC value
- l: Mininum length to be considered a valid read
- length_limit: Maximum length to be considered a valid read (0 means no limit)
- g: Poly G tail trimming
- x: Poly N tail trimming that occurs before the Poly G tail

## Step 2: Align

**HISAT2**

- Useful for larger datasets
- Needs less powerful hardware needed
- Faster and lower memory usage

**STAR**

*Use most of the time*

- More accurate and robust
- Requires more memory and more powerful hardware
- Slower

### Step 2a: Generate genomic indices

*STAR Command*

`docker run -v $(pwd):/data mgibio/star:2.7.0f STAR --runMode genomeGenerate --genomeDir /data/star_index --genomeFastaFiles /data/genome.fa --sjdbGTFfile /data/annotations.gtf --sjdbOverhang ReadLength-1 --runThreadN 4`

*HISAT2 Command*

`docker run -v $(pwd):/data dceoy/hisat2 hisat2-build -p 4 /data/genome.fa /data/genome`

- p: Number of processing units used for job (multithreading)
- genome.fa: FASTA file containing genome
- genome: Base name for genomic index files

### Step 2b: Align reads to genome

*STAR Command*

```
docker run -v $(pwd):/data mgibio/star:2.7.0f STAR --runMode alignReads --runThreadN 4 --genomeDir /data/star_index/ --readFilesIn /data/SRR6875349.trim.R1.fq.gz SRR6875349.trim.R2.fq.gz --readFilesCommand zcat --twopassMode None --genomeLoad  NoSharedMemory --outTmpDir tmp --outReadsUnmapped Fastx --outSAMreadID SRR6875349 --outSAMtype BAM Unsorted --outSAMattributes NH MD --chimOutType Junctions SeparateSAMold --chimOutJunctionFormat 1 --chimSegmentReadGapMax 3 --chimJunctionOverhangMin 12 --chimSegmentMin 12

mv SJ.out.tab SRR6875349.js.txt

sort -k1,1 -k2,2n Chimeric.out.junction > SRR6875349.chimera.txt

mv Unmapped.out.mate1 SRR6875349.unmapped.R1.fq

if [[ -f Unmapped.out.mate2 ]] ; then
    mv Unmapped.out.mate2 SRR6875349.unmapped.R2.fq
fi

mv Log.final.out SRR6875349.alignerout.txt

mv Aligned.out.bam out.bam

docker run -v $(pwd):/data biocontainers/samtools:v1.9-4-deb_cv1 samtools sort -@ 4 -O BAM -o SRR6875349.bam out.bam

docker run -v $(pwd):/data biocontainers/samtools:v1.9-4-deb_cv1 samtools index -@ 4 SRR6875349.bam
```

*HISAT2 Command*

```
docker run -v $(pwd):/data nanozoo/hisat2:2.2.1.commit7e01700--5e923e8 hisat2 -p 4 --rg-id SRR6875349 --rg LB:tx --rg PL:illumina --rg PU:barcode --rg SM:SRR6875349 --dta -x hisat2_index/genome -1 SRR6875349.trim.R1.fq.gz -2 SRR6875349.trim.R2.fq.gz -S out.sam --summary-file SRR6875349.alignerout.txt

docker run -v $(pwd):/data nanozoo/hisat2:2.2.1.commit7e01700--5e923e8 samtools view -1 --threads 4 -o out.bam out.sam

docker run -v $(pwd):/data nanozoo/hisat2:2.2.1.commit7e01700--5e923e8 samtools sort -@ 4 -O BAM -o SRR6875349.bam out.bam

docker run -v $(pwd):/data nanozoo/hisat2:2.2.1.commit7e01700--5e923e8 samtools index -@ 4 SRR6875349.bam
```

## Step 3: Feature Counts

`docker run -v $(pwd):/data pegi3s/feature-counts featureCounts -s 0 -M --fraction -J --ignoreDup -T 4 -p -g gene_id -a /data/gencode.v48.annotation.gtf -o /data/SRR6875349.cts.txt /data/SRR6875349.bam`

- s: Standedness of RNA data; Mostly will use unstranded (reads count towards + and - strands); 0 = unstranded, 1 = stranded, 2 = reversely stranded
- M: Multi-mapping reads are counted (reads can map to multiple features)
- fraction: Adds fractional counts for multi-mapped reads (1/x where x = total alignments reported for a specific read)
- J: Count reads aligned to exon-exon junctions
- ignoreDup: Ignores duplicate reads (due to amplification)
- T: Number of processing units used (multithreading)
- p: Specifies input data is paired-end
- g: Feature to group by (from GTF annotation file)
- a: GTF annotation file
- o: Output file containing counts

### Step 3a: Transcript Counts

```
mkdir -p SRR6875349_stringtie

grep '#\|transcript_id' gencode.v48.annotation.gtf | grep -v 'transcript_id \"\"' > transcripts.gtf

docker run -v $(pwd):/data mgibio/stringtie:v2.2.1-focal stringtie /data/SRR6875349.bam -p 4 -G /data/transcripts.gtf -B -e -o /data/SRR6875349.denovo.gtf -A /data/SRR6875349.fpkm.txt

mv *.ctab SRR6875349_stringtie

mv SRR6875349.denovo.gtf SRR6875349_stringtie

mv SRR6875349.fpkm.txt SRR6875349_stringtie
```

## Step 4: R Script for statistical analysis of feature counts

This will perform DEG analysis and output relevant figures and tables.

**See `nextflow.config` to setup filesytem for running this pipeline.**