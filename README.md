# Waltz

A fast, efficient bam pileup and application modules based on it, like coverage metrics, genotyping, signature finding etc.

This software was developed at the Innovation Lab, Center for Molecular Oncology, Memorial Sloan Kettering Cancer Center.
<br/>


Waltz has 2 main modules:  
1. **Bam metrics**: Generate various useful metrics for a given bam file
2. **Genotyping**: Determine the fragment count and allele fraction of given mutations in given bam file


## Java
Java 1.8 or above is required.

## Dependencies (bundled with the release jar)

1. BioinfoUtils
2. HTSJDK
3. Google Guava
4. Apache Commons IO

<br/>

### 1. Bam Metrics

#### Generate bam level metrics

java -server -Xms4g -Xmx4g -cp Waltz.jar org.mskcc.juber.waltz.countreads.CountReads bam-file coverageThreshold canonical-transcripts-bed-file intervals-bed-file

where  
coverageThreshold is the average coverage above which a contiguous region should be considered covered (suggested value: 5)  
canonical-transcripts-bed-file is the bed file with all exons in across the genomes (included above)  
intervals-bed-file is the bed file of chosen genomic intervals  


This produces 3 files:  
.covered-regions: regions of contiguous coverage, annotated with canonical transcripts. Useful for checking what regions are actually covered in the bam file. Columns: chr, start, end, length, average total coverage in the contiguous region.

.read-counts: bam-level stats. Columns: bam file name, total reads, unmapped reads, total mapped reads, unique mapped reads, duplicate fraction, total on-target reads, unique on-target reads, total on-target rate, unique on-target rate

.fragment-sizes: fragment size distribution. Columns: fragment-size, total frequency, unique frequency

 
#### Generate metrics specific to given genomic regions

java -server -Xms4g -Xmx4g -cp Waltz.jar org.mskcc.juber.waltz.Waltz PileupMetrics mappinngQualityThreshold bam-file reference-fasta intervals-bed-file

This produces 4 different files:
-pileup.txt: per-position fragment count for different alleles. Columns: chr, position, ref, depth (including N's), fragment counts for A, C, G, T, insertions, deletions, soft clip start, soft clip end, hard clip start, hard clip end

-pileup-without-duplicates.txt: similar to above but only unique fragments are counted

-intervals.txt: stats per genomic interval. Columns: chr, start, end, interval name, interval length, peak coverage, average coverage, GC fraction, number of fragments mapped

-intervals-without-duplicates.txt: similar to above but only unique fragments are considered


#### Collect metrics across samples

Run aggregate-bam-metrics.sh script in the folder where the above output files are present to collect metrics across samples.

This produces 3 main files with self-explanatory headers.
read-counts.txt: collection of metrics from *.read-counts files

waltz-coverage.txt: per sample coverage calculated across chosen genomic intervals

fragment-sizes.txt: fragment size distributions for all samples

 


### 2. Genotyping

java -server -Xms4g -Xmx4g -cp Waltz.jar org.mskcc.juber.waltz.Waltz Genotyping mappinngQualityThreshold bam-file reference-fasta intervals-bed-file mutations-maf-file

where
mutations-maf-file is a file in maf format specifying the mutations to be profiled in the given bam. Required fields are Chromosome, Start_Position, Variant_Type, Reference_Allele and Tumor_Seq_Allele2

This will produce a -genotypes.maf file with 4 addtional columns at the end: Waltz_total_t_depth, Waltz_total_t_alt_count, Waltz_MD_t_depth and Waltz_MD_t_alt_count. All sample-specific columns will be made empty while all the mutation-specific information will be retained. Tumor_Sample_Barcode will contain the name of the sample being genotyped.

#### Collect genotypes across multiple samples

Run aggregate-genotypes.sh script in the folder where the -genotypes.maf files are present to collect genotyping information across multiple samples. The output is a genotypes.maf file. 





















