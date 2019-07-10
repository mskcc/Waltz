# Waltz

A fast, efficient bam pileup and application modules based on it, like coverage metrics, genotyping, signature finding etc.

This software was developed at the Innovation Lab, Center for Molecular Oncology, Memorial Sloan Kettering Cancer Center.
<br/>
<br/>
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

java -server -Xms4g -Xmx4g -cp Waltz.jar org.mskcc.juber.waltz.countreads.CountReads bam-file coverageThreshold transcripts-file bed-file

#### Generate metrics specific to given regions

java -server -Xms4g -Xmx4g -cp Waltz.jar org.mskcc.juber.waltz.Waltz PileupMetrics mappinngQualityThreshold bam-file reference-fasta bed-file


### 2. Genotyping

java -server -Xms4g -Xmx4g -cp Waltz.jar org.mskcc.juber.waltz.Waltz Genotyping mappinngQualityThreshold bam-file reference-fasta bed-file mutations-maf-file












