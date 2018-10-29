#! /bin/bash

# $1 - module to run: PileupMetrics, Genotyping
# $2 - minimum mapping quality
# $3 - bam file
# $4 - reference fasta file
# $5 - interval bed file
# $6 - module argument
	# empty for Metrics module
	# maf file with mutations to be profiled
	# comma-separated list of signatures to look for for SignatureFinding module

java=/opt/common/CentOS_6/java/jdk1.8.0_31/bin/java

export TMPDIR=/ifs/work/scratch

$java -server -Xms4g -Xmx4g -cp ~/software/Waltz.jar org.mskcc.juber.waltz.Waltz $@
