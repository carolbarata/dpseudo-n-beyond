# <i>Drosophila pseudoobscura</i> data analysis pipeline

This set of scripts includes raw read trimming, read mapping, mapped read parsing, variant called, variant filtering and cross-referencing across my <i>D. pseudoobscura</i> time series.

This pipeline was run on the St Andrews Bioinformatics cluster (marvin) and any submission scripts are only relevant for marvin users or other linux SGE-managed servers.

S -> submission script

M -> main pipeline script

NB: includes both shell, R and python scripts

## Software requirements

This pipeline requires installing the following software:
```
    Trimmomatic-0.38
    FastQC
    bwa-0.7.17
    novoalign
    java - jdk1.8.0_191
    samtools-1.9
    GenomeAnalysisTK-3.8.1
    sambamba
    bcftools
    bgzip
    freebayes
    python3 or higher (including pandas library)  
```

### Mapping

    S: `full_mapping_newref_submissionscript.sh`
	M: `main_mapping_pipeline.sh`
		`trim_raw_reads_with_trimmomatic.sh #trims raw reads w/ Trimmomatic`
		`get_initial_and_trimmed_read_number.sh #counts no. of reads lost after trimming`
		`get_initial_and_trimmed_read_lengths.sh #fetches read length distns before&after trimming`
		`raw_vs_trimmed_read_length_plot.R #plots read length distns`
		`parsing_mapped_reads_with_samtools.sh #parses sam file, realigns reads around indels, computes mapping stats`
		`mapq_from_samfiles.sh #extracts mapping quality from sam file`
		`coverage_distn_plot.R #plots read coverage distn`
		`raw_vs_trimmed_mapped_read_distn_plot.R #plots raw vs trimmed alignment qual scores`

```
$ ./main_mapping_pipeline.sh -h
    main_mapping_pipeline.sh [-h] [-i1 -i2 -u -dt -df -dm -dn -ds -dg -i -r -o] -- where:
 -h show this help text
 -p full path to pipeline directory
 -i1 full path to forward fastq read file
 -i2 full path to reverse fastq read file
 -u whether trimmed reads should be mapped to the ref (TRUE/FALSE)
 -dt full path to Trimmomatic directory (NA if reads are not to be trimmed)
 -df full path to FastQC directory (NA if reads are not to go through QC)
 -dm full path to bwa mem directory (NA if bwa mem is not to be used)
 -dn full path to NovoAlign directory (NA if NovoAlign is not to be used)
 -dj full path to Java 1.8 directory
 -ds full path to samtools directory
 -dg full path to GATK directory
 -db full path to sambamba directory
 -i whether the reference genome is already indexed or not (TRUE/FALSE)
 -r full path to reference genome
 -o full path to output directory (e.g. /full/path/to/outdir/)
```

### Variant calling

    S: `snp_calling_submissionscript.sh`
	M: `main_snp_calling_pipeline.sh`

```
$ ./main_snp_calling_pipeline.sh
    main_snp_calling_pipeline.sh [-h] [-id -b -s -fb -f -o -r -i] -- where:
 -h show this help text
 -id sample ID/name
 -m name of program used for read mapper (bwa or novoalign)
 -b full path to bcftools directory
 -s full path to samtools directory
 -z full path to bgzip directory
 -fb full path to freebayes directory -f full path to input bam file
 -o full path to output directory
 -r full path to reference fasta file
 -i whether the reference genome is already indexed or not (TRUE or FALSE)
```

### Extract relevant info from VCFs

    S: `data_extracting_submissionscript.sh`
	M: `main_data_extracting_pipeline.sh`

```
$ ./main_data_extracting_pipeline.sh -h
    main_data_extracting_pipeline.sh [-h] [-id -b -s -fb -f -o -r -i] -- where:
 -h show this help text
 -id sample ID/name
 -i full path to all input VCF files
 -b full path to bcftools directory
 -z full path to bgzip directory
 -o full path to output directory
 -c full path to file containing SNP counts per sample
```

### Annotate SNP data files with chromosome names

    S/M: `chromosome_name_annotation_submissionscript.sh`

### Subsetting data into chromosome level files

	S/M: `data_subsetting_submissionscript.sh`

### Filter and intersect bwa/novoalign data

	S: `data_filtering_submissionscript.sh`
	M: `main_sample_snp_filtering_pipeline.py`

```
$ python3 ./main_sample_snp_filtering_pipeline.py -h
usage: main_sample_snp_filtering_pipeline.py [-h] -b B -n N -t T -s S -o O

SNP filtering pipeline

optional arguments:
  -h, --help            show this help message and exit
  -b B, --input_bwa_file B
                        Full path to input tab separated data file for bwa mem
                        mapped reads
  -n N, --input_novoalign_file N
                        Full path to input tab separated data file for
                        novoalign mapped reads
  -t T, --time_point T  Time point in your time series
  -s S, --sample_size S
                        Pooled sequencing sample size
  -o O, --output_file O
                        Full path to output file
```
	
### Intersect SNPs across replicates & time points #future commit

	S: `data_intersecting_submissionscript.sh`
	M: `main_sample_intersecting_pipeline.py`

```
$ python3 ./main_sample_intersecting_pipeline.py -h
usage: main_sample_intersecting_pipeline.py [-h] -f F -r R -g G -c C -o O -q Q

SNP filtering pipeline

optional arguments:
  -h, --help            show this help message and exit
  -f F, --list_of_input_files F
                        Full path to each input tab separated data file
  -r R, --replicates R  Total number of replicate experimental populations
  -g G, --generations G
                        Total number of time points in your time series
  -c C, --chromosome C  Chromosome name
  -o O, --output_directory O
                        Full path to output directory
  -q Q, --qualtity_control_file Q
                        Full path to qualtity control output file
```

