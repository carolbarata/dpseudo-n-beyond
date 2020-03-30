# Drosophila pseudoobscura data analysis pipeline

This set of scripts includes raw read trimming, read mapping, mapped read parsing, variant called, variant filtering and cross-referencing across my D. pseudoobscura time series.

S -> submission script
M -> main pipeline script

NB: includes both shell, R and python scripts

## Mapping
    S: full_mapping_newref_submissionscript.sh
	M: main_mapping_pipeline.sh
		trim_raw_reads_with_trimmomatic.sh #trims raw reads w/ Trimmomatic
		get_initial_and_trimmed_read_number.sh #counts no. of reads lost after trimming
		get_initial_and_trimmed_read_lengths.sh #fetches read length distns before&after trimming 
		raw_vs_trimmed_read_length_plot.R #plots read length distns
		parsing_mapped_reads_with_samtools.sh #parses sam file
											  #realigns reads around indels
											  #computes mapping stats
		mapq_from_samfiles.sh #extracts mapping quality from sam file
		coverage_distn_plot.R #plots read coverage distn
		raw_vs_trimmed_mapped_read_distn_plot.R #plots raw vs trimmed alignment qual scores 

## Variant calling
    S: snp_calling_submissionscript.sh
	M: main_snp_calling_pipeline.sh

## Extract relevant info from VCFs
    S: data_extracting_submissionscript.sh
	M: main_data_extracting_pipeline.sh

## Annotate SNP data files with chromosome names
    S/M: chromosome_name_annotation_submissionscript.sh

## Subsetting data into chromosome level files
	S/M: data_subsetting_submissionscript.sh

## Filter and intersect bwa/novoalign data
	S: data_filtering_submissionscript.sh
	M: main_sample_snp_filtering_pipeline.py
	
## Intersect SNPs across replicates & time points #future commit
	S: data_intersecting_submissionscript.sh
	M: main_sample_intersecting_pipeline.py

