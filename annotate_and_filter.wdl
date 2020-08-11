version 1.0

import "https://raw.githubusercontent.com/alexanderhsieh/EM-mosaic-pipeline/master/tasks/tasks_annotation.wdl" as annotation

workflow annotate_and_filter {

	input {
		File input_variants

		## annotation
		String ref_ver
		File cache_dir
		File rr_map 
		File rr_seg 
		File rr_lcr 

		## apply filters
		String output_prefix
		String CAF_outprefix
		Int cohort_size
		Int expected_dnsnvs 
		Int case_cutoff 

	}

	##########################################################################
	## Annotate raw de novo variants
	##########################################################################
	
	# convert txt to vcf
	call annotation.txt_to_vcf {
		input:
			variants = input_variants
	}

	# generate VEP annotations
	call annotation.run_vep {
		input:
			ref = ref_ver,
			vcf = txt_to_vcf.out,
			vcf_idx = txt_to_vcf.idx,
			cache_dir = cache_dir
	}

	# Parse and append VEP columns to original vcf file
	call annotation.add_vep_cols {
		input:
			original_variants = input_variants,
			vep_vcf = run_vep.vep_out
	}

	# flag GATK RankSum values (similar to samtools PV4)
	call annotation.flag_PV4 {
		input:
			infile = add_vep_cols.out
	}

	# flag SB (strand bias) 
	call annotation.flag_SB {
		input:
			infile = flag_PV4.out
	}

	# flag FDR (FDR-based min altdp) 
	call annotation.flag_FDR {
		input:
			infile = flag_SB.out
	}

	# flag RR (repeat region) 
	call annotation.flag_RR {
		input:
			infile = flag_FDR.out,
			map = rr_map,
			seg = rr_seg,
			lcr = rr_lcr
	}

	# flag VC (variant cluster) 
	call annotation.flag_VC {
		input:
			infile = flag_RR.out
	}

	#########################################################
	## Steps from Apply Filters
	#########################################################
	# estimate cohort AF from raw de novos file
	call annotation.estimate_cohort_AF {
		input:
			infile = flag_VC.out,
			cohort_size = cohort_size
	}

	# run CAF (cohort allele frequency)
	call annotation.flag_CAF {
		input:
			infile = flag_VC.out,
			caf_file = estimate_cohort_AF.out
	}



	#########################################################
	## parse filter flags, summarize filtering, output variants passing all filters
	#########################################################
	#run update_filter_column script to combine filter flags into single column
	call annotation.update_filt_col as update1 {
		input:
			infile = flag_CAF.out
	}

	#run outlier filter
	## NOTE: REQUIRES CLEAN DE NOVOS TO ACCURATELY IDENTIFY OUTLIERS
	call annotation.flag_outlier {
		input:
			infile = update1.outfile,
			cohort_size = cohort_size,
			exp = expected_dnsnvs,
			cutoff = case_cutoff
	}

	# 2nd run to update filter column with outlier flag information
	call annotation.update_filt_col as update2 {
		input:
			infile = flag_outlier.out
	}

	#run script to summarize counts of variants flagged by each filter
	call annotation.summarize_counts {
		input:
			infile = update2.outfile
	}

	#run script to write out variants passing all filters, to be used as input to EM-mosaic
	call annotation.print_pass_vars {
		input:
			infile = update2.outfile
	}

	output {
		File denovos_all = update2.outfile # full denovos table with all annotations
		File denovos_PASS = print_pass_vars.outfile # only denovos passing all filters, to be passed to EM-mosaic scoring step
		File filter_summary = summarize_counts.outfile # summary table detailing how many variants failed each filter step
	}


}