# Numbat Functions (14)

#' Perform retrieve segmentation operation
#'
#' @param rds_file File path
#' @return Function result
#' @export
retrieve_segmentation <- function(rds_file){
		joint_post <- readRDS(rds_file)[["joint_post"]]
		return(joint_post)
	}

#' Calculate or compute values
#'
#' @param df Input data frame or dataset
#' @return Calculated values
#' @export
calculate_arm_percent <- function(df) {
	
	df <- df |>
		dplyr::rename(seqnames = chrom, start = chromstart, end = chromend) |> 
		dplyr::select(-width) |> 
		# dplyr::filter(cnv_state == "del") |> 
		as_granges()
	
	# arms_df <- 
	"http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz" |> 
		read_tsv(col_names = c("chrom", "chromStart", "chromEnd", "name", "gieStain")
		) |>
		mutate(arm = substring(name, 1, 1)) |>
		group_by(chrom, arm) |>
		summarise(
			start = min(chromStart),
			end = max(chromEnd),
			length = end - start
		) |>
		dplyr::mutate(seqnames = chrom) |>
		dplyr::mutate(armstart = start, armend = end) |> 
		as_granges() |>
		join_overlap_intersect(df) |>
		dplyr::mutate(percent_start = (start-armstart)/(length)) |> 
		dplyr::mutate(percent_end = (end-armstart)/(length)) |> 
		as_tibble() |> 
		dplyr::select(chr = seqnames, arm, start = percent_start, end = percent_end, cnv_state) |> 
		identity()
	
}

#' Calculate or compute values
#'
#' @return Calculated values
#' @export
calculate_6p_percent_affected <- function(){
	excel_file <- "results/table_s12.xlsx"
	
	test1 <- excel_sheets(excel_file) |> 
		set_names() |> 
		map(~readxl::read_excel(excel_file, sheet = .x)) |> 
		map(calculate_arm_percent) |>
		identity()
	
	# test3 <- 
	test1[c("SRX10264524", "SRX10264525", "SRX14116944", "SRX22868105")] |> 
		dplyr::bind_rows(.id = "tumor") |> 
		dplyr::filter(chr == "chr6", arm == "p") |> 
		dplyr::mutate(start = percent(start, accuracy = 1)) |> 
		dplyr::mutate(end = percent(end, accuracy = 1)) |> 
		down_csv()
	
}

#' Perform binom test clone percent operation
#'
#' @param mydf Parameter for mydf
#' @return Function result
#' @export
binom_test_clone_percent <- function(mydf) {
	
	whole_tumor_prop = sum(mydf[[3]])/sum(mydf[[2]], mydf[[3]])
	
	success_var = colnames(mydf)[[3]]
	
	# browser()
	
	test_df <-
		mydf |> 
		dplyr::filter(.data[[success_var]] > 0) |> 
		mutate(trials = rowSums(across(where(is.numeric)), na.rm = TRUE)) |> 
		dplyr::rowwise() |> 
		dplyr::mutate(test_result = list(binom.test(.data[[success_var]], trials, p = whole_tumor_prop, alternative = "two.sided"))) |>
		dplyr::mutate(test_result = broom::tidy(test_result)) |>
		tidyr::unnest(test_result) |>
		dplyr::group_by(clusters) |>
		dplyr::select(clusters, p.value) |> 
		dplyr::right_join(mydf, by = "clusters") |> 
		dplyr::filter(.data[[success_var]] > 0) |> 
		adorn_percentages(denominator = "col",, any_of(colnames(mydf)[2:3])) |>
		adorn_pct_formatting(digits = 2,,, any_of(colnames(mydf)[2:3])) |> 
		identity()
}

#' Perform tabulate clone percent operation
#'
#' @param seu_path File path
#' @return Function result
#' @export
tabulate_clone_percent <- function(seu_path = "output/seurat/SRX11133594_filtered_seu.rds") {

	sample_id <- str_extract(seu_path, "SR[RX][0-9]+")
		
	seu <- readRDS(seu_path)
	
	# seu <- seu[,!str_detect(seu$clusters, "other")]
	
	input_sample <- 
		janitor::tabyl(seu@meta.data, clusters, scna)
	
	clone_percent_df <- 
		input_sample |> 
		binom_test_clone_percent() |> 
		dplyr::mutate(p.adjust = p.adjust(p.value, method = "BH")) |>
		dplyr::mutate(
			signif =
				symnum(p.adjust,
							 corr = FALSE,
							 cutpoints = c(0, .001, .01, .05, .1, 1),
							 symbols = c("***", "**", "*", ".", " ")
				)
		) |>
		identity()
	
}

