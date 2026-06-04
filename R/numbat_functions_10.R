# Numbat Functions (10)

#' Perform tabulate clone percent integrated operation
#'
#' @param seu_path File path
#' @return Function result
#' @export
tabulate_clone_percent_integrated <- function(seu_path = "output/seurat/integrated_1q/integrated_seu_1q_complete.rds") {
	
	seu_integrated <- readRDS(seu_path)
	
	seu_integrated <- seu_integrated[,!str_detect(seu_integrated$clusters, "other")]
	
	integrated_and_split_samples <- 
		janitor::tabyl(seu_integrated@meta.data, clusters, scna, batch)
	
	integrated_sample <- 
		list("integrated" = janitor::tabyl(seu_integrated@meta.data, clusters, scna))
	
	table_input <- c(integrated_sample, integrated_and_split_samples)
	
	clone_percent_df <- 
		table_input |> 
		map(binom_test_clone_percent) |> 
		dplyr::bind_rows(.id = "sample_id") |>
		dplyr::mutate(p.adjust = p.adjust(p.value)) |>
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

#' Perform make table s05 operation
#'
#' @param table_path File path
#' @return Function result
#' @export
make_table_s05 <- function(table_path = "results/table_s05.xlsx") {
	
	clone_percent_1q = tabulate_clone_percent_integrated("output/seurat/integrated_1q/integrated_seu_1q_complete.rds")
	clone_percent_16q = tabulate_clone_percent_integrated("output/seurat/integrated_16q/integrated_seu_16q_complete.rds")
	clone_percent_2p = tabulate_clone_percent_integrated("output/seurat/integrated_2p/seurat_2p_integrated_duo.rds")
	
	list("1q+" = clone_percent_1q, 
			 "16q-" = clone_percent_16q,
			 "2p+" = clone_percent_2p) |> 
		writexl::write_xlsx(table_path)
	
	return(table_path)
}

#' Perform get chr arms operation
#'
#' @param build Character string (default: "hg19")
#' @return Function result
#' @export
get_chr_arms <- function(build = "hg19"){
	glue::glue("http://hgdownload.cse.ucsc.edu/goldenpath/{build}/database/cytoBand.txt.gz") |> 
		read_tsv(col_names = c("chrom", "chromStart", "chromEnd", "name", "gieStain")
		) |>
		mutate(arm = substring(name, 1, 1)) |>
		group_by(chrom, arm) |>
		summarise(
			start = min(chromStart), 
			end = max(chromEnd),
			length = end - start
		) |> 
		dplyr::select(chrom, arm, length) |> 
		dplyr::filter(chrom %in% paste0("chr", 1:22)) |> 
		dplyr::mutate(chrom = as.numeric(str_remove(chrom, "chr"))) |> 
		dplyr::arrange(chrom) |> 
		tidyr::unite("chr", chrom:arm, sep = "") |> 
		tibble::deframe()
	
}

#' Perform get arms ranges operation
#'
#' @return Function result
#' @export
get_arms_ranges <- function() {
	arms_df <- readr::read_tsv("http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz",
														 col_names = c("chrom", "chromStart", "chromEnd", "name", "gieStain")
	) |>
		mutate(arm = substring(name, 1, 1)) |>
		group_by(chrom, arm) |>
		summarise(
			start = min(chromStart),
			end = max(chromEnd),
			length = end - start
		) |>
		dplyr::mutate(chrom = stringr::str_remove(chrom, "chr")) |>
		dplyr::mutate(seqnames = chrom) |>
		dplyr::filter(seqnames %in% c(1:22, "X", "Y")) |> 
		dplyr::mutate(chrom_arm = str_c(chrom, arm)) |> 
		plyranges::as_granges()
}

#' Perform annotables to grange operation
#'
#' @return Function result
#' @export
annotables_to_grange <- function(){
	annotables::grch38 |> 
		dplyr::mutate(seqnames = chr) |> 
		plyranges::as_granges()
}

