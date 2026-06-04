# Plot Functions (109)

#' Perform reorder within operation
#'
#' @param x Parameter for x
#' @param by Parameter for by
#' @param within Parameter for within
#' @param fun Parameter for fun
#' @param sep Character string (default: "___")
#' @param ... Additional arguments passed to other functions
#' @return Function result
#' @export
# Performance optimizations applied:
# - multiple_joins: Combine multiple joins into single join operation where possible

reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
	new_x <- paste(x, within, sep = sep)
	stats::reorder(new_x, by, FUN = fun)
}

#' Perform scale x reordered operation
#'
#' @param ... Additional arguments passed to other functions
#' @param sep Character string (default: "___")
#' @return Function result
#' @export
scale_x_reordered <- function(..., sep = "___") {
	reg <- paste0(sep, ".+$")
	ggplot2::scale_x_discrete(labels = function(x) gsub(reg, "", x), ...)
}

#' Perform scale y reordered operation
#'
#' @param ... Additional arguments passed to other functions
#' @param sep Character string (default: "___")
#' @return Function result
#' @export
scale_y_reordered <- function(..., sep = "___") {
	reg <- paste0(sep, ".+$")
	ggplot2::scale_y_discrete(labels = function(x) gsub(reg, "", x), ...)
}

#' Perform make fig s16 table s04 operation
#'
#' @param plot_path File path
#' @param table_path File path
#' @param ... Additional arguments passed to other functions
#' @return ggplot2 plot object
#' @export
rb_scna_frequency_in_tcga_by_cancer_type <- function(plot_path = "results/fig_s02.pdf", table_path = "results/table_s04.csv", ...) {
	tcga_abbreviations <- 
		tibble::tribble(
			~study_abbreviation,                                                        ~study_name,
			"LAML",                                           "Acute Myeloid Leukemia",
			"ACC",                                         "Adrenocortical carcinoma",
			"BLCA",                                     "Bladder Urothelial Carcinoma",
			"LGG",                                         "Brain Lower Grade Glioma",
			"BRCA",                                        "Breast invasive carcinoma",
			"CESC", "Cervical squamous cell carcinoma and endocervical adenocarcinoma",
			"CHOL",                                               "Cholangiocarcinoma",
			"LCML",                                     "Chronic Myelogenous Leukemia",
			"COAD",                                             "Colon adenocarcinoma",
			"CNTL",                                                         "Controls",
			"ESCA",                                             "Esophageal carcinoma",
			"FPPP",                                              "FFPE Pilot Phase II",
			"GBM",                                          "Glioblastoma multiforme",
			"HNSC",                            "Head and Neck squamous cell carcinoma",
			"KICH",                                               "Kidney Chromophobe",
			"KIRC",                                "Kidney renal clear cell carcinoma",
			"KIRP",                            "Kidney renal papillary cell carcinoma",
			"LIHC",                                   "Liver hepatocellular carcinoma",
			"LUAD",                                              "Lung adenocarcinoma",
			"LUSC",                                     "Lung squamous cell carcinoma",
			"DLBC",                  "Lymphoid Neoplasm Diffuse Large B-cell Lymphoma",
			"MESO",                                                     "Mesothelioma",
			"MISC",                                                    "Miscellaneous",
			"OV",                                "Ovarian serous cystadenocarcinoma",
			"PAAD",                                        "Pancreatic adenocarcinoma",
			"PCPG",                               "Pheochromocytoma and Paraganglioma",
			"PRAD",                                          "Prostate adenocarcinoma",
			"READ",                                            "Rectum adenocarcinoma",
			"SARC",                                                          "Sarcoma",
			"SKCM",                                          "Skin Cutaneous Melanoma",
			"STAD",                                           "Stomach adenocarcinoma",
			"TGCT",                                      "Testicular Germ Cell Tumors",
			"THYM",                                                          "Thymoma",
			"THCA",                                                "Thyroid carcinoma",
			"UCS",                                           "Uterine Carcinosarcoma",
			"UCEC",                             "Uterine Corpus Endometrial Carcinoma",
			"UVM",                                                   "Uveal Melanoma"
		) |> 
		janitor::clean_names()

# asdfff
	tcga_cohorts <- 
		tcga_gistic_available() |> 
		dplyr::filter(CopyNumberLevel == "all") |> 
		dplyr::select(FullName, Cohort) |> 
		dplyr::mutate(FullName = glue("{FullName} ({Cohort})")) |> 
		tibble::deframe() |>
		identity()

	possibly_check_tcga_peaks <- possibly(check_tcga_peaks)

	sig_peaks <- 
		tcga_cohorts |> 
		map(possibly_check_tcga_peaks) |>
    dplyr::bind_rows(.id = "cohort") |>
    tidyr::pivot_longer(cols = c("X1q", "X2p", "X6p", "X16q"), names_to = "arm", values_to = "is_peak") |>
    dplyr::mutate(arm = dplyr::case_when(
      arm == "X1q" ~ "1q+",
      arm == "X2p" ~ "2p+",
      arm == "X6p" ~ "6p+",
      arm == "X16q" ~ "16q-",
    )) |>
    identity()

	taylor_freq0 <- read_xlsx("data/taylor_et_al_2018_genomic_and_functional_approaches_to_understanding_cancer_aneuploidy/1-s2.0-S1535610818301119-mmc2.xlsx", skip = 1) |> 
		clean_names() |> 
		identity()
	
	taylor_freq1 <- 
		taylor_freq0 |> 
		tidyr::pivot_longer(starts_with("x"), names_to = "arm", values_to = "change") |> 
		dplyr::mutate(arm = str_remove(arm, "x")) |> 
		dplyr::filter(arm %in% c("1q", "2p", "6p", "16q")) |> 
		dplyr::mutate(rb_scna = case_when(
			(arm == "1q" & change == 1) ~ 1,
			(arm == "2p" & change == 1) ~ 1,
			(arm == "6p" & change == 1) ~ 1,
			(arm == "16q" & change == -1) ~ -1,
			.default = 0
		)) |> 
		dplyr::mutate(arm = case_when(
			arm == "1q" ~ "1q+",
			arm == "2p" ~ "2p+",
			arm == "6p" ~ "6p+",
			arm == "16q" ~ "16q-",
		)) |> 
		# dplyr::summarize(percent_affected = sum(rb_scna)) |>
		dplyr::mutate(arm = factor(arm, levels = c("1q+", "2p+", "6p+", "16q-"))) |> 
		dplyr::left_join(tcga_abbreviations, by = c("type" = "study_abbreviation")) |> 
		dplyr::group_by(type, study_name, arm) |> 
		dplyr::summarize(percent_affected = abs(sum(rb_scna)/dplyr::n())) |>
		dplyr::arrange(arm, desc(percent_affected)) |> 
    dplyr::left_join(sig_peaks, by = c("type" = "abbreviation", "arm" = "arm")) |>
		identity()
	
	taylor_freq1 |> 
  mutate(type = ifelse(is_peak, paste0(type, "*"), type)) |> 
  ggplot(aes(x = percent_affected, y = reorder_within(type, percent_affected, arm), fill = type)) + 
		geom_col() + 
		scale_y_reordered() +
		facet_wrap(~arm, scales = "free_y") + 
		theme_minimal() + 
		theme(axis.text.x = element_text(angle = 45, hjust = 1),
					axis.title.y = element_blank())
	
	ggsave(plot_path, ...)
	
	write_csv(taylor_freq1, table_path)
	
	return(list(plot_path, table_path))
}

#' Perform differential expression analysis
#'
#' @param seu Seurat object
#' @param enrichments Parameter for enrichments
#' @param diffexes Parameter for diffexes
#' @param myset Character string (default: "HALLMARK_E2F_TARGETS")
#' @param mycomparison Character string (default: "g1_6 v. g1_0")
#' @return Differential expression results
#' @export
heatmap_from_cluster_diffex <- function(seu, enrichments, diffexes, myset = "HALLMARK_E2F_TARGETS", mycomparison ="g1_6 v. g1_0") {
	genes_in_set <- 
		enrichments[[mycomparison]] |> 
		# dplyr::filter(ID=="HALLMARK_G2M_CHECKPOINT") |> 
		DOSE::setReadable(OrgDb = org.Hs.eg.db, keyType = "ENTREZID") |> 
		slot("result") |> 
		dplyr::filter(ID==myset) |> 
		dplyr::mutate(geneID = str_split(geneID, "/")) |> 
		dplyr::pull(geneID)
	
	comparisons = str_split_1(mycomparison, " v. ")
	
	genes_in_1q <- find_genes_by_arm(genes_in_set[[1]]) |> 
		dplyr::filter(seqnames == "01" & arm == "q") |> 
		dplyr::pull(symbol)
	
	seu[,seu$clusters %in% comparisons] |> 
		seu_complex_heatmap(features = genes_in_set[[1]], group.by = "clusters")
	
	# diffexes[[mycomparison]] |> 
	# 	tibble::rownames_to_column("symbol") |> 
	# 	# dplyr::filter(symbol %in% genes_in_set[[1]]) |> 
	# 	dplyr::filter(symbol %in% genes_in_1q) |> 
	# 	identity()
}

