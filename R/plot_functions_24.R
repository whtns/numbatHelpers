# Plot Functions (123)
#' Create a plot visualization
#'
#' @param seu_path File path
#' @param mymarkers Parameter for mymarkers
#' @param plot_type Parameter for plot type
#' @param group_by Parameter for group by
#' @return ggplot2 plot object
#' @export
plot_markers_in_sample <- function(seu_path, mymarkers, plot_type = FeaturePlot, group_by = group_by) {
  #
  sample_id <- str_extract(seu_path, "SR[RX][0-9]+")

  numbat_dir <- "numbat_sridhar"

  dir_create(glue("results/{numbat_dir}"))
  dir_create(glue("results/{numbat_dir}/{sample_id}"))

  seu <- readRDS(seu_path)

  mymarkers <- mymarkers[mymarkers %in% rownames(seu)]

  if (identical(plot_type, VlnPlot)) {
    feature_plots_first <- plot_type(seu, features = mymarkers, group.by = group_by, combine = FALSE, pt.size = 0)

    names(feature_plots_first) <-
      map_chr(feature_plots_first, ~ .x[["labels"]][["title"]])

    max_ys <- map(feature_plots_first, ~ layer_scales(.x)$y$get_limits()) %>%
      map(2) %>%
      identity()

    feature_plots <- map2(feature_plots_first, max_ys, ~ {
      .x +
        # expand_limits(y = c(0, .y*2.5)) +
        stat_compare_means(comparisons = list(c(1, 2)), method = "t.test", label.y = .y * 0.9) +
        # stat_compare_means() +
        geom_boxplot(width = 0.2) +
        theme(
          legend.position = "none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank()
        ) +
        # stat_compare_means(method = "anova", label.y= 0.4) +
        NULL
    })
  } else if (identical(plot_type, FeaturePlot)) {
    feature_plots <- plot_type(seu, features = mymarkers, combine = FALSE) %>%
      set_names(mymarkers)
  }


  feature_plots <- map(feature_plots, ~ (.x + labs(title = sample_id)))

  return(feature_plots)
}

