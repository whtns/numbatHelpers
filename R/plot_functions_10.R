# Plot Functions (102)

#' Create a plot visualization
#'
#' @param seu_list Parameter for seu list
#' @param scna_of_interest Character string (default: "1q")
#' @param ... Additional arguments passed to other functions
#' @return ggplot2 plot object
#' @export
clone_cc_plots_by_scna <- function(seu_list, scna_of_interest = "1q", ...) {
	test0 <- map(seu_list, plot_clone_cc_plots, scna_of_interest = scna_of_interest, ...)
	print(test0)
	
	myplot <- wrap_plots(test0) + 
		plot_layout(ncol = 1)
	
	plot_path <- ggsave(glue("results/{scna_of_interest}_clone_distribution.pdf"), myplot, h = length(test0)*16/5, w = length(test0)*8/5)
	
	return(plot_path)
}

#' Perform add image slide operation
#'
#' @param ppt Parameter for ppt
#' @param image1 Parameter for image1
#' @param image2 Parameter for image2
#' @param h1 Character string (default: "heading 1")
#' @param h2 Character string (default: "heading 2")
#' @return Function result
#' @export
add_image_slide <- function(ppt, image1, image2, h1 = "heading 1", h2 = "heading 2"){
	# 
	ppt <- add_slide(ppt, layout = "Blank")
	
	ppt <- ph_with(
		x = ppt,
		value = h1,
		location = ph_location(left = 1, top = 0)
	)
	
	ppt <- ph_with(
		x = ppt,
		value = external_img(src = image1),
		location = ph_location(left = 1, top = 3)
	)
	
	ppt <- ph_with(
		x = ppt,
		value = h2,
		location = ph_location(left = 7, top = 0)
	)
	
	ppt <- ph_with(
		x = ppt,
		value = external_img(src = image2),
		location = ph_location(left = 7, top = 3)
	)
	
	return(ppt)
}

#' Create a plot visualization
#'
#' @param plots1 Parameter for plots1
#' @param plots2 Parameter for plots2
#' @param file File path
#' @param ... Additional arguments passed to other functions
#' @return ggplot2 plot object
#' @export
paired_plots_to_pptx <- function(plots1, plots2, file, ...){
	# 
	
	ppt <- read_pptx("~/blank_template.pptx")
	
	
#' Create a plot visualization
#'
#' @param ppt Parameter for ppt
#' @param plot1 Parameter for plot1
#' @param plot2 Parameter for plot2
#' @param h1 Character string (default: "heading 1")
#' @param h2 Character string (default: "heading 2")
#' @return ggplot2 plot object
#' @export
add_plots_slide <- function(ppt, plot1, plot2, h1 = "heading 1", h2 = "heading 2"){
		# 
		ppt <- add_slide(ppt, layout = "Blank")
		
		ppt <- ph_with(
			x = ppt,
			value = h1,
			location = ph_location(left = 1, top = 0)
		)
		
		plot_path <- ggsave("plot_image.png", plot = plot1, width = 6, height = 4)
		
		ppt <- ph_with(
			x = ppt,
			value = external_img(src = plot_path),
			location = ph_location(left = 1, top = 3)
		)
		
		ppt <- ph_with(
			x = ppt,
			value = h2,
			location = ph_location(left = 7, top = 0)
		)
		
		plot_path <- ggsave("plot_image.png", plot = plot1, width = 6, height = 4)
		
		ppt <- ph_with(
			x = ppt,
			value = external_img(src = plot_path),
			location = ph_location(left = 7, top = 3)
		)
		
		return(ppt)
	}
	
	for(index in seq_along(plots1)){
		ppt <- add_plots_slide(ppt, plots1[[index]], plots2[[index]], ...)
	}
	
	print(ppt, target = file)
	
	return(file)
	
}
#' Perform clustering analysis
#'
#' @param x Character string (default: "0")
#' @param y Character string (default: "1")
#' @param df_list Parameter for df list
#' @return Calculated values
#' @export
calculate_avg_cluster_distance <- function(x = "0", y = "1", df_list) {
	# Calculate pairwise distances
	
	cluster_x = df_list[[x]]
	cluster_y = df_list[[y]]
	
	distances <- dist(rbind(cluster_x, cluster_y))
	
	# Convert distance object to matrix
	dist_matrix <- as.matrix(distances)
	
	# Extract only the distances between X and Y points
	cross_distances <- dist_matrix[1:nrow(cluster_x), 
																 (nrow(cluster_x) + 1):ncol(dist_matrix)]
	
	# Calculate average distance
	avg_distance <- mean(cross_distances)
	
	return(avg_distance)
}

