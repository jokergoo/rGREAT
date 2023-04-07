
# == title
# Make volcano plot
#
# == param
# -object A `GreatJob-class` object returned by `submitGreatJob`.
# -ontology A single ontology names. Valid values are in `availableOntologies`. 
# -min_region_hits Minimal number of input regions overlapping to the geneset associated regions.
# -x_values Which values for the x-axis.
# -y_values Which values for the y-axis.
# -main Title of the plot.
#
# == details
# Since the enrichment is an over-representation test, it is only the half volcano.
setMethod(f = "plotVolcano",
	signature = "GreatJob",
	definition = function(object, ontology, min_region_hits = 5, 
	x_values = c("fold_enrichment", "z-score"),
	y_values = c("p_value", "p_adjust"),
	main = NULL) {

	tb = getEnrichmentTables(object, ontology = ontology, verbose = FALSE)[[1]]

	x_values = match.arg(x_values)[1]
	y_values = match.arg(y_values)[1]

	n = param(object, "n_region")
	if(is.null(main)) {
		main = paste0("Volcano plot for ", n, " regions, ontology: ", ontology)
	}
	plot_volcano(
		observed_region_hits = tb$Binom_Observed_Region_Hits,
		genome_fraction = tb$Binom_Genome_Fraction,
		fold_enrichment = tb$Binom_Fold_Enrichment,
		z_score = (tb$Binom_Observed_Region_Hits - n*tb$Binom_Genome_Fraction)/sqrt(n*tb$Binom_Genome_Fraction*(1-tb$Binom_Genome_Fraction)),
		p_value = tb$Binom_Raw_PValue,
		p_adjust = tb$Binom_Adjp_BH,
		min_region_hits = min_region_hits,
		x_values = x_values,
		y_values = y_values,
		main = main
	)

})

# == title
# Make volcano plot
#
# == param
# -object A `GreatObject-class` object returned by `great`.
# -min_region_hits Minimal number of input regions overlapping to the geneset associated regions.
# -x_values Which values for the x-axis.
# -y_values Which values for the y-axis.
# -main Title of the plot.
#
# == details
# Since the enrichment is an over-representation test, it is only the half volcano.
setMethod(f = "plotVolcano",
	signature = "GreatObject",
	definition = function(object, min_region_hits = 5, 
	x_values = c("fold_enrichment", "z-score"),
	y_values = c("p_value", "p_adjust"),
	main = NULL) {

	x_values = match.arg(x_values)[1]
	y_values = match.arg(y_values)[1]

	tb = object@table
	n = length(object@gr)
	if(is.null(main)) {
		main = paste0("Volcano plot for ", n, " regions")
	}
	plot_volcano(
		observed_region_hits = tb$observed_region_hits,
		genome_fraction = tb$genome_fraction,
		fold_enrichment = tb$fold_enrichment,
		z_score = (tb$observed_region_hits - n*tb$genome_fraction)/sqrt(n*tb$genome_fraction*(1-tb$genome_fraction)),
		p_value = tb$p_value,
		p_adjust = tb$p_adjust,
		min_region_hits = min_region_hits,
		x_values = x_values,
		y_values = y_values,
		main = main
	)

})

plot_volcano = function(
	observed_region_hits,
	genome_fraction,
	fold_enrichment,
	z_score,
	p_value,
	p_adjust,
	min_region_hits = 5, 
	x_values = "fold_enrichment",
	y_values = "p_value",
	main = NULL) {

	col_fun = colorRamp2(seq(min_region_hits, min_region_hits + min(50, max(observed_region_hits)), length = 11), rev(brewer.pal(11, "Spectral")))
	size_fun = function(x) x^0.3*2 + 0.1

	if(x_values == "fold_enrichment") {
		x = log2(fold_enrichment)
		xlab = "log2 fold enrichment: log2(obs/exp)"
	} else {
		x = z_score
		xlab = "z-score: (obs-exp)/sd"
	}

	
	if(y_values == "p_value") {
		y = -log10(p_value)
		ylab = "-log10(p_value)"
	} else {
		y = -log10(p_adjust)
		ylab = "-log10(p_adjust)"
	}
	y[is.infinite(y)] = max(y[is.finite(y)])

	l = observed_region_hits >= min_region_hits
	if(!any(l)) {
		plot(NULL, xlim = c(0, 1), ylim = c(0, 1), ann = FALSE, axes = FALSE)
		text(0.5, 0.5, qq("No term left with min_region_hits >= @{min_region_hits}"))
		return(invisible(NULL))
	}

	plot(x[l], y[l], pch = 16,
		col = col_fun(observed_region_hits[l]),
		cex = size_fun(genome_fraction[l]),
		xlab = xlab, ylab = ylab)
	if(is.null(main)) {
		title("Volcano plot")
	} else {
		title(main)
	}
	abline(h = -log10(0.05), col = "grey", lty = 2)
	text(par("usr")[2], -log10(0.05), qq("@{y_values} = 0.05"), adj = c(1.05, 1.5), cex = 0.8)
	if(x_values == "fold_enrichment") abline(v = 1, col = "grey", lty = 2)

	max_region_hits = min_region_hits + min(50, max(observed_region_hits))
	breaks = seq(min_region_hits, max_region_hits, length = 5)
	breaks = round(breaks)
	size = legend("topleft", pch = 16, cex = 0.8, legend = breaks, col = col_fun(breaks), 
		title = "Region hits", bty = "n")

	rg = range(genome_fraction[l])
	size_rg = size_fun(rg)
	breaks = ((seq(size_rg[1], size_rg[2], length = 5) - 0.1)/2)^(1/0.3)

	size = legend(x = size$rect$left, y = size$rect$top - size$rect$h*1.05,
		pch = 16, pt.cex = size_fun(breaks), cex = 0.8,
		legend = paste0(sprintf("%.2f", breaks*100), "%"), 
		col = "black", 
		title = "Genome fraction", bty = "n")

	text(size$rect$left, size$rect$top - size$rect$h, qq("Region hits >= @{min_region_hits}"), adj = c(-0.05, 2), cex = 0.8)

}

