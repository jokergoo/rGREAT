

check_pkg = function(pkg, bioc = FALSE, github = NULL) {
	if(requireNamespace(pkg, quietly = TRUE)) {
		return(NULL)
	} else {

		if(!interactive()) {
			if(bioc) {
				stop_wrap(qq("You need to manually install package '@{pkg}' from Bioconductor."))
			} else {
				stop_wrap(qq("You need to manually install package '@{pkg}' from CRAN."))
			}
		}

		if(bioc) {
			answer = readline(qq("Package '@{pkg}' is required but not installed. Do you want to install it from Bioconductor? [y|n] "))
		} else {
			answer = readline(qq("Package '@{pkg}' is required but not installed. Do you want to install it from CRAN? [y|n] "))
		}

		if(bioc) {
			if(tolower(answer) %in% c("y", "yes")) {
				if(!requireNamespace("BiocManager", quietly = TRUE)) {
					install.packages("BiocManager")
				}
				suppressWarnings(BiocManager::install(pkg, update = FALSE))

				if(!requireNamespace(pkg, quietly = TRUE)) {
					if(is.null(github)) {
						stop_wrap(qq("Cannot find '@{pkg}' from Bioconductor."))
					} else {
						answer = readline(qq("Not on Bioconductor. Install '@{pkg}' from GitHub: '@{github}/@{pkg}'? [y|n] "))
						if(tolower(answer) %in% c("y", "yes")) {
							BiocManager::install(paste0(github, "/", pkg), update = FALSE)
						} else {
							stop_wrap(qq("You need to manually install package '@{pkg}' from CRAN."))
						}
					}
				}
			} else {
				stop_wrap(qq("You need to manually install package '@{pkg}' from Bioconductor."))
			}
		} else {
			if(tolower(answer) %in% c("y", "yes")) {
				suppressWarnings(install.packages(pkg))
				if(!requireNamespace(pkg, quietly = TRUE)) {
					if(is.null(github)) {
						stop_wrap(qq("Cannot find '@{pkg}' from CRAN"))
					} else {
						answer = readline(qq("Not on CRAN. Install '@{pkg}' from GitHub: '@{github}/@{pkg}'? [y|n] "))
						if(tolower(answer) %in% c("y", "yes")) {
							BiocManager::install(paste0(github, "/", pkg), update = FALSE)
						} else {
							stop_wrap(qq("You need to manually install package '@{pkg}' from CRAN."))
						}
					}
				}
			} else {
				stop_wrap(qq("You need to manually install package '@{pkg}' from CRAN."))
			}
		}
	}
}


stop_wrap = function (...) {
    x = paste0(...)
    x = paste(strwrap(x), collapse = "\n")
    stop(x, call. = FALSE)
}

warning_wrap = function (...) {
    x = paste0(...)
    x = paste(strwrap(x), collapse = "\n")
    warning(x, call. = FALSE)
}

message_wrap = function (...) {
    x = paste0(...)
    x = paste(strwrap(x), collapse = "\n")
    message(x)
}


get_table_from_orgdb = function(pattern, orgdb) {
	check_pkg(orgdb)
	all_tb = ls(envir = getNamespace(orgdb))
	i_tb = which(grepl(pattern, all_tb))
	if(length(i_tb)) {
		getFromNamespace(all_tb[i_tb], orgdb)
	} else {
		NULL
	}	
}

# == title
# Generate random regions
#
# == param
# -genome UCSC genome version, e.g. "hg19".
# -nr Number of regions.
# -seqlengths Alternatively, you can also specify a named vector of seqlengths (chromosome lengths).
# -width_fun A function which defines the distribution of region widths.
#
# == details
# The number of regions per chromosome is proportional to the chromsome length.
#
# == example
# gr = randomRegions(genome = "hg19")
# quantile(width(gr))
randomRegions = function (genome = NULL, nr = 1000, seqlengths = NULL,
	width_fun = function(n) runif(n, min = 1000, max = 10000)) {

	if(is.null(genome) && is.null(seqlengths)) {
		stop_wrap("Either `genome` or `seqlengths` should be specified.")
	}

	if(!is.null(genome)) {
	    chr.len = getChromInfoFromUCSC(genome)
	    chromosome = names(chr.len)
	} else {
		chr.len = seqlengths
		chromosome = names(chr.len)
	}

    dl = lapply(seq_along(chr.len), function(i) {
        k = round(nr * chr.len[i]/sum(chr.len))
        breaks = sort(sample(chr.len[i], k))
        
        rr = data.frame(chr = rep(chromosome[i], length(breaks)),
        	       start = breaks,
        	       end = breaks + round(width_fun(length(breaks))))

        rr = rr[rr[, 3] <= chr.len[i], , drop = FALSE]
        rr
    })
    df = NULL
    for (i in seq_along(dl)) {
        df = rbind(df, dl[[i]])
    }

    gr = GRanges(seqnames = df[, 1], ranges = IRanges(df[, 2], df[, 3]))
    seqlevels(gr) = unique(chromosome)
    # if(!is.null(genome)) genome(gr) = genome

    return(gr)
}

# == title
# Generate random regions from a BioMart genome
#
# == param
# -biomart_dataset A BioMart dataset. Values should be in ``BioMartGOGeneSets::supportedOrganisms``.
# -nr Number of regions.
# -... Pass to `randomRegions`.
#
# == details
# The number of regions per chromosome is proportional to the chromsome length.
#
# == example
# if(FALSE) {
#     # Giant panda
#     gr = randomRegionsFromBioMartGenome("amelanoleuca_gene_ensembl")
# }
randomRegionsFromBioMartGenome = function(biomart_dataset, nr = 1000, ...) {
	genes = getGenesFromBioMart(biomart_dataset)
	sl = tapply(end(genes), seqnames(genes), max)
	sl = structure(as.vector(sl), names = names(sl))
	sl = sl[!is.na(sl)]
	randomRegions(nr = nr, seqlengths = sl, ...)
}

guess_gene_id_type = function(id) {
    l = grepl("^\\d+$", id)
    if (sum(l)/length(l) > 0.5) {
        return("ENTREZ")
    }
    l = grepl("^ENS.*G", id)
    if (sum(l)/length(l) > 0.5) {
        return("ENSEMBL")
    }
    l = grepl("^(NC|NG|NM|NR|NP|XM|XR|XP|WP)_\\d+", id)
    if (sum(l)/length(l) > 0.5) {
        return("REFSEQ")
    }
    return(NULL)
}


# == title
# Read gmt gene sets file
#
# == param
# -x The file name of a .gmt file.
# -from Gene ID type in the original gmt file. Value can only take values in 'ENTREZ/SYMBOL/ENSEMBL/REFSEQ'.
# -to Gene ID type that you want to convert to. Value can only take values in 'ENTREZ/SYMBOL/ENSEMBL/REFSEQ'.
# -orgdb The name of an OrgDb database.
#
# == value
# A named list of vectors.
#
# == example
# read_gmt(url("http://dsigdb.tanlab.org/Downloads/D2_LINCS.gmt"))
read_gmt = function(x, from = NULL, to = NULL, orgdb = NULL) {
	ln = readLines(x)
	lt = strsplit(ln, "\t")

	nm = sapply(lt, function(x) x[1])
	v = lapply(lt, function(x) x[-(1:2)])
	names(v) = nm

	if(!is.null(from) && !is.null(to) && !is.null(orgdb)) {

		if(!from %in% c("ENTREZ", "SYMBOL", "ENSEMBL", "REFSEQ")) {
			stop_wrap("`from` can only take values in 'ENTREZ/SYMBOL/ENSEMBL/REFSEQ'.")
		}
		if(!to %in% c("ENTREZ", "SYMBOL", "ENSEMBL", "REFSEQ")) {
			stop_wrap("`to` can only take values in 'ENTREZ/SYMBOL/ENSEMBL/REFSEQ'.")
		}

		if(from == "ENTREZ") {
			map = get_table_from_orgdb(paste0(to, "$"), orgdb)
			map = unlist(as.list(map))
			v = lapply(v, function(x) {
				x2 = map[x]
				x2 = x2[!is.na(x2)]
				unique(x2)
			})
		} else if(to == "ENTREZ") {
			map = get_table_from_orgdb(paste0(from, "2EG$"), orgdb)
			map = unlist(as.list(map))
			v = lapply(v, function(x) {
				x2 = map[x]
				x2 = x2[!is.na(x2)]
				unique(x2)
			})
		} else {
			map1 = get_table_from_orgdb(paste0(from, "2EG$"), orgdb)
			map1 = unlist(as.list(map1))

			map2 = get_table_from_orgdb(paste0(to, "$"), orgdb)
			map2 = unlist(as.list(map2))
			v = lapply(v, function(x) {
				x2 = map2[map1[x]]
				x2 = x2[!is.na(x2)]
				unique(x2)
			})
		}
	}
	v = v[sapply(v, length) > 0]
	v
}


get_url = function(url) {
	nm = basename(url)

	f = paste0(tempdir(), "/", nm)
	if(!file.exists(f)) {
		e = try(suppressWarnings(download.file(url, destfile = f, quiet = TRUE)), silent = TRUE)
		if(inherits(e, "try-error")) {
			stop_wrap(qq("Cannot download file from @{url}."))
		}
	}
	f
}



gs_dataframe2list = function(df) {
    n1 = length(unique(df[, 1]))
    n2 = length(unique(df[, 2]))
    if(n1 < n2) {
        split(df[, 2], df[, 1])
    } else {
        split(df[, 1], df[, 2])
    }
}

# chr_len_db = tapply(BIOC_ANNO_PKGS$txdb, BIOC_ANNO_PKGS$genome_version_in_txdb, function(x) {
# 	x = x[1]
# 	check_pkg(x, bioc = TRUE)
# 	txdb = getFromNamespace(x, ns = x)
# 	seqlengths(txdb)
# })

# chr_len_db = lapply(chr_len_db, function(x) {
# 	x = x[nchar(names(x)) < 10]
# 	if(grepl("$chr", names(x)[1])) {
# 		names(x) = paste0("chr", names(x))
# 	}
# 	x
# })


# manually curated
BIOMART_GENOME = c(vpacos_gene_ensembl = "vicPac1", bbbison_gene_ensembl = "bisBis1",
	dnovemcinctus_gene_ensembl = "dasNov3", gmorhua_gene_ensembl = "gadMor3",
	sbboliviensis_gene_ensembl = "saiBol1", ppaniscus_gene_ensembl = "panPan2",
	ogarnettii_gene_ensembl = "otoGar3", cbrenneri_eg_gene = "caePb2",
	cbriggsae_eg_gene = "cb4", celegans_gene_ensembl = "ce11", celegans_eg_gene = "ce11",
	cjaponica_eg_gene = "caeJap1", cremanei_eg_gene = "caeRem3",
	fcatus_gene_ensembl = "felCat9", ggallus_gene_ensembl = "galGal6",
	ptroglodytes_gene_ensembl = "panTro5", cgchok1gshd_gene_ensembl = "criGriChoV2",
	lchalumnae_gene_ensembl = "latCha1", btaurus_gene_ensembl = "bosTau9",
	mfascicularis_gene_ensembl = "macFas5", clfamiliaris_gene_ensembl = "canFam4",
	ttruncatus_gene_ensembl = "turTru2", dananassae_eg_gene = "droAna3",
	derecta_eg_gene = "droEre2", dgrimshawi_eg_gene = "droGri2",
	dmelanogaster_gene_ensembl = "dm6", dmelanogaster_eg_gene = "dm6",
	dmojavensis_eg_gene = "droMoj3", dpersimilis_eg_gene = "droPer1",
	dpseudoobscura_eg_gene = "dp4", dsechellia_eg_gene = "droSec1",
	dsimulans_eg_gene = "droSim1", dvirilis_eg_gene = "droVir3",
	dwillistoni_eg_gene = "droWil1", dyakuba_eg_gene = "droYak2",
	lafricana_gene_ensembl = "loxAfr3", cmilii_gene_ensembl = "calMil1",
	mpfuro_gene_ensembl = "musFur1", trubripes_gene_ensembl = "fr3",
	nleucogenys_gene_ensembl = "nomLeu3", acchrysaetos_gene_ensembl = "aquChr2",
	rroxellana_gene_ensembl = "rhiRox1", ggorilla_gene_ensembl = "gorGor4",
	acarolinensis_gene_ensembl = "anoCar2", cporcellus_gene_ensembl = "cavPor3",
	eeuropaeus_gene_ensembl = "eriEur1", ecaballus_gene_ensembl = "equCab3",
	hsapiens_gene_ensembl = "hg38", pcapensis_gene_ensembl = "proCap1",
	dordii_gene_ensembl = "dipOrd1", pmarinus_gene_ensembl = "petMar2",
	gfortis_gene_ensembl = "geoFor1", pvampyrus_gene_ensembl = "pteVam1",
	mlucifugus_gene_ensembl = "myoLuc2", mmusculus_gene_ensembl = "mm39",
	hgfemale_gene_ensembl = "hetGla2", oniloticus_gene_ensembl = "oreNil2",
	mdomestica_gene_ensembl = "monDom5", cpbellii_gene_ensembl = "chrPic1",
	sscrofa_gene_ensembl = "susScr11", oprinceps_gene_ensembl = "ochPri2",
	oanatinus_gene_ensembl = "ornAna1", ppacificus_eg_gene = "priPac1",
	ocuniculus_gene_ensembl = "oryCun2", rnorvegicus_gene_ensembl = "rn7",
	scerevisiae_eg_gene = "sacCer3", scerevisiae_gene_ensembl = "sacCer3",
	oaries_gene_ensembl = "oviAri3", saraneus_gene_ensembl = "sorAra1",
	choffmanni_gene_ensembl = "choHof1", itridecemlineatus_gene_ensembl = "speTri2",
	gaculeatus_gene_ensembl = "gasAcu1", csyrichta_gene_ensembl = "tarSyr2",
	sharrisii_gene_ensembl = "sarHar1", tnigroviridis_gene_ensembl = "tetNig2",
	tbelangeri_gene_ensembl = "tupBel1", tcastaneum_eg_gene = "triCas2",
	xtropicalis_gene_ensembl = "xenTro9", mgallopavo_gene_ensembl = "melGal5",
	csabaeus_gene_ensembl = "chlSab2", neugenii_gene_ensembl = "macEug2",
	tguttata_gene_ensembl = "taeGut1", drerio_gene_ensembl = "danRer11"
)


# https://bioinformatics.stackexchange.com/questions/4390/expand-granges-object-different-amounts-upstream-vs-downstream
expandRange = function(x, upstream=2000, downstream=1000) {
  strand_is_minus = strand(x) == "-"
  on_plus = which(!strand_is_minus)
  on_minus = which(strand_is_minus)
  start(x)[on_plus] = start(x)[on_plus] - upstream
  start(x)[on_minus] = start(x)[on_minus] - downstream
  end(x)[on_plus] = end(x)[on_plus] + downstream
  end(x)[on_minus] = end(x)[on_minus] + upstream
  x
}
