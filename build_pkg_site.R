writeLines("destination: docs", "_pkgdown.yml")
if(file.exists(".Rbuildignore")) {
	ln = readLines(".Rbuildignore")
	if(!any(ln == "^_pkgdown\\.yml$")) {
		ln = c(ln, "^_pkgdown\\.yml$")
	}
	if(!any(ln == "^docs$")) {
		ln = c(ln, "^docs$")
	}
	if(!any(ln == "^pkgdown$")) {
		ln = c(ln, "^pkgdown$")
	}
	writeLines(ln, ".Rbuildignore")
} else {
	writeLines("
^_pkgdown\\.yml$
^docs$
^pkgdown$", ".Rbuildignore")
}

pkgdown::build_site()
