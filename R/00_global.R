DEFAULT_VERSION = "4.0.4"


BASE_URL_LIST = c("default" = "http://great.stanford.edu/public/cgi-bin",
	         "4.0.4" = "http://great.stanford.edu/public-4.0.4/cgi-bin",
	         "3.0.0" = "http://great.stanford.edu/public-3.0.0/cgi-bin",
	         "2.0.2" = "http://great.stanford.edu/public-2.0.2/cgi-bin")
BASE_URL_LIST["3.0"] = BASE_URL_LIST["3"] = BASE_URL_LIST["3.0.0"]
BASE_URL_LIST["2.0"] = BASE_URL_LIST["2"] = BASE_URL_LIST["2.0.2"]
BASE_URL_LIST["4.0"] = BASE_URL_LIST["4"] = BASE_URL_LIST["4.0.4"]

GENOME = list("4.0.4" = c("hg38", "hg19", "mm9", "mm10"),
	           "3.0.0" = c("hg19", "mm9", "mm10", "danRer7"),
	           "2.0.2" = c("hg19", "hg18", "mm9", "danRer7"))
GENOME$default = GENOME[[DEFAULT_VERSION]]

GENOME["4.0"] = GENOME["4"] = GENOME["4.0.4"]
GENOME["3.0"] = GENOME["3"] = GENOME["3.0.0"]
GENOME["2.0"] = GENOME["2"] = GENOME["2.0.2"]

rGREAT_env = new.env()
rGREAT_env$LAST_REQUEST_TIME = 0



# == title
# Global parameters for rGREAT
#
# == param
# -... Arguments for the parameters, see "details" section
# -RESET Reset to default values.
# -READ.ONLY Please ignore.
# -LOCAL Pllease ignore.
# -ADD Please ignore.
# 
# == details
# There are following parameters:
# 
# -``verbose`` Whether to show messages.
#
# == example
# great_opt
great_opt = function(..., RESET = FALSE, READ.ONLY = NULL, LOCAL = FALSE, ADD = FALSE) {}
great_opt = setGlobalOptions(
	verbose = TRUE,
	test = list(
		.value = FALSE,
		.visible = FALSE
	)
)
