
BASE_URL_LIST = c("default" = "http://great.stanford.edu/public/cgi-bin",
	         "3.0.0" = "http://great.stanford.edu/public-3.0.0/cgi-bin",
	         "2.0.2" = "http://great.stanford.edu/public-2.0.2/cgi-bin")
BASE_URL_LIST["3.0"] = BASE_URL_LIST["3"] = BASE_URL_LIST["3.0.0"]
BASE_URL_LIST["2.0"] = BASE_URL_LIST["2"] = BASE_URL_LIST["2.0.2"]

SPECIES = list("default" = c("hg19", "mm9", "mm10", "danRer7"),
	           "3.0.0" = c("hg19", "mm9", "mm10", "danRer7"),
	           "2.0.2" = c("hg19", "hg18", "mm9", "danRer7"))
SPECIES["3.0"] = SPECIES["3"] = SPECIES["3.0.0"]
SPECIES["2.0"] = SPECIES["2"] = SPECIES["2.0.2"]

rGREAT_env = new.env()
rGREAT_env$LAST_REQUEST_TIME = 0

