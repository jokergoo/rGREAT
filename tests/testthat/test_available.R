context("Test someting")

library(rGREAT)

oe = try(httpHEAD(gsub("public/cgi-bin$", "", rGREAT:::BASE_URL_LIST[["default"]]), .opts = list(connecttimeout = 10)), silent = TRUE)
connected = TRUE
if(inherits(oe, "try-error")) {
	connected = FALSE
}


# test will be added later
test_that("Test something", {
    
    expect_that(1, is_identical_to(1))

})
