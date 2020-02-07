context("Test clmi")


test_that("clmi throws errors correctly", {
    lod.var <- toy_data$lod
    expect_error(clmi("poll ~ case_cntrl + smoking + gender, toy_data", lod, 1))
    expect_error(clmi(poll ~ case_cntrl + smoking + gender, toy_data, lodu, 1))
    expect_error(clmi(poll ~ case_cntrl + smoking + gender, toy_data, NULL, 1))

    expect_error(clmi(a * poll ~ case_cntrl + smoking + gender, toy_data,
                      lod, 1))

    a <- 2
    fn <- function(x) a * x
    expect_true(!is.null(clmi(fn(poll) ~ case_cntrl + smoking + gender, toy_data,
                      lod, 1)))

    expect_error(clmi(poll ~ case_cntrl + smoking + gender, "toy_data", lod, 1))
    expect_error(clmi(poll ~ case_cntrl + smoking + gender, toy_data, lod, "a"))
    expect_error(clmi(poll ~ case_cntrl + smoking + gender, toy_data, lod, 1, 0))
    expect_error(clmi(poll ~ case_cntrl + smoking + gender, toy_data, lod,
                      1, c(2, 4)))

    df <- toy_data
    df$lod <- NULL
    expect_error(clmi(poll ~ case_cntrl + smoking + gender, df, lod, 1))

    # Test for improper lod values
    df <- toy_data
    df$poll <- log(df$poll)
    expect_error(clmi(poll ~ case_cntrl + smoking + gender, df, lod, 1))

    # Test that missing covariates throw an error
    df <- toy_data
    df$smoking[1] <- NA
    expect_error(clmi(poll ~ case_cntrl + smoking + gender, df, lod, 1))

    # Test that factors cause an error
    df <- toy_data
    df$gender <- as.factor(df$gender)
    expect_error(clmi(poll ~ case_cntrl + smoking + gender, df, lod, 1))
})
