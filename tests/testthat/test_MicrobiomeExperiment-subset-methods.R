context("MicrobiomeExperiment subset methods")

test_that("[ subsetting works", {
    library(metagenomeSeq)
    library(data.table)
    library(S4Vectors)
    data(mouseData)
    counts <- MRcounts(mouseData)

    rowData <- TreeIndex(fData(mouseData))
    me <- MicrobiomeExperiment(SimpleList(counts = counts), rowData = rowData)

    # subset rows and columns
    subset_me <- MbExp[1:100,1:5]

    expect_true(is(subset_me, "MicrobiomeExperiment"))
    expect_equal(nrow(subset_me), 100)
    expect_equal(ncol(subset_me), 5)
    expect_equal(nrow(subset_me@rowData@.hierarchy), 100)
    expect_equal(nrow(subset_me@colData), 5)
    expect_equal(nrow(assays(subset_me)$counts), 100)
    expect_equal(ncol(assays(subset_me)$counts), 5)
})
