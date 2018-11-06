## Class constructor
.MicrobiomeExperiment <- setClass("MicrobiomeExperiment",
                                  contains = "SummarizedExperiment",
                                  representation(rowData = "TreeIndex"))

#' The MicrobiomeExperiment representation class
#'
#' SummarizedExperiment-like class for microbiome data. rowData is
#' a MicrobiomeFeatures object LINK so it includes: a taxonomy table
#' (DataFrame), #' optional phylogentic tree (phylo object), and a sequence
#' database.
#'
#' It supports (most) of the interface to phyloseq objects
#'
#' @include TreeIndex-class.R
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#'
#' @examples
#'
#' library(metagenomeFeatures)
#' data(mock_mgF)
#'
#' sampleNames <- letters[1:4]
#' pd <- DataFrame(a=letters[1:4], b=1:4)
#' numcounts <- nrow(mock_mgF) * 4
#' counts <- matrix(sample(1:1000,numcounts,replace=TRUE), nr=nrow(mock_mgF), nc=4)
#'
#' MicrobiomeExperiment(assays=SimpleList(counts=counts),
#' rowData=mock_mgF,
#' colData=pd
#' )
#'
#' @aliases MicrobiomeExperiment-class
#' @export
MicrobiomeExperiment <- function(assays = SimpleList(),
                                 rowData = TreeIndex(),
                                 ...) {
    if (is.data.frame(rowData))
        rowData <- TreeIndex(rowData)

    SummarizedExperiment <-
        if (!is(assays, "SummarizedExperiment"))
            SummarizedExperiment(assays = assays, ...)
    else
        assays

    .MicrobiomeExperiment(SummarizedExperiment, rowData = rowData)
}

#' @importFrom methods callNextMethod
#' @importFrom SummarizedExperiment assays rowData colData
#  @importFrom BiocGenerics normalize
#' @importFrom S4Vectors metadata
#' @export
setMethod("[", signature(x = "MicrobiomeExperiment"),
          function(x, i, j) {
              obj <- callNextMethod()
              counts <- assays(x)$counts

              if (!missing(i)) {getIndices
                  i <- as.vector(i)
                  rowData <- x@rowData[i, ]
                  counts <- counts[i, ]
              }

              if (!missing(j)) {
                  j <- as.vector(j)
                  colData <- x@colData[j, ]
                  counts <- counts[, j]
              }

              MicrobiomeExperiment(SimpleList(counts = counts), rowData = rowData)
          })

#' @export
setMethod("rowData", signature("MicrobiomeExperiment"),
          function(x) {
              x@rowData
          })

#' @export
setMethod("show", signature("MicrobiomeExperiment"),
          function(object) {
              cat("class: MicrobiomeExperiment \n", sep=" ")
              cat("dim:", nrow(object), ncol(object), "\n", sep=" ")
              cat("metadata:\n")
              cat(show(metadata(object)))
              cat("rowData:\n")
              cat(show(rowData(object)), "\n")
              cat("colData:\n")
              cat(show(colData(object)))
          })

#' @export
setGeneric("aggregateAt", signature = "x",
           function(x, ...) standardGeneric("aggregateAt"))

#' @export
setMethod("aggregateAt", "MicrobiomeExperiment",
          function(x, samples=NULL, selectedLevel=3, selectedNodes=NULL, aggFun=colSums, start=1, end=1000, format="MicrobiomeExperiment") {

              if(is.null(samples)) {
                  samples <- colnames(x)
              }

              groups <- getIndices(rowData(x), selectedLevel=3, selectedNodes=NULL, start=1, end=1000, format="list")
              counts <-  assays(x)$counts

              newMat <- array(NA,dim=c(length(groups), ncol(x)))
              for(i in seq_along(groups)) {
                  indices <- as.integer(strsplit(groups[[i]], ",")[[1]])
                  newMat[i,] = aggFun(counts[indices,])
              }

              rownames(newMat) <- names(groups)
              colnames(newMat) <- colnames(x)


              if(format == "MicrobiomeExperiment") {
                subtree <- getIndices(rowData(x), selectedLevel=3, selectedNodes=NULL, start=1, end=1000, format="TreeIndex")
                return(MicrobiomeExperiment(SimpleList(counts = newMat), rowData = subtree))
              }
              else {
                return(newMat)
              }
          })
