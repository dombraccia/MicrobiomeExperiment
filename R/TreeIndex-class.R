#' Graph implementation to query hierarchical feature data
#'
#' Used to manage aggregation and range queries from the Metaviz app UI.
#' @import S4Vectors
setClass("TreeIndex",
         slots = list(
             feature_order = "ANY",
             leaf_of_table = "ANY",
             hierarchy_tree = "ANY",
             node_ids_table = "ANY",
             nodes_table = "ANY",
             .hierarchy="ANY"
         ),
         contains = c("DataFrame")
)

.generate_hierarchy_tree <- function(hierarchy, feature_order) {
    fd <- hierarchy
    for( i in seq(ncol(fd))){
        fd[,i] = as.character(fd[,i])
    }
    hierarchy <- fd

    replacing_na_obj_fData <- hierarchy[,feature_order]

    nas_replaced <- .replaceNAFeatures(replacing_na_obj_fData, feature_order)
    obj_fData <- as.data.table(nas_replaced)
    cols <- feature_order[1:length(feature_order)-1]
    order <- rep(1, length(feature_order)-1)
    ordered_fData <- setorderv(obj_fData, cols = cols, order = order)

    otu_indexes <- seq(1:length(ordered_fData[,get(feature_order[length(feature_order)])]))
    ordered_fData <- ordered_fData[, otu_index:=otu_indexes]
    ordered_fData_df <- as.data.frame(ordered_fData)

    if(length(unique(ordered_fData_df[,1])) > 1){
        allFeatures <- rep("AllFeatures", nrow(ordered_fData_df))
        ordered_fData_df <- cbind(allFeatures, ordered_fData_df)
        feature_order <- unlist(c("allFeatures", feature_order))
    }

    ordered_fData_df
}

.generate_node_ids <- function(hierarchy_tree, feature_order) {
    table_node_ids <- hierarchy_tree
    id_list <- sapply(feature_order, function(level) {
        depth <- which(feature_order == level)
        temp_level <- data.table(table_node_ids[, c(level, "otu_index")])
        temp_level_count <- temp_level[, .(leaf_index = .I[which.min(otu_index)], count = .N), by=eval(level)]

        level_features <- as.character(table_node_ids[[level]])
        for(i in seq_len(nrow(temp_level_count))) {
            row <- temp_level_count[i,]
            if(depth==1 && i == 1){
                id <- paste(depth-1, 0, sep="-")
            } else{
                id <- paste(depth-1, paste(digest(row[,1], algo="crc32"), i, sep=""), sep="-")
            }
            level_features <- replace(level_features, which(level_features == row[[level]]), id)
        }
        level_features
    })

    node_ids_dt <- as.data.table(id_list)
    node_ids_dt$otu_index <- as.character(table_node_ids$otu_index)

    node_ids_table <- node_ids_dt
}

.generate_nodes_table <- function(hierarchy_tree, node_ids_table, feature_order) {
    lineage_DF <- as.data.frame(node_ids_table)
    lineage_table <- node_ids_table
    lineage_DF[,feature_order[1]] <- lineage_table[,get(feature_order[1])]

    for(i in seq(2,length(feature_order))){
        lineage_DF[,feature_order[i]] <- paste(lineage_DF[,feature_order[i-1]], lineage_table[,get(feature_order[i])], sep=",")
    }
    lineage_DT <- as.data.table(lineage_DF)

    root_parents <- rep("None", length(node_ids_table[,get(feature_order[1])]))
    nodes_tab <- data.frame(id = node_ids_table[,get(feature_order[1])], parent = root_parents,
                            lineage = node_ids_table[,get(feature_order[1])],
                            node_label = hierarchy_tree[,1], level = rep(0, length(hierarchy_tree[,1])))

    for(i in seq(2, length(feature_order))){
        temp_nodes_tab <- data.frame(id = node_ids_table[,get(feature_order[i])],
                                     parent = node_ids_table[,get(feature_order[i-1])],
                                     lineage = lineage_DT[,get(feature_order[i])],  node_label = hierarchy_tree[,i],
                                     level = rep(i-1, length(hierarchy_tree[,i])))

        nodes_tab <- rbind(nodes_tab[rownames(unique(nodes_tab[,c("id","parent")])),], temp_nodes_tab[rownames(unique(temp_nodes_tab[,c("id","parent")])),])
    }

    ret_table <- as.data.table(nodes_tab)
    ret_table <- ret_table[,id:=as.character(id)]
    ret_table <- ret_table[,parent:=as.character(parent)]
    ret_table <- ret_table[,lineage:=as.character(lineage)]
    ret_table <- ret_table[,node_label:=as.character(node_label)]
    ret_table <- ret_table[,level:=as.integer(level)]

    ret_table <- ret_table[order(parent)]
    parent_list <- ret_table[,parent]
    orders <- rep(1, length(parent_list))

    for(j in seq(2, length(parent_list))){
        if(parent_list[j] == parent_list[j-1]){
            orders[j] = orders[j-1]+1
        }
    }
    ret_table[,order:=orders]

    ret_table
}

.generate_leaf_of_table <- function(hierarchy_tree, node_ids_table, nodes_table, feature_order) {
    temp_hiearchy_DT <- as.data.table(hierarchy_tree)
    num_features <- length(feature_order)
    hiearchy_cols <- colnames(hierarchy_tree)

    melt_res <- melt(temp_hiearchy_DT, id.vars = c(feature_order[num_features], "otu_index"),
                     measure.vars = c(hiearchy_cols[1:(length(hiearchy_cols)-1)]))
    label_table <- melt_res[,c(1,2,4)]
    setnames(label_table, c("leaf", "otu_index","node_label"))

    label_table <- label_table[,leaf:=as.character(leaf)]
    label_table <- label_table[,otu_index:=as.character(otu_index)]

    lineage_DF <- as.data.frame(node_ids_table)
    lineage_table <- node_ids_table
    lineage_DF[,feature_order[1]] <- lineage_table[,get(feature_order[1])]

    for(i in seq(2,length(feature_order))){
        lineage_DF[,feature_order[i]] <- paste(lineage_DF[,feature_order[i-1]], lineage_table[,get(feature_order[i])], sep=",")
    }
    lineage_DT <- as.data.table(lineage_DF)

    melt_res_lineage <- melt(lineage_DT, id.vars = c(feature_order[num_features], "otu_index"), measure.vars = c(hiearchy_cols[1:(length(hiearchy_cols))-1]))

    lineage_leaf_of_table <- unique(melt_res_lineage[,c(2,4)])
    setnames(lineage_leaf_of_table, c("otu_index","lineage"))

    lineage_leaf_of_table <- lineage_leaf_of_table[,otu_index:=as.character(otu_index)]

    lineage_df <- as.data.frame(lineage_leaf_of_table)
    leaf_node_label <- as.data.frame(label_table)[,c("leaf", "node_label")]

    ret_table <- as.data.table(cbind(lineage_df, leaf_node_label))

    leaf_of_table <- ret_table

    leaf_of_table <- merge(unique(nodes_table[,mget(c("lineage", "id"))]),
                           unique(leaf_of_table) , by="lineage")
    leaf_of_table[,id:=as.character(id)]
}

.replaceNAFeatures = function(replacing_na_obj_fData, feature_order) {

    for(i in seq(1, length(feature_order))){
        na_indices <- which(is.na(replacing_na_obj_fData[,feature_order[i]]))
        for(j in seq(1, length(na_indices))){
            if(i > 1) {
                replacing_na_obj_fData[,feature_order[i]][na_indices[j]] <- paste("Not_Annotated", feature_order[i], replacing_na_obj_fData[,feature_order[1]][na_indices[j]], sep="_")
            } else {
                replacing_na_obj_fData[,feature_order[i]][na_indices[j]] <- paste("Not_Annotated", feature_order[i], sep="_")
            }
        }
        na_indices <- which(replacing_na_obj_fData[,feature_order[i]] == "NA")
        for(j in seq(1, length(na_indices))){
            if(i > 1){
                replacing_na_obj_fData[,feature_order[i]][na_indices[j]] <- paste("Not_Annotated", feature_order[i], replacing_na_obj_fData[,feature_order[1]][na_indices[j]], sep="_")
            } else{
                replacing_na_obj_fData[,feature_order[i]][na_indices[j]] <- paste("Not_Annotated", feature_order[i], sep="_")
            }
        }
        null_indices <- which(replacing_na_obj_fData[,feature_order[i]] == "NULL")
        for(j in seq(1, length(null_indices))){
            if(i > 1){
                replacing_na_obj_fData[,feature_order[i]][null_indices[j]] <- paste("Not_Annotated", feature_order[i], replacing_na_obj_fData[,feature_order[1]][null_indices[j]], sep="_")
            } else{
                replacing_na_obj_fData[,feature_order[i]][null_indices[j]] <- paste("Not_Annotated", feature_order[i], sep="_")
            }
        }
    }

    replacing_na_obj_fData
}


#' Method to initialize TreeIndex Class
#' @export
TreeIndex <- function(hierarchy = NULL, feature_order = NULL){

    if(is.null(hierarchy)) {
        return(new("TreeIndex",
           feature_order = data.frame(),
           leaf_of_table = data.frame(),
           hierarchy_tree = data.frame(),
           node_ids_table = data.frame(),
           nodes_table = data.frame()
        ))
    }

    if(is.null(feature_order)) {
        feature_order <- colnames(hierarchy)
    }

    hierarchy_tree <- .generate_hierarchy_tree(hierarchy, feature_order)
    node_ids_table <- .generate_node_ids(hierarchy_tree, feature_order)
    nodes_table <- .generate_nodes_table(hierarchy_tree, node_ids_table, feature_order)
    leaf_of_table <- .generate_leaf_of_table(hierarchy_tree, node_ids_table, nodes_table, feature_order)

    new("TreeIndex",
        feature_order = feature_order,
        leaf_of_table = leaf_of_table,
        hierarchy_tree = hierarchy_tree,
        node_ids_table = node_ids_table,
        nodes_table = nodes_table,
        .hierarchy = hierarchy
    )
}

#' @export
setMethod("[", "TreeIndex",
          function(x, i, j, ..., drop = FALSE) {
              sHierarchy <- x@.hierarchy[i, j, ..., drop=drop]
              x@hierarchy_tree <- .generate_hierarchy_tree(sHierarchy, x@feature_order)
              x@node_ids_table <- .generate_node_ids(x@hierarchy_tree, x@feature_order)
              x@nodes_table <- .generate_nodes_table(x@hierarchy_tree, x@node_ids_table, x@feature_order)
              x@leaf_of_table <- .generate_leaf_of_table(x@hierarchy_tree, x@node_ids_table, x@nodes_table, x@feature_order)
              # x@.hierarchy <- sHierarchy

              new("TreeIndex",
                  feature_order = x@feature_order,
                  leaf_of_table = x@leaf_of_table,
                  hierarchy_tree = x@hierarchy_tree,
                  node_ids_table = x@node_ids_table,
                  nodes_table = x@nodes_table,
                  .hierarchy = sHierarchy
              )
          }
)

#' @export
setGeneric("getNodes", signature = "x",
           function(x, ...) standardGeneric("getNodes"))

#' @export
setMethod("getNodes", "TreeIndex",
          function(x, selectedLevel=3, start=1, end=1000) {
              nodes_at_level <- x@nodes_table[level==selectedLevel, ]
              nodes_at_level_ids <- nodes_at_level[,id]
              unique(data.frame(ids=nodes_at_level_ids, names=nodes_at_level$node_label))
          })

#' @export
setGeneric("getNodeStates", signature = "x",
           function(x, ...) standardGeneric("getNodeStates"))

#' @export
setMethod("getNodeStates", "TreeIndex",
          function(x) {
              return(list("removed"=0, "expanded"=1, "aggregate"=2))
          })

#' @export
setGeneric("getIndices", signature = "x",
           function(x, ...) standardGeneric("getIndices"))

#' @export
setMethod("getIndices", "TreeIndex",
          function(x, selectedLevel=3, selectedNodes=NULL, start=1, end=1000, format="list") {
              nodes_at_level <- x@nodes_table[level==selectedLevel, ]
              nodes_at_level_ids <- nodes_at_level[,id]

              if(!is.null(selectedNodes) && !(length(selectedNodes) == 0)){
                  nodes_at_level_selections <- rep(2, length(nodes_at_level_ids))
                  names(nodes_at_level_selections) <- nodes_at_level_ids
                  selectedNodes <- c(selectedNodes, nodes_at_level_selections)

                  expand_selections <- which(selectedNodes == 1)
                  if(!is.null(expand_selections) && length(expand_selections) > 0){
                      expand_selection_indices = which(names(selectedNodes) %in% names(expand_selections))
                      selectedNodes <- selectedNodes[-expand_selection_indices]
                  }

                  expanded_children <- x@nodes_table[parent %in% names(expand_selections),id]

                  child_lineage <- x@nodes_table[id %in% names(selectedNodes),]
                  remove_selections <- which(selectedNodes == 0)
                  if(length(remove_selections) > 0){
                      kept_nodes <- child_lineage[!grepl(paste(paste(names(remove_selections), collapse=",|"), ",",sep=""), lineage),]
                      kept_nodes <- kept_nodes[!(id %in% names(remove_selections)),]
                  } else {
                      kept_nodes <- child_lineage
                  }

                  agg_selections <- which(selectedNodes == 2)
                  if(length(agg_selections) > 0){
                      kept_nodes <- kept_nodes[!grepl(paste(paste(names(agg_selections), collapse=",|"), ",",sep=""), lineage),]
                  }
                  kept_nodes <- as.character(kept_nodes[,id])
                  nodes_at_level <- x@nodes_table[id %in% c(kept_nodes,expanded_children),]
              }

              leaf_order_table <- as.data.table(x@hierarchy_tree[,c(x@feature_order[length(x@feature_order)], "otu_index")])
              setnames(leaf_order_table, c("leaf", "otu_index"))
              leaf_order_table <- leaf_order_table[,leaf:=as.character(leaf)]
              leaf_order_table <- leaf_order_table[otu_index >= start & otu_index <= end]

              leaf_indices <- merge(leaf_order_table, merge(nodes_at_level, x@leaf_of_table, by="id"), by="leaf")
              setorderv(leaf_indices, "otu_index.x")
              leaf_indices$lineage.y <- NULL
              leaf_indices$otu_index.y <- NULL
              leaf_indices$node_label.y <- NULL
              colnames(leaf_indices) <- c("leaf", "otu_index", "id", "parent", "lineage", "node_label", "level", "order")

              if(format == "aggTable") {
                  return(leaf_indices)
              }
              else if(format == "dataframe"){
                  groups <- leaf_indices[, .(indices=paste0(otu_index, collapse=","), leaf_nodes=paste0(leaf, collapse=",")), by=.(id, parent, lineage, node_label, level, order)]
                  return(groups)
              }
              else if(format == "list") {
                  groups <- leaf_indices[, .(indices=paste0(otu_index, collapse=","), leaf_nodes=paste0(leaf, collapse=",")), by=.(id, parent, lineage, node_label, level, order)]
                  nodes <- as.list(groups$indices)
                  names(nodes) <- groups$node_label
                  return(nodes)
              }
              return(leaf_indices)
          }
)

setAs("DataFrame", "TreeIndex", function(from) {
    ## TODO: add coercion method for DataFrame to TreeIndex
    from@.hierarchy
})

#' @keywords internal
#' @importFrom methods callNextMethod
#' @export
setMethod("show", "TreeIndex", function(object) {
    cat(
        "Tree Index",
        "with height:", length(object@feature_order),
        "\n Tree levels:", paste(object@feature_order, collapse=" -> "),
        "\n Leaf nodes:", nrow(object@.hierarchy),
        sep=" "
    )
})
