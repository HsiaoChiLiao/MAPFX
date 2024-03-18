#' @author Hsiao-Chi Liao
#' 
#' @importFrom grDevices dev.off jpeg
#' @importFrom graphics abline
#' @importFrom methods is
#' @importFrom stats as.formula contr.sum contrasts<- dexp dnorm median model.matrix optim pexp pnorm quantile sd setNames
#' @importFrom utils head read.csv
#' 
#' 
initialize <-
function(
    path_to_fcs=path_to_fcs, #FCSpath,
    path_to_output=path_to_output, #Outpath,
    verbose=TRUE
    ){
        
    ## The paths below have to point to directories. If they do not exist the script will create them to store outputs
    path_to_intermediary <- file.path(path_to_output,"/intermediary")
    path_to_graph <- file.path(path_to_output,"/graph")
    path_to_downstream <- file.path(path_to_output,"/downstream")

    ## create folders
    dir.create(path_to_intermediary)
    dir.create(path_to_graph)
    dir.create(path_to_downstream)
        
    paths <- c(
    input=file.path(path_to_fcs,"/"),
    intermediary=file.path(path_to_intermediary,"/"),
    graph=file.path(path_to_graph,"/"),
    downstream=file.path(path_to_downstream,"/")
    )
    
    paths <- vapply(paths, path.expand, "path")
        
    return(list(paths=paths))
    }
