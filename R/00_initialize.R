initialize <-
function(
    path_to_fcs=FCSpath,
    path_to_output=Outpath,
    verbose=TRUE
      ){
        
        # if(path_to_intermediary_results==tempdir()){
        #   tmpdir = path_to_intermediary_results
        #   if(dir.exists(tmpdir)){
        #     i = 1
        #     while(dir.exists(file.path(tmpdir, i))){
        #       i = i + 1
        #     }
        #     path_to_intermediary_results = file.path(tmpdir, i)
        #   }
        #   if(verbose){
        #     message("Using ", tmpdir, " temporary directory to store intermediary results as no non-temporary directory has been specified")
        #   }
        # }
        
        ## The paths below have to point to directories. If they do not exist the script will create them to store outputs
        # path_to_subsetted_fcs <- file.path(path_to_intermediary_results,"subsetted_fcs") ## Subsetted here means for instance "gated on live singlets CD45+"
        path_to_intermediary <- file.path(path_to_output,"/intermediary")
        path_to_graph <- file.path(path_to_output,"/graph")
        path_to_downstream <- file.path(path_to_output,"/downstream")
        # ## A CSV file to map files to PE targets
        # path_to_annotation_file <- file.path(path_to_intermediary_results,"annotation.csv") ## Has to be a comma-separated csv file, with two columns. The first column has to be the name of the FCS files and the second the marker bound to the PE reporter. The first line (column names) have to be "file" and "target". You can leave the unlabelled PEs empty
        
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
