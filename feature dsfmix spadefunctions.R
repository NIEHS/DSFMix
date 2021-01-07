####### INTEGRATED ALGORITM FOR SPADE FOREST for scRNA-seq #############

##### the first portion of SPADE.driver, used to downsample (if necessary) and cluster data into k (default 200) clusters. Should not really generate any graphs or things, quite yet. Note: only implemented for one file at a time currently.

######>
# default parameters are set to keep as many cells as possible.
# run_spade=TRUE indicates that normal spade graphs will be produced along with tSNE graphs and split-spade trees.
# run_silhouette=TRUE will perform silhouette.
# run_spade_only=FALSE is for testing purposes and stops the run after producing normal spade trees.
# gated=TRUE indicates there is a marker named "type" in the FCS file which contains integers corresponding to the identity of the cell in question.
# remove_single_cell_clusters=TRUE will cause clusters with only 1 cell in them to be retained. This is important if number of cells is comparable to k. This can be safely turned on only for large cell counts, otherwise it is better to leave off as in default.
# perplexity by default is k/20 (Not needed for spade forest)
# costtol, e.g. when running tSNE the run stops when decrease in error is less than 10e-(costtol). A higher value for costtol increases run time and the separation between populations. (Not needed for spade forest)
# cofactor must be specified if it is to be passed to t-SNE (Not needed for spade forest)


###### Algorithm steps ###############

#1. check that fewer than 50000 cells are in FCS file. If yes, stop.
#2. if "gated" parameter is false, create a new FCS file "firstout" with an extra column called "type" which will contain all 1s (if "gated" is true, the type column should already exist and should contain integers corresponding to cell types)
#3. if a directory with more than 1 file is specified, stop.
#4. if "panels" are specified, stop.
#5. remove existing Density and Cluster Columns from the file.
#6. calculate new density values.
#7. if "downsample" is true, downsample as normal (this should be false). If "downsample" is false, skip it, but add columns to match.
#8. Cluster files using hierarchichical clustering as in SPADE, then generate a minimum spanning tree and graph files.
#9. if "run_silhouette" is true, perform silhouette filtering first by reassigning clusters and then by purifying each cluster (this should be set to false though).
#10. Make a new copy of the FCS file with a cluster assignment column
#11. Calculate how well cell types were segregated into clusters. Find the cell type which is most common in the cluster and the portion of the cells in the cluster which belong to that cell type. Assign this cluster to the "max_type" column in a new FCS file. Assign the percentage of this type in the cluster to the "max_type" column. Calculate the standard deviation within each cluster for each column and average it. Assign this to the "average_sd" column.
#12. If downsampling was performed within the integrated tool, perform upsampling.
#13. Preparing data for t-SNE run: calculate medians for each cluster, saving as CSV files (filteredclustered.csv contains all columns, filteredclusteredprepared contains only the clustering columns).
#14. Create a new FCS file which has been reduced to clusters, adding the original cell_counts, as well as an integer indicating cluster assignments
#15. If "run_spade" is true, copy all relevant files (FCS files, graph files, etc.) to a new folder and run original SPADE on these files. Note this is performed on data BEFORE it was reduced to clusters
#16. Run TAHC: Hierarchical clusering for spade forest See. published paper for reference
#17. separate regions. Cut the MST graph into dcut number of subtrees
  #18. Creating pdfs for each marker of separate subtrees from original MST representing disjoint tree motifs.


###### INTEGRATED SPADE.DRIVER FUNCTIONS ##############

dsfmixmain <-
function (files, file_pattern = "*.fcs", out_dir = ".", cluster_cols,
panels = NULL, comp = FALSE, arcsinh_cofactor = NULL, transforms = flowCore::arcsinhTransform(a = 0, b = 0.2), downsampling_target_number = NULL, downsampling_target_pctile = 1, downsampling_target_percent = 1, downsampling_exclude_pctile = 0, k = 200, clustering_samples = NULL, layout = igraph:::layout.kamada.kawai, pctile_color = c(0, 1), fcs_channel_mappings_json = NULL, graph_cols = NULL, run_spade=TRUE, do_real_filtering = FALSE, run_spade_only=FALSE, gated=FALSE,forestcluster=NULL, dcut = 3, leafsort=TRUE,minsize=20, remove_single_cell_clusters=FALSE, perplexity = NULL, costtol=7, cofactor=5, downsample=FALSE,reg=NULL,forestlayout=igraph:::layout_in_circle)

{
    
    #1. check that fewer than 50000 cells are in FCS file. If yes, stop.
    
    # check that fewer than 50000 cells are in the file.
    ff1 <- spade:::SPADE.read.FCS(files)
   # if(nrow(exprs(ff1)) >= 50000) {
    #    stop("this implementation is not designed for more than 50000 cells. Continuing to run will produce crippling errors when upsampling.")
    #}
    
    #2. if "gated" parameter is false, create a new FCS file "firstout" with an extra column called "type" which will contain all 1s (if "gated" is true, the type column should already exist and should contain integers corresponding to cell types)
    if (!gated) {
        
        ff1 <- spade:::SPADE.read.FCS(files)
        in_data <- exprs(ff1)
        params <- flowCore::parameters(ff1)
        desc <- description(ff1)
        
        pd <- pData(params)
        typevec <- 1:nrow(in_data)
        typevec[1:nrow(in_data)] <- 1
        
        
        
        firstout <- paste(out_dir, basename(files), ".type.fcs", sep="")
        
        channel_number <- ncol(in_data) + 1
        channel_id <- paste("$P", channel_number, sep = "")
        channel_name <- "type"
        channel_range <- 25
        plist <- matrix(c(channel_name, channel_name, channel_range,
        0, channel_range - 1))
        rownames(plist) <- c("name", "desc", "range", "minRange",
        "maxRange")
        colnames(plist) <- c(channel_id)
        pd <- rbind(pd, t(plist))
        pData(params) <- pd
        out_data <- cbind(in_data, type = typevec)
        out_frame <- flowFrame(out_data, params, description = desc)
        keyval <- list()
        keyval[[paste("$P", channel_number, "B", sep = "")]] <- "32"
        keyval[[paste("$P", channel_number, "R", sep = "")]] <- toString(channel_range)
        keyval[[paste("$P", channel_number, "E", sep = "")]] <- "0,0"
        keyval[[paste("$P", channel_number, "N", sep = "")]] <- channel_name
        keyval[[paste("$P", channel_number, "S", sep = "")]] <- channel_name
        keyword(out_frame) <- keyval
        
        write.FCS(out_frame, firstout)
        files <- firstout
        
    }
    
    
    
    # this is part of the silhouette run. if true clusters will be split into two clusters and the one with the higher silhouette width will be kept.
    
    
    
    
    # if a directory was specified, find FCS files within it. However usage of multiple FCS files will cause errors with this implementation. Combine them before running, and indicate their source with the "type" column.
    
    if (length(files) == 1 && file.info(files)$isdir) {
        files <- dir(spade:::SPADE.strip.sep(files), full.names = TRUE,
        pattern = glob2rx(file_pattern))
    }
    if (length(files) == 0) {
        stop("No input files found")
    }
    
    # combine files if more than one was found. However this simply pools the cells and does not indicate their sources.
    #3. if a directory with more than 1 file is specified, stop.
    
    if (length(files) > 1) {
        stop("only implemented for one file. Combine prior to running this script")
        
        oldf <- files
        
        oldf1 <- spade:::SPADE.read.FCS(oldf[1])
        
        params <- flowCore::parameters(oldf1)
        des <- description(oldf1)
        
        newdata <- NULL
        outloc <- paste(out_dir, "combined.data.fcs", sep="")
        denseloc <- paste(out_dir, "combined.data.density.fcs", sep="")
        downsamploc <- paste(out_dir, "combined.data.downsample.fcs", sep="")
        
        for(f in oldf) {
            f <- spade:::SPADE.read.FCS(f)
            newdata <- rbind(newdata, exprs(f))
            
        }
        new_frame <- flowFrame(newdata, params, des)
        write.FCS(new_frame, outloc)
        
        
        
        SPADE.addDensityToFCS(denseloc, outloc, cols = cluster_cols, comp=comp, transforms = transforms)
        SPADE.downsampleFCS(downsamploc, outloc, downsampling_exclude_pctile, downsampling_target_pctile, downsampling_target_number, downsampling_target_percent)
        
        exclude_pctile <- 0
        target_pctile <- 1
        target_number <- NULL
        target_percent <- 1
        
        
        files <- outloc
        message("files were combined successfully")
    }
    
    
    
    # panels. current implementation does not support panels.
    # much of this section is unchanged from original SPADE
    #4. if "panels" are specified, stop.
    
    out_dir <- spade:::SPADE.normalize.out_dir(out_dir)
    if (!is.null(panels)) { ###delete panel variable
        if (!is.list(panels))
        stop("Invalid panels argument, see function documentation")
        lapply(panels, function(x) {
            if (!is.list(x) || !all(c("panel_files", "median_cols") %in%
            names(x)))
            stop("Invalid panel found, see function documentation for proper panel structure")
            if (!all(x$panel_files %in% basename(files)))
            stop("Panel files must be a subset of analysis files")
            if (!is.null(x$reference_files) && !all(x$reference_files %in%
            x$panel_files))
            stop("Panel reference files must be a subset of panel files")
        })
    }
    if (!is.null(arcsinh_cofactor)) {
        warning("arcsinh_cofactor is deprecated, use transform=flowCore::arcsinhTransform(...) instead")
        transforms <- flowCore::arcsinhTransform(a = 0, b = 1/arcsinh_cofactor)
    }
    #5. remove existing Density and Cluster Columns from the file.
    
    for (f in files) {
        SPADE.removeExistingDensityAndClusterColumns(f)
    }
    
    #6. calculate new density values.
    #7. if "downsample" is true, downsample as normal (this should be false). If "downsample" is false, skip it, but add columns to match.
    
    # return(list(files, cluster_cols, transforms))
    density_files <- c()
    sampled_files <- c()
    for (f in files) {
        message("Prepare to Downsample file: ", f)
        f_density <- paste(out_dir, basename(f), ".density.fcs",
        sep = "")
        f_sampled <- paste(out_dir, basename(f), ".downsample.fcs",
        sep = "")
        
        if(downsample) {
            
            message("Begin Downsampling file: ", f)
            spade::SPADE.addDensityToFCS(f, f_density, cols = cluster_cols,
            comp = comp, transforms = transforms)
            spade::SPADE.downsampleFCS(f_density, f_sampled, exclude_pctile = 	downsampling_exclude_pctile,
            target_pctile = downsampling_target_pctile, target_number = downsampling_target_number,
            target_percent = downsampling_target_percent)
            density_files <- c(density_files, f_density)
            sampled_files <- c(sampled_files, f_sampled)
            
        } else {
            
            print("skipping downsampling")
            
            ff <- spade:::SPADE.read.FCS(f)
            params <- flowCore::parameters(ff)
            pd <- pData(params)
            in_data <- exprs(ff)
            desc <- description(ff)
            
            dd <- 1:nrow(in_data)
            
            channel_number <- ncol(in_data) + 1
            channel_id <- paste("$P", channel_number, sep = "")
            channel_name <- "density"
            channel_range <- nrow(in_data) + 1
            plist <- matrix(c(channel_name, channel_name, channel_range,
            0, channel_range - 1))
            rownames(plist) <- c("name", "desc", "range", "minRange",
            "maxRange")
            colnames(plist) <- c(channel_id)
            pd <- rbind(pd, t(plist))
            pData(params) <- pd
            out_data <- cbind(in_data, density = dd)
            out_frame <- flowFrame(out_data, params, description = desc)
            keyval <- list()
            keyval[[paste("$P", channel_number, "B", sep = "")]] <- "32"
            keyval[[paste("$P", channel_number, "R", sep = "")]] <- toString(channel_range)
            keyval[[paste("$P", channel_number, "E", sep = "")]] <- "0,0"
            keyval[[paste("$P", channel_number, "N", sep = "")]] <- channel_name
            keyval[[paste("$P", channel_number, "S", sep = "")]] <- channel_name
            keyword(out_frame) <- keyval
            write.FCS(out_frame, f_density)
            write.FCS(out_frame, f_sampled)
            density_files <- c(density_files, f_density)
            sampled_files <- c(sampled_files, f_sampled)
            
            
        }
    }

    
    params <- NULL
    desc <- NULL
    
    h <- spade:::SPADE.read.FCS(sampled_files)
    params <- flowCore::parameters(h)
    params1 <- params
    desc <- description(h)
    
    #8. Cluster files using hierarchichical clustering as in SPADE, then generate a minimum spanning tree and graph files.
    message("Clustering files...")
    cells_file <- paste(out_dir, "clusters.fcs", sep = "") # will contain downsampled cells, their cluster assignments, and only clustering markers.
    clust_file <- paste(out_dir, "clusters.table", sep = "")
    graph_file <- paste(out_dir, "mst.gml", sep = "")
    
    # cluster cells and generate files for upsampling. slight modifications have been made to the SPADE.FCSToTree and SPADE.cluster functions so they do not drop clusters with a single cell.
    print("code still running.....")
    gg = spadefcstotree2(sampled_files, cells_file, graph_file, clust_file, cols = cluster_cols, transforms = transforms, k = k, desired_samples = clustering_samples, comp = comp, remove_single_cell_clusters=remove_single_cell_clusters)
    
    print("code still running.....")
    # return(gg)
    clust=gg[[1]]
    k=gg[[2]]
    print(paste("Number of spade clusters =",k))
    # centers = gg[[3]]
    # first cluster assignments have been obtained.
    # keep the original downsampled file, just in case.
    
    downsampled_files <- sampled_files
    origdown <- spade:::SPADE.read.FCS(downsampled_files)
    
    downsampled_out <- paste(out_dir, "downsampled.original.fcs", sep="")
    write.FCS(origdown, file=downsampled_out)
    
    hclustered = clust$hclust
    assign = clust$assign
    #d = clust$distance
    
    print("cluster object obtained")
    
    # cluster assignments
    groups <- assign
    #9. if "run_silhouette" is true, perform silhouette filtering first by reassigning clusters and then by purifying each cluster (this should be set to false though).
    
    # if silhouette is to be run....
    if(do_real_filtering) {
        
        # run silhouette as normal. Then optionally purify the re-assigned clusters using more clustering, on a per-cluster basis.
        
       # silhouette_out = paste(out_dir, "silhouette.pdf")
       # silout = ccast_silhouette2(d=d, fileout=silhouette_out, nCluster.init = k, assign=assign, iter.max = 5, diagnosis = TRUE)
        
        #groups = silout$optFitCluster[["cluster"]]
        
        in_fcs <- spade:::SPADE.read.FCS(downsampled_out)
        groups2 <- NULL
        
        data <- c()
        in_data <- exprs(in_fcs)
        params <- flowCore::parameters(in_fcs)
        pd <- pData(params)
        
        if (is.null(cluster_cols)) {
            cols <- as.vector(pd$name)
        } else {
            cols<-cluster_cols
        }
        idxs <- match(cols, pd$name)
        if (any(is.na(idxs))) {
            stop("Invalid column specifier")
        }
        data <- rbind(data, in_data[, idxs, drop = FALSE])
        colnames(data) <- pd$name[idxs]
        print(paste("dimension of unfiltered data =", dim(data)))
        # we run silhouette 5 times to re-arrange clusters of cells. Now it is time to go into each cluster one by one. We will hierarchichically cluster the cells in each cluster into two groups, and then keep the cluster with less variance. This will make each cluster exceptionally pure.
            
            original.ids <- 1:nrow(data)
            
            # iterate through each cluster
            is.na(groups) <- which(groups == 0)
            
            all.ids <- NULL
            
            cluster.count <- 1 # because sometimes clusters are lost.
            # This ensures the numberings in "groups2" are correct
            message("filtering clusters")
            for (i in c(1:max(groups, na.rm = TRUE))) {
                
                # obs contains the original row indices of the group of cells
                obs <- which(groups == i)
                #obs2 <- cbind(obs, 1:length(obs))
                # cluster.cell.set <- NULL
                if (length(obs) > 3) {
                    cells <- data[obs,]
                    
                    # now we need to obtain the row indices of
                    
                    idset <- NULL
                    #d <- dist(cells)
                    #h <- hclust2(d)
                    h<- hclust2(objects=as.matrix(cells), thresholdGini=0.2)
                    treec <- cutree(h, k=2)
                
                    # need to separate each cluster in two. One bad, one good.
                    # First perform hierarchical clustering on that cluster with k=2.
                    # If there is a homogeneous group, it should all be found in one of the clusters.
                    # Silhouette can be used to evaluate the clusters and choose the better one.
                    
                    idset=which(treec %in% Mode(treec))
                    ##all.ids <- c(all.ids, obs[which(idset %in% obs2[,2])])
                    all.ids <- c(all.ids, obs[idset])
                    
                    j <- 1:length(idset)
                    j[1:length(idset)] <- cluster.count
                    
                    groups2 <- c(groups2, j)
                    cluster.count <- cluster.count + 1
                }
                
            }
            # reassemble data
            new.data <- in_data[all.ids,]
            k <- max(groups2)
            groups<- groups2
        
        save(new.data,file="new.data.rdata")
        print(paste("dimension of filtered data =", dim(new.data)))
        # generate a new graph file
        
        centers <- NULL
            is.na(groups2) <- which(groups2 == 0)
            for (i in c(1:max(groups2, na.rm = TRUE))) {
                obs <- which(groups2 == i)
                if (length(obs) > 1) {
                    centers <- rbind(centers, colMeans(new.data[obs,cluster_cols , drop = FALSE]))
                    # groups2[obs] <- nrow(centers)
                } else {
                        is.na(groups2) <- obs
                }
                
            }
        rownames(centers)=1:max(groups2)
        save(centers,file="centers.rdata")
        spade:::SPADE.writeGraph(spade:::SPADE.clustersToMST(centers), graph_file)

     #  makegraph(new.data[,cluster_cols],  graph_file, run_silhouette, transforms, groups2)
        
    }
    
    ## creates a new FCS file from the downsampled data but with the silhouette-modified clusters. So that upsampling can happen properly
    
    downsampledsilloc<-paste(out_dir, "downsampled.silclustered.fcs", sep="")
    

    print("trying to create new downsampled clustered file")
    
    # adding in a marker which contains cluster assignment
    #10. Make a new copy of the FCS file with a cluster assignment column
    pd <- pData(params1)
    
    print(paste("do_real_filtering is",do_real_filtering))
    if(!do_real_filtering) {
    
        in_data <- exprs(spade:::SPADE.read.FCS(downsampled_files))
        new.data <- in_data
             k <- max(groups)
        
         print(paste("dimension of filtered data =", dim(new.data)))
         # generate a new graph file
         centers <- NULL
             is.na(groups) <- which(groups == 0)
             for (i in c(1:max(groups, na.rm = TRUE))) {
                 obs <- which(groups == i)
                 if (length(obs) > 1) {
                     centers <- rbind(centers, colMeans(new.data[obs,cluster_cols , drop = FALSE]))
                     # groups2[obs] <- nrow(centers)
                 } else {
                         is.na(groups) <- obs
                 }

             }
        rownames(centers)=1:max(groups)
         save(new.data,file="new.data.rdata")
         save(groups,file="groups.rdata")
         save(centers,file="centers.rdata")
    } else {
        in_data <- new.data
    }
    # return(list(new.data, groups))
    
    channel_number <- ncol(in_data) + 1
    channel_id <- paste("$P", channel_number, sep = "")
    channel_name <- "cluster"
    channel_range <- k + 1
    plist <- matrix(c(channel_name, channel_name, channel_range,
    0, channel_range - 1))
    rownames(plist) <- c("name", "desc", "range", "minRange",
    "maxRange")
    colnames(plist) <- c(channel_id)
    pd <- rbind(pd, t(plist))
    pData(params1) <- pd
    out_data <- cbind(in_data, cluster = groups)
    out_frame <- flowFrame(out_data, params1, description = desc)
    keyval <- list()
    keyval[[paste("$P", channel_number, "B", sep = "")]] <- "32"
    keyval[[paste("$P", channel_number, "R", sep = "")]] <- toString(channel_range)
    keyval[[paste("$P", channel_number, "E", sep = "")]] <- "0,0"
    keyval[[paste("$P", channel_number, "N", sep = "")]] <- channel_name
    keyval[[paste("$P", channel_number, "S", sep = "")]] <- channel_name
    keyword(out_frame) <- keyval
    write.FCS(out_frame, downsampledsilloc)
    
    message("created downsampled file with correct clusters")
    message("No need to try to add 'percentage max cell type'")
    
    # upsampling as in SPADE
    #12. If downsampling was performed within the integrated tool, perform upsampling.
    sampled_files <- downsampledsilloc
    #if(downsample) {
     #   sampled_files <- c()
      #  for (f in density_files) {
      #      message("Upsampling file: ", f)
      #      f_sampled <- paste(f, ".cluster.fcs", sep = "")
      #      SPADE.addClusterToFCS2(f, f_sampled, cells_file, cols = cluster_cols,
       #     transforms = transforms, comp = comp)
       #     sampled_files <- c(sampled_files, f_sampled) #the FCS file with density and clusters added
      #  }
    #} else {
      #  sampled_files <- downsampledsilloc
   # }
    #message("file upsampled.")
    
    # now the cluster assignments are properly placed as the final column of "sampled_files". Groups is not necessary, but finding the medians is.
    
    #13. Preparing data for SPADE forest run: calculate medians for each cluster, saving as CSV files (filteredclustered.csv contains all columns, filteredclusteredprepared contains only the clustering columns).
    
    #q will contain the matrix of medians for input into SPADE forest
    q <- convertdata(sampled_files, k)
    save(q,file="q.rdata")
    # ready matrices
    
    output1 = paste(out_dir, "filteredclustered.csv", sep="")
   output2 = paste(out_dir, "filteredclusteredprepared.csv", sep="")
   qnotrans = q
   cell_count <- qnotrans[,ncol(qnotrans)]
   qnotrans <- qnotrans[,-ncol(qnotrans)]
    
    #make a new FCS file which will contain the properly clustered and transformed data
   # add the "cell_count" column to FCS
    #14. Create a new FCS file which has been reduced to clusters, adding the original cell_counts, as well as an integer indicating cluster assignments
    
   newframeloc<-paste(out_dir, "newframe.fcs", sep="")
    
   pd <- pData(params)
   in_data <- qnotrans
   channel_number <- ncol(in_data) + 1
   channel_id <- paste("$P", channel_number, sep = "")
   channel_name <- "cell_count"
   channel_range <- max(cell_count) + 1
   plist <- matrix(c(channel_name, channel_name, channel_range,
   0, channel_range - 1))
  rownames(plist) <- c("name", "desc", "range", "minRange","maxRange")
  colnames(plist) <- c(channel_id)
    pd <- rbind(pd, t(plist))
    pData(params) <- pd
    out_data <- cbind(in_data, cell_count = cell_count)
    
    out_frame <- flowFrame(out_data, params, description = desc)
    keyval <- list()
    keyval[[paste("$P", channel_number, "B", sep = "")]] <- "32"
    keyval[[paste("$P", channel_number, "R", sep = "")]] <- toString(channel_range)
    keyval[[paste("$P", channel_number, "E", sep = "")]] <- "0,0"
    keyval[[paste("$P", channel_number, "N", sep = "")]] <- channel_name
    keyval[[paste("$P", channel_number, "S", sep = "")]] <- channel_name
    keyword(out_frame) <- keyval
    
    write.FCS(out_frame, newframeloc)
    
    
    # add cluster column
    in_fcs <- spade:::SPADE.read.FCS(newframeloc)
    params2 <- flowCore::parameters(in_fcs)
    desc <- description(in_fcs)
    pd <- pData(params2)
    in_data <- exprs(in_fcs)
    
    cluster_data <- 1:nrow(in_data)
    channel_number <- ncol(in_data) + 1
    channel_id <- paste("$P", channel_number, sep = "")
    channel_name <- "cluster"
    channel_range <- max(cluster_data) + 1
    plist <- matrix(c(channel_name, channel_name, channel_range,
    0, channel_range - 1))
    rownames(plist) <- c("name", "desc", "range", "minRange",
    "maxRange")
    colnames(plist) <- c(channel_id)
    pd <- rbind(pd, t(plist))
    pData(params2) <- pd
    out_data <- cbind(in_data, cluster = cluster_data)
    out_frame <- flowFrame(out_data, params2, description = desc)
   keyval <- list()
   keyval[[paste("$P", channel_number, "B", sep = "")]] <- "32"
   keyval[[paste("$P", channel_number, "R", sep = "")]] <- toString(channel_range)
   keyval[[paste("$P", channel_number, "E", sep = "")]] <- "0,0"
   keyval[[paste("$P", channel_number, "N", sep = "")]] <- channel_name
   keyval[[paste("$P", channel_number, "S", sep = "")]] <- channel_name
   keyword(out_frame) <- keyval
   save(out_data,file="out_data.rdata")
   write.FCS(out_frame, newframeloc)
    # run original SPADE if specified
    
    if(run_spade) {
        # relevant files are copied into a new directory and run from there
        #15. If "run_spade" is true, copy all relevant files (FCS files, graph files, etc.) to a new folder and run original SPADE on these files. Note this is performed on data BEFORE it was reduced to clusters
        copyfiles <- c(density_files, sampled_files, graph_file, cells_file, clust_file)
        
        spadeDir <- "orig_spade/"
        if (! file.exists(paste(out_dir, spadeDir, sep=""))) {
            dir.create(paste(out_dir, spadeDir, sep=""))
        }
        
        new_out_dir <- paste(out_dir, spadeDir, sep="")
        file.copy(copyfiles, new_out_dir, overwrite=TRUE)
        
        new_files<-paste(new_out_dir, basename(files), sep="")

            new_density_files<-paste(new_out_dir, basename(density_files), sep="")
            print(new_density_files)
            new_sampled_files<-paste(new_out_dir, basename(sampled_files), sep="")
            print(new_sampled_files)
            new_graph_file<-paste(new_out_dir, basename(graph_file), sep="")
             print(new_graph_file)
            new_cells_file<-paste(new_out_dir, basename(cells_file), sep="")
            print( new_cells_file)
            new_clust_file<-paste(new_out_dir, basename(clust_file), sep="")
            print( new_clust_file)
        
        if( file.exists(graph_file)) {
            file.remove(graph_file)
        }
        
        
        r1 <- original_spade(new_files, file_pattern, new_out_dir, cluster_cols, panels, comp, arcsinh_cofactor, transforms, downsampling_target_number, downsampling_target_pctile, downsampling_target_percent, downsampling_exclude_pctile, k, clustering_samples, layout, pctile_color, fcs_channel_mappings_json, groups, new_cells_file, new_graph_file, new_clust_file, new_density_files, new_sampled_files)
        
        ##### Run spade forest ####
        if(run_spade_only) {
            stop("finished because only SPADE was run")
        } else {
    
    
    qcol = q[,cluster_cols]
    
    write.table(q, file=output1, row.names=FALSE, sep=",")
    write.table(qcol, file=output2, row.names=FALSE, sep=",")
    
    ###print(q)
    
    # relevant files are copied into a spadeforest directory and run from there
    #16. If "not run_spade_only" is true, copy all relevant files (FCS files, graph files, etc.) to a new spade forest folder and run spade forest analysis. Note this is performed on data BEFORE it was reduced to clusters
    copyfiles2 <- c(new_density_files, new_sampled_files, new_graph_file, new_cells_file, new_clust_file)
    
    spadeforestDir <- "orig_spade/spadeforest/"
    if (! file.exists(paste(out_dir, spadeforestDir, sep=""))) {
        dir.create(paste(out_dir, spadeforestDir, sep=""))
    }
    
    new_out_dir2 <- paste(out_dir, spadeforestDir, sep="")
    file.copy(copyfiles2, new_out_dir2, overwrite=TRUE)
    
    new_files2<-paste(new_out_dir2, basename(files), sep="")
    
    new_density_files2<-paste(new_out_dir2, basename(density_files), sep="")
    print(new_density_files2)
    new_sampled_files2<-paste(new_out_dir2, basename(sampled_files), sep="")
    print(new_sampled_files2)
    new_graph_file2<-paste(new_out_dir2, basename(graph_file), sep="")
    print(new_graph_file2)
    new_cells_file2<-paste(new_out_dir2, basename(cells_file), sep="")
    print( new_cells_file2)
    new_clust_file2<-paste(new_out_dir2, basename(clust_file), sep="")
    print( new_clust_file2)
    
    if( file.exists(new_graph_file2)) {
        file.remove(new_graph_file2)
    }
    
    graph_file = new_graph_file2
    print(graph_file)
    
    print("Now performing TAHC")
    
    #16. Run TAHC: Hierarchical clusering for spade forest See. published paper for reference
    
    t = dotahc(filein=new_out_dir, cols = cluster_cols, graph_cols = graph_cols, dirout=new_out_dir2, gated=gated,leafsort=leafsort)
    save(t,file=paste(new_out_dir2,"dendog_t.rdata",sep=""))
    # return(t)
    
    # separate regions of the graph into individual parts
    
    #17. separate regions. Cut the graph into dcut number of trees
   
   if (!is.null(reg)) {
       reg=reg
      }else{
       reg <- separate.regions(t,dcut,minsize)
   }
   
   ####
    print("Now generating spade forest")
    save(reg,file=paste(new_out_dir2,"reg.rdata",sep=""))
    print("Creating separate trees from original MST")
    glayout=makegraph3(input_graph_file=new_graph_file,output_graph_file=graph_file,reg=reg,centers,forestcluster,treedat=out_data,forestlayout=forestlayout)
    refgraph= igraph:::read.graph(new_graph_file,format="gml")
    
    
    #18. Creating pdfs for each marker of separate trees from original MST representing disjoint tree motifs
    r2 <- spade_forest(new_files2, file_pattern, new_out_dir2, cluster_cols, panels, comp, arcsinh_cofactor, transforms, downsampling_target_number, downsampling_target_pctile, downsampling_target_percent, downsampling_exclude_pctile, k, clustering_samples, layout=forestlayout, pctile_color, fcs_channel_mappings_json, groups, new_cells_file2, new_graph_file2, new_clust_file2, new_density_files2, new_sampled_files2,refgraph)
    
    
         }
    
    }

}
    
########## END SPADEDRIVERC #############
    ### Other codes #####

SPADE.plot.trees0 <-
function (graph, files, file_pattern = "*anno.Rsave", out_dir = ".",
layout = SPADE.layout.arch, attr_pattern = "percent|medians|fold|cvs",
scale = NULL, pctile_color = c(0.0, 1), normalize = "global",
size_scale_factor = 1, edge.color = "grey", bare = FALSE,
palette = "bluered")
{
    if (!is.igraph(graph)) {
        stop("Not a graph object")
    }
    if (!is.null(scale) && (!is.vector(scale) || length(scale) !=
    2)) {
        stop("scale must be a two element vector")
    }
    if (!is.vector(pctile_color) || length(pctile_color) != 2) {
        stop("pctile_color must be a two element vector with values in [0,1]")
    }
    if (length(files) == 1 && file.info(SPADE.strip.sep(files))$isdir) {
        files <- dir(SPADE.strip.sep(files), full.names = TRUE,
        pattern = glob2rx(file_pattern))
    }
    out_dir <- SPADE.normalize.out_dir(out_dir)
    load_attr <- function(save_file) {
        anno <- NULL
        l <- load(save_file)
        stopifnot(l == "anno")
        return(anno)
    }
    boundaries <- NULL
    if (normalize == "global") {
        boundaries <- c()
        all_attrs <- c()
        for (f in files) {
            attrs <- load_attr(f)
            for (i in grep(attr_pattern, colnames(attrs))) {
                n <- colnames(attrs)[i]
                all_attrs[[n]] <- c(all_attrs[[n]], attrs[, i])
            }
        }
        for (i in seq_along(all_attrs)) {
            boundaries[[names(all_attrs)[i]]] <- quantile(all_attrs[[i]],
            probs = pctile_color, na.rm = TRUE)
        }
    }
    if (is.function(layout))
    graph_l <- layout(graph)
    else graph_l <- layout
    if (palette == "jet")
    palette <- colorRampPalette(c("#00007F", "blue", "#007FFF",
    "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    else if (palette == "bluered")
    palette <- colorRampPalette(c("blue", "#007FFF", "cyan",
    "#7FFF7F", "yellow", "#FF7F00", "red"))
    else stop("Please use a supported color palette.  Options are 'bluered' or 'jet'")
    colorscale <- palette(100)
    for (f in files) {
        attrs <- load_attr(f)
        vsize <- attrs$percenttotal
        vsize <- vsize/(max(vsize, na.rm = TRUE)^(1/size_scale_factor)) *
        3 + 2
        vsize[is.na(vsize) | (attrs$count == 0)] <- 1
        for (i in grep(attr_pattern, colnames(attrs))) {
            name <- colnames(attrs)[i]
            attr <- attrs[, i]
            if (!is.null(scale))
            boundary <- scale
            else if (normalize == "global") {
                boundary <- boundaries[[name]]
            }
            else boundary <- quantile(attr, probs = pctile_color,
            na.rm = TRUE)
            if (length(grep("^medians|percent|cvs", name)))
            boundary <- c(min(boundary), max(boundary))
            else boundary <- c(-max(abs(boundary)), max(abs(boundary)))
            boundary <- round(boundary, 2)
            if (boundary[1] == boundary[2]) {
                boundary <- c(boundary[1] - 1, boundary[2] +
                1)
            }
            grad <- seq(boundary[1], boundary[2], length.out = length(colorscale))
            color <- colorscale[findInterval(attr, grad, all.inside = TRUE)]
            color[is.na(attr) | (attrs$count == 0)] <- "grey"
            if (grepl("^percenttotalratiolog$", name)) {
                color[is.na(attr) & attrs$count > 0] <- tail(colorscale,
                1)
            }
            fill_color <- color
            is.na(fill_color) <- is.na(attr)
            frame_color <- color
            pdf(paste(out_dir, basename(f), ".", name, ".pdf",
            sep = ""))
            graph_aspect <- ((max(graph_l[, 2]) - min(graph_l[,
            2]))/(max(graph_l[, 1]) - min(graph_l[, 1])))
            plot(graph, layout = graph_l, vertex.shape = "circle",
            vertex.color = fill_color, vertex.frame.color = frame_color,
            edge.color = edge.color, vertex.size = vsize,
            vertex.label = NA, edge.arrow.size = 0.25, edge.arrow.width = 1,
            asp = graph_aspect)
            if (!bare) {
                if (length(grep("^medians", name)))
                name <- sub("medians", "Median of ", name)
                else if (length(grep("^fold", name)))
                name <- sub("fold", "Arcsinh diff. of ", name)
                else if (grepl("^percenttotal$", name))
                name <- sub("percent", "Percent freq. of ",
                name)
                else if (grepl("^percenttotalratiolog$", name))
                name <- "Log10 of Ratio of Percent Total of Cells in Each Cluster"
                else if (grepl("^cvs", name))
                name <- sub("cvs", "Coeff. of Variation of ",
                name)
                if (grepl("_clust$", name))
                name <- sub("_clust", "\n(Used for tree-building)",
                name)
                title(main = paste(strsplit(basename(f), ".fcs")[[1]][1],
                sub = name, sep = "\n"))
                subplot(image(grad, c(1), matrix(1:length(colorscale),
                ncol = 1), col = colorscale, xlab = ifelse(is.null(scale),
                paste("Range:", pctile_color[1], "to", pctile_color[2],
                "pctile"), ""), ylab = "", yaxt = "n", xaxp = c(boundary,
                1)), x = "right,bottom", size = c(1, 0.2))
            }
            dev.off()
        }
    }
}

#############################
    SPADE.plot.trees2<-
    function (graph, files, file_pattern = "*anno.Rsave", out_dir = ".",
    layout = igraph:::layout.kamada.kawai, attr_pattern = "medians",
    scale = NULL, pctile_color = c(0, 1), normalize = "local",
    size_scale_factor = 1, edge.color = "grey", bare = FALSE,
    palette = "bluered",refgraph)
    {
        if (!is.igraph(graph)) {
            stop("Not a graph object")
        }
        if (!is.null(scale) && (!is.vector(scale) || length(scale) !=
        2)) {
            stop("scale must be a two element vector")
        }
        if (!is.vector(pctile_color) || length(pctile_color) != 2) {
            stop("pctile_color must be a two element vector with values in [0,1]")
        }
        if (length(files) == 1 && file.info(spade:::SPADE.strip.sep(files))$isdir) {
            files <- dir(spade:::SPADE.strip.sep(files), full.names = TRUE,
            pattern = glob2rx(file_pattern))
        }
        out_dir <- spade:::SPADE.normalize.out_dir(out_dir)
        load_attr <- function(save_file) {
            anno <- NULL
            l <- load(save_file)
            stopifnot(l == "anno")
            return(anno)
        }
        boundaries <- NULL
        if (normalize == "global") {
            boundaries <- c()
            all_attrs <- c()
            for (f in files) {
                attrs <- load_attr(f)
                for (i in grep(attr_pattern, colnames(attrs))) {
                    n <- colnames(attrs)[i]
                    all_attrs[[n]] <- c(all_attrs[[n]], attrs[, i])
                }
            }
            for (i in seq_along(all_attrs)) {
                boundaries[[names(all_attrs)[i]]] <- quantile(all_attrs[[i]],
                probs = pctile_color, na.rm = TRUE)
            }
        }
        if (is.function(layout)) {
        graph_l <- layout(graph)
        graph_r <- layout(refgraph)
            } else {
            graph_l <- layout
            graph_r <- layout
            }
        idref=match(V(graph)$name,V( refgraph)$name)
        if (palette == "jet")
        palette <- colorRampPalette(c("#00007F", "blue", "#007FFF",
        "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
        else if (palette == "bluered")
        palette <- colorRampPalette(c("blue", "#007FFF", "cyan",
        "#7FFF7F", "yellow", "#FF7F00", "red"))
        else stop("Please use a supported color palette.  Options are 'bluered' or 'jet'")
        colorscale <- palette(100)
        for (f in files) {
            attrs <- load_attr(f)
            vsize <- attrs$percenttotal
            vsize <- vsize/(max(vsize, na.rm = TRUE)^(1/size_scale_factor)) *
            3 + 2
            vsize[is.na(vsize) | (attrs$count == 0)] <- 1
            for (i in grep(attr_pattern, colnames(attrs))) {
                name <- colnames(attrs)[i]
                attr <- attrs[, i]
                if (!is.null(scale))
                boundary <- scale
                else if (normalize == "global") {
                    boundary <- boundaries[[name]]
                }
                else boundary <- quantile(attr, probs = pctile_color,
                na.rm = TRUE)
                if (length(grep("^medians|percent|cvs", name)))
                boundary <- c(min(boundary), max(boundary))
                else boundary <- c(-max(abs(boundary)), max(abs(boundary)))
                boundary <- round(boundary, 2)
                #print(paste("boundary=",boundary))
                if (boundary[1] == boundary[2]) {
                    boundary <- c(boundary[1] - 1, boundary[2] +
                    1)
                }
                grad <- seq(boundary[1], boundary[2], length.out = length(colorscale))
                color <- colorscale[findInterval(attr, grad, all.inside = TRUE)]
                color[is.na(attr) | (attrs$count == 0)] <- "grey"
                if (grepl("^percenttotalratiolog$", name)) {
                    color[is.na(attr) & attrs$count > 0] <- tail(colorscale,
                    1)
                }
                fill_color <- color[idref]
                is.na(fill_color) <- is.na(attr)
                frame_color <- color[idref]
                pdf(paste(out_dir, basename(f), ".", name, ".pdf",
                sep = ""))
                graph_aspect <- ((max(graph_l[, 2]) - min(graph_l[,
                2]))/(max(graph_l[, 1]) - min(graph_l[, 1])))
                plot(graph, layout = graph_l, vertex.shape = "circle",
                vertex.color = fill_color, vertex.frame.color = frame_color,
                edge.color = edge.color, vertex.size = vsize,
                vertex.label = NA, edge.arrow.size = 0.25, edge.arrow.width = 1,
                asp = graph_aspect)
                if (!bare) {
                    if (length(grep("^medians", name)))
                    name <- sub("medians", "Median of ", name)
                    else if (length(grep("^fold", name)))
                    name <- sub("fold", "Arcsinh diff. of ", name)
                    else if (grepl("^percenttotal$", name))
                    name <- sub("percent", "Percent freq. of ",
                    name)
                    else if (grepl("^percenttotalratiolog$", name))
                    name <- "Log10 of Ratio of Percent Total of Cells in Each Cluster"
                    else if (grepl("^cvs", name))
                    name <- sub("cvs", "Coeff. of Variation of ",
                    name)
                    if (grepl("_clust$", name))
                    name <- sub("_clust", "\n(Used for tree-building)",
                    name)
                    title(main = paste(strsplit(basename(f), ".fcs")[[1]][1],
                    sub = name, sep = "\n"))
                    spade:::subplot(image(grad, c(1), matrix(1:length(colorscale),
                    ncol = 1), col = colorscale, xlab = ifelse(is.null(scale),
                    paste("Range:", pctile_color[1], "to", pctile_color[2],
                    "pctile"), ""), ylab = "", yaxt = "n", xaxp = c(boundary,
                    1)), x = "right,bottom", size = c(1, 0.2))
                }
                dev.off()
            }
        }
    }
 
    subplot<-
    function (fun, x, y = NULL, size = c(1, 1), vadj = 0.5, hadj = 0.5,
    inset = c(0, 0), type = c("plt", "fig"), pars = NULL)
    {
        old.par <- par(no.readonly = TRUE)
        on.exit(par(old.par))
        type <- match.arg(type)
        if (missing(x))
        x <- locator(2)
        if (is.character(x)) {
            if (length(inset) == 1)
            inset <- rep(inset, 2)
            x.char <- x
            tmp <- par("usr")
            x <- (tmp[1] + tmp[2])/2
            y <- (tmp[3] + tmp[4])/2
            if (length(grep("left", x.char, ignore.case = TRUE))) {
                x <- tmp[1] + inset[1] * (tmp[2] - tmp[1])
                if (missing(hadj))
                hadj <- 0
            }
            if (length(grep("right", x.char, ignore.case = TRUE))) {
                x <- tmp[2] - inset[1] * (tmp[2] - tmp[1])
                if (missing(hadj))
                hadj <- 1
            }
            if (length(grep("top", x.char, ignore.case = TRUE))) {
                y <- tmp[4] - inset[2] * (tmp[4] - tmp[3])
                if (missing(vadj))
                vadj <- 1
            }
            if (length(grep("bottom", x.char, ignore.case = TRUE))) {
                y <- tmp[3] + inset[2] * (tmp[4] - tmp[3])
                if (missing(vadj))
                vadj <- 0
            }
        }
        xy <- xy.coords(x, y)
        if (length(xy$x) != 2) {
            pin <- par("pin")
            tmp <- cnvrt.coords(xy$x[1], xy$y[1], "usr")$plt
            x <- c(tmp$x - hadj * size[1]/pin[1], tmp$x + (1 - hadj) *
            size[1]/pin[1])
            y <- c(tmp$y - vadj * size[2]/pin[2], tmp$y + (1 - vadj) *
            size[2]/pin[2])
            xy <- cnvrt.coords(x, y, "plt")$fig
        }
        else {
            xy <- cnvrt.coords(xy, , "usr")$fig
        }
        par(pars)
        if (type == "fig") {
            par(fig = c(xy$x, xy$y), new = TRUE)
        }
        else {
            par(plt = c(xy$x, xy$y), new = TRUE)
        }
        fun
        tmp.par <- par(no.readonly = TRUE)
        return(invisible(tmp.par))
    }
    
dotrees<-
function(graph_file, output_dir, treeoutput=NULL, scale=NULL) {
    
    if(is.null(treeoutput)) {
        treeoutput<-""
    }
    
    
    # read the graph, generate a layout, and plot the tree, using a global scale
    
    print("drawing tree")
    mst_graph = igraph:::read.graph(graph_file,format="gml")
    layout = read.table(paste(output_dir,"layout.table",sep=.Platform$file.sep))
    
    print("scale is ")
    print(scale)
    t = SPADE.plot.trees(graph = mst_graph, files=output_dir, out_dir=paste(output_dir, paste(treeoutput, "pdf", sep=""), sep=.Platform$file.sep),layout=as.matrix(layout), attr_pattern = "medians",scale=scale,normalize = "local",pctile_color = c(0, 1))
    
    return(t)
}


dotrees2<-
function(graph_file, layout_dir,output_dir, treeoutput=NULL, scale=NULL,refgraph) {
    
    if(is.null(treeoutput)) {
        treeoutput<-""
    }
    
    
    # read the graph, generate a layout, and plot the tree, using a global scale
    
    print("drawing tree")
    mst_graph = igraph:::read.graph(graph_file,format="gml")
    ##layout = read.table(layout_dir)
    layout=igraph:::layout.kamada.kawai(mst_graph)*50
    print("scale is ")
    print(scale)
    t = SPADE.plot.trees2(graph = mst_graph, files=output_dir, out_dir=paste(output_dir,"pdf",sep=.Platform$file.sep), attr_pattern ="medians", layout=layout, scale=scale,normalize = "local",refgraph=refgraph)
    
    return(t)
}

makegraph<-
function(new.data, graph_file, run_silhouette=TRUE, transforms=NULL, groups2=NULL) {
    
    
    # after having generated coordinates for a graph, draw a minimum spanning tree between them and save it
        
        centers <- NULL
        
        is.na(groups2) <- which(groups2 == 0)
        for (i in c(1:max(groups2, na.rm = TRUE))) {
            obs <- which(groups2 == i)
            if (length(obs) > 1) {
                centers <- rbind(centers, colMeans(new.data[obs, , drop = FALSE]))
                # groups2[obs] <- nrow(centers)
            } else {
  	   		        is.na(groups2) <- obs
            }
            
        }
    
    
    spade:::SPADE.writeGraph(spade:::SPADE.clustersToMST(centers), graph_file)
    
}

makegraph2 <- function(input_graph_file,output_graph_file,reg,centers) {
    g = igraph:::read.graph(input_graph_file,format="gml")
    ##g <- igraph:::read.graph(paste(input_graph_file,"mst.gml",sep=.Platform$file.sep),format="gml")
    # after having generated coordinates for a subgraph, draw a minimum spanning tree between them and save it
    graphs=list()
    graphs2=list()
    V1=list()
    cmat=list()
    for (r in 1:length(reg)) {
        idx=reg[[r]]
        ###centers = (new.data[,c("X1", "X2")])
        graphs2[[r]] <- induced.subgraph(graph=g,vids=idx)
        V1[[r]]=as.numeric(as_ids(V(graphs2[[r]])))
         cmat[[r]]=centers[V1[[r]],,drop=FALSE]
        graphs[[r]] <- spade:::SPADE.clustersToMST(cmat[[r]])
        
    }
    save(graphs2,file="graphs2.rdata")
    save(graphs,file="graphs.rdata")
    save(V1,file="V1.rdata")
    save(cmat,file="cmat.rdata")
    ##layouts <- lapply(graphs, igraph:::layout.kamada.kawai)
    layouts <- lapply(graphs, igraph:::layout_in_circle)
    lay <- merge_coords(graphs, layouts)
    g1 <- disjoint_union(graphs)
    ##igraph:::write.graph(g1,file=input_graph_file,format="gml")
    #centers <-lay*50
   # centers <- centers
    ##spade:::SPADE.writeGraph(spade:::SPADE.clustersToMST(centers), output_graph_file)
    
    g2=mst(g1)
    spade:::SPADE.writeGraph(g2, output_graph_file)
    
    return(g2)
    
}

makegraph3 <- function(input_graph_file,output_graph_file,reg,centers,forestcluster,treedat,forestlayout) {
    g = igraph:::read.graph(input_graph_file,format="gml")
    if(!is.null(forestcluster)){
       #load("out_data.rdata")
       reg=list()
       for ( fr in 1:length(forestcluster)) {
       reg[[fr]]=treedat[,"cluster"][treedat[,"timeclust"]%in%forestcluster[[fr]]]
       
       }
       save(reg,file="reg.rdata")
       graphs=list()
       graphs2=list()
       V1=list()
       cmat=list()
       for (r in 1:length(reg)) {
           idx=reg[[r]]
           ###centers = (new.data[,c("X1", "X2")])
           graphs2[[r]] <- induced.subgraph(graph=g,vids=idx)
           V1[[r]]=as.numeric(as_ids(V(graphs2[[r]])))
            cmat[[r]]=centers[V1[[r]],,drop=FALSE]
           graphs[[r]] <- spade:::SPADE.clustersToMST(cmat[[r]])
      }
    } else {
    graphs=list()
    graphs2=list()
    V1=list()
    cmat=list()
    for (r in 1:length(reg)) {
        idx=reg[[r]]
        ###centers = (new.data[,c("X1", "X2")])
        graphs2[[r]] <- induced.subgraph(graph=g,vids=idx)
        V1[[r]]=as.numeric(as_ids(V(graphs2[[r]])))
         cmat[[r]]=centers[V1[[r]],,drop=FALSE]
        graphs[[r]] <- spade:::SPADE.clustersToMST(cmat[[r]])
        
      }
    }
    save(graphs2,file="graphs2.rdata")
    save(graphs,file="graphs.rdata")
    save(V1,file="V1.rdata")
    save(cmat,file="cmat.rdata")
    ##layouts <- lapply(graphs, igraph:::layout.kamada.kawai)
    layouts <- lapply(graphs,forestlayout )
    lay <- merge_coords(graphs, layouts)
    g1 <- disjoint_union(graphs)

    #g2=mst(g1)
    g2=g1
    spade:::SPADE.writeGraph(g2, output_graph_file)
    return(g2)

}

#### Comparing the clusters in the original weighted graph and the MST
## N is a confusion matrix with N_ij = # of nodes in the real cluster i that appears in the detected cluster j
# n is number of nodes in graph, C1=number of real clusters C2=number of detected clusters
nmi<-function(N,C1,C2) {
    n=sum(N)
    num1=0
    for ( i in 1:C1) {
        num2=0
        for (j in 1:C2) {
            if(N[i,j]==0) {
                num2=num2+0
                print(num2)
            }else{
            num2=num2+(N[i,j]*log(N[i,j]*n/(sum(N[i,])*sum(N[,j]))))
            print(num2)
            }
        }
        num1=num1+num2
    }
    num=-2*num1
    den1=0
    for ( i in 1:C1){
        if(sum(N[i,])==0) {
            den1=den1+0
            print(paste("den1",den1))
         }else {
        den1=den1+(sum(N[i,])*log(sum(N[i,])/n))
        print(paste("den1",den1))
         }
    }
    den2=0
    for ( j in 1:C2){
        if(sum(N[,j])==0) {
            den2=den2+0
            print(paste("den2",den2))
        }else {
        den2=den2+(sum(N[,j])*log(sum(N[,j])/n))
        print(paste("den2",den2))
        }
    }
    den=den1+den2
    nmiscore= num/den
    print(paste("num",num))
    print(paste("den",num))
    
  return(nmiscore)
}

###run s=nmi(N=A,C1=5,C2=4)
treegen2<-
function(files, file_pattern = "*.fcs", out_dir = ".", cluster_cols, panels = NULL, comp = FALSE, arcsinh_cofactor = NULL, transforms = flowCore::arcsinhTransform(a = 0, b = 0.2), downsampling_target_number = NULL, downsampling_target_pctile = 1, downsampling_target_percent = 1, downsampling_exclude_pctile = 0, k = k, clustering_samples = 50000, layout=igraph:::layout.kamada.kawai , pctile_color = c(0.0, 1), fcs_channel_mappings_json = NULL, graph_cols = NULL, graph_file, sampled_files, annoout, find_global_scale=FALSE) {
    
    # find_global_scale=TRUE indicates that the only goal is to return marker medians, not generate a tree. many processes will be skipped. Most of this is unchanged from regular SPADE
    
    if(!find_global_scale) {
        
        
        annoout <- spade:::SPADE.normalize.out_dir(annoout)
        
        graph <- read.graph(graph_file, format = "gml")
        
        num.clusters <- vcount(graph)
        layout_table <- layout(graph)
        if (deparse(substitute(layout)) != "SPADE.layout.arch")
        layout_table = layout_table * 50
        write.table(layout_table, paste(annoout, file = "layout.table",
        sep = "/"), row.names = FALSE, col.names = FALSE)
        attr_values <- list()
    }
    
    if (is.null(panels)) {
        
        panels <- list(list(panel_files = basename(files), median_cols = NULL))
    }
    print("111")
    for (p in panels) {
        reference_medians <- NULL
        if (!is.null(p$reference_files)) {
            reference_files <- sapply(as.vector(p$reference_files),
            function(f) {
                sampled_files[match(f, basename(files))[1]]
            })
            
            if(find_global_scale) {
                clust.num<-k
            } else {
                graph <- read.graph(graph_file, format = "gml")
                clust.num <- vcount(graph)
            }
            
            reference_medians <- spmm(reference_files,
            clust.num, cols = p$fold_cols, transforms = transforms,
            cluster_cols = cluster_cols, comp = comp)
            
            print("p$reference_files exist")
            
        }
        print("222")
        print(p)
        ## this part is the part that is run, this is important. The file going in has already been reduced to k cells. The FCS file, the number of nodes, a null vector called median_cols, and other less important stuff get passed in.
        #for (f in as.vector(p$panel_files)) {
        for (f in as.vector(p[[1]][1])) {
            # f <- sampled_files[match(f, basename(files))[1]]
            
            f <- paste(out_dir, f, sep="/")
            print(f)
            ## because there is only one cell in each cluster, the coefficient of variance for each cell is NA. It will be replaced with zero instead.
            
            message("Computing medians for file: ", f)
            
            if(find_global_scale) {
                
                # return(f)
                
                print("333")
                
                clust.num<-k
            } else {
                graph <- read.graph(graph_file, format = "gml")
                clust.num <- vcount(graph)
            }
            
            anno <- spmm(f, clust.num, cols = p$median_cols,
            transforms = transforms, cluster_cols = cluster_cols,
            comp = comp)
            print("444")
            if(find_global_scale) {
                attr_values<-NULL
                anno <- SPADE.flattenAnnotations(anno)
                for (c in colnames(anno)) {
                    attr_values[[c]] <- c(attr_values[[c]], anno[,c])
                }
                
                return(anno)
            }
            
            # anno contains the medians from here on
            
            print("555")
            if (!is.null(reference_medians)) { # i don't think this is performed, because reference_medians is null unless reference files were specified
                graph <- read.graph(graph_file, format = "gml")
                print("reference_medians are not null")
                message("Computing fold change for file: ", f)
                fold_anno <- spmm(f, vcount(graph),
                cols = p$fold_cols, transforms = transforms,
                cluster_cols = cluster_cols, comp = comp)
                fold <- fold_anno$medians - reference_medians$medians
                raw_fold <- fold_anno$raw_medians/reference_medians$raw_medians
                ratio <- log10(fold_anno$percenttotal/reference_medians$percenttotal)
                colnames(ratio) <- c("percenttotalratiolog")
                is.na(ratio) <- fold_anno$count == 0 | reference_medians$count ==
                0
                anno <- c(anno, list(percenttotalratiolog = ratio,
                fold = fold, raw_fold = raw_fold))
                
                
            }
            # this part creates a new .gml file which is the old one PLUS the medians. Probably want to make sure this is correct. The values written don't seem to correspond exactly?
            print("666")
            graph <- read.graph(graph_file, format = "gml")
            print(graph)
            SPADE.write.graph(SPADE.annotateGraph2(graph, layout = layout,
            anno = anno), paste(f, ".medians.gml", sep = ""),
            format = "gml")
            anno <- SPADE.flattenAnnotations(anno)
            for (c in colnames(anno)) {
                attr_values[[c]] <- c(attr_values[[c]], anno[,c])
            }
            print("777")
            
            # now attr_values is vector which contains the annotation data
            
            if(find_global_scale) {
                
                globout <- paste(f, "global.anno.Rsave", sep = ".")
                save(anno, file = globout)
                ###return(globout)
            }
            
            print("888")
            ##save(anno, file = paste(annoout, "1.anno.Rsave", sep = "/"))
            
            # save(anno, file = paste(f, "anno.Rsave", sep = "."))
        }
    }
    
    attr_ranges <- t(sapply(attr_values, function(x) {
        #quantile(x, probs = c(0, pctile_color, 1), na.rm = TRUE)
        quantile(x, probs = c(0, 1), na.rm = TRUE)
    }))
    rownames(attr_ranges) <- sapply(rownames(attr_ranges), function(x) {
        gsub("[^A-Za-z0-9_]", "", x)
    })
    print("999")
    write.table(attr_ranges, paste(out_dir, "global_boundaries.table",
    sep = "/"), col.names = FALSE)
    ruleDir = system.file(paste("tools", "PopulationRules", "",
    sep = .Platform$file.sep), package = "spade")
    print("10000")
    if (!is.null(fcs_channel_mappings_json)) {
        library("rjson")
        library("hash")
        fcs_channel_mapping = fromJSON(fcs_channel_mappings_json)
        cat("Evaluating Population Rules....\n")
        population_rule_mappings = SPADE.createPopulationMapping(sampled_files[1],
        ruleDir, fcs_channel_mapping)
        for (filename in names(population_rule_mappings)) {
            for (sampled_file in sampled_files) {
                cat(paste("Rule:", filename, "\n", sep = " "))
                message("Evaluating population rule: ", filename,
                " for file: ", sampled_file)
                SPADE.evaluateCellTypeRule(out_dir, sampled_file,
                ruleCols = population_rule_mappings[[filename]],
                ruleDir = ruleDir, ruleFile = filename)
            }
        }
    }
    print("11111111")
    message("Producing tables...")
    print(paste(annoout, "tables", sep = "/"))
    dir.create(paste(annoout, "tables", sep = "/"), recursive = TRUE, showWarnings=FALSE)
    #  dir.create(paste(out_dir, "tables", sep = "/"), recursive = TRUE,
    #      showWarnings = FALSE)
    
    files <- dir(annoout, full.names = TRUE, pattern = glob2rx("*.anno.Rsave"))
    
    # files <- dir(out_dir, full.names = TRUE, pattern = glob2rx("*.anno.Rsave"))
    params <- unique(as.vector(sapply(files, function(f) {
        load(f)
        colnames(anno)
    })))
    print("1")
    
    dir.create(paste(annoout, "tables", "byAttribute", sep = "/"), recursive = TRUE, showWarnings = FALSE)
    
    # dir.create(paste(out_dir, "tables", "byAttribute", sep = "/"),
    #     recursive = TRUE, showWarnings = FALSE)
    
    print("2")
    print(files)
    for (p in params) {
        pivot <- c()
        names <- c()
        for (f in files) {
            load(f)
            f = basename(f)
            if (p %in% colnames(anno)) {
                #print("3")
                pivot <- cbind(pivot, anno[, p])
                names <- c(names, f)
                ##print(pivot)
            }
        }
        names <- gsub("[[:alnum:][:punct:]]+/output/([[:alnum:][:punct:]]+).fcs.density.fcs.cluster.fcs.anno.Rsave",
        "\\1", names)
        #print(names)
        pivot <- cbind(1:nrow(pivot), pivot)
        colnames(pivot) <- c("name", names)
        #print("4")
        
        if (!is.null(pivot) && ncol(pivot) > 0) {
            
            write.csv(pivot, file = paste(annoout, "tables/byAttribute/",
            p, "_table", ".csv", sep = ""), row.names = FALSE)
            
            # write.csv(pivot, file = paste(out_dir, "tables/byAttribute/",
            # p, "_table", ".csv", sep = ""), row.names = FALSE)
        }
    }
    byNodeData = list()
    
    dir.create(paste(annoout, "tables", "bySample", sep = "/"),
    recursive = TRUE, showWarnings = FALSE)
    
    
    
    #    dir.create(paste(out_dir, "tables", "bySample", sep = "/"),
    #       recursive = TRUE, showWarnings = FALSE)
    for (f in files) {
        load(f)
        f = basename(f)
        pivot <- anno
        # print("5")
        names <- colnames(pivot)
        pivot <- cbind(1:nrow(pivot), pivot)
        colnames(pivot) <- c("ID", names)
        name <- gsub("output/([[:alnum:][:punct:]]+).fcs.density.fcs.cluster.fcs.anno.Rsave",
        "\\1", f)
        write.csv(pivot, file = paste(annoout, "tables/bySample/",
        name, "_table", ".csv", sep = ""), row.names = FALSE)
        
        #        write.csv(pivot, file = paste(out_dir, "tables/bySample/",
        #            name, "_table", ".csv", sep = ""), row.names = FALSE)
    }
    
    
    dir.create(paste(annoout, "tables", "byNodeID", sep = "/"),
    recursive = TRUE, showWarnings = FALSE)
    
    
    
    #   dir.create(paste(out_dir, "tables", "byNodeID", sep = "/"),
    #        recursive = TRUE, showWarnings = FALSE)
    print("6")
    print(head(pivot[,1:3]))
    for (node in rownames(pivot)) {
        tableData = list()
        for (f in files) {
            load(f)
            f = basename(f)
            tableData[[f]] = unlist(anno[node, , drop = T])
        }
        # print("7")
        tableData = do.call("cbind", tableData)
        write.csv(tableData, file = paste(annoout, "tables/byNodeID/",
        node, "_table", ".csv", sep = ""), row.names = TRUE,
        quote = F)
    }
    #print("8")
    ##return(anno)
    
}



forestgen<-
function(files, file_pattern = "*.fcs", out_dir = ".", cluster_cols, panels = NULL, comp = FALSE, arcsinh_cofactor = NULL, transforms = flowCore::arcsinhTransform(a = 0, b = 0.2), downsampling_target_number = NULL, downsampling_target_pctile = 1, downsampling_target_percent = 1, downsampling_exclude_pctile = 0, k = k, clustering_samples = 50000, layout=igraph:::layout.kamada.kawai , pctile_color = c(0.0, 1), fcs_channel_mappings_json = NULL, graph_cols = NULL, graph_file, sampled_files,find_global_scale=FALSE) {
    
    # find_global_scale=TRUE indicates that the only goal is to return marker medians, not generate a tree. many processes will be skipped. Most of this is unchanged from regular SPADE
    if(!find_global_scale) {
        
        
        annoout <- spade:::SPADE.normalize.out_dir(out_dir)
        
        graph <- read.graph(graph_file, format = "gml")
        
        num.clusters <- vcount(graph)
        layout_table <- layout(graph)
        if (deparse(substitute(layout)) != "SPADE.layout.arch")
        layout_table = layout_table * 50
        write.table(layout_table, paste(annoout, file = "layout.table",
        sep = "/"), row.names = FALSE, col.names = FALSE)
        attr_values <- list()
    }
    
    if (is.null(panels)) {
        
        panels <- list(list(panel_files = basename(files), median_cols = NULL))
    }
    print("111")
    for (p in panels) {
        reference_medians <- NULL
        if (!is.null(p$reference_files)) {
            reference_files <- sapply(as.vector(p$reference_files),
            function(f) {
                sampled_files[match(f, basename(files))[1]]
            })
            
            if(find_global_scale) {
                clust.num<-k
            } else {
                graph <- read.graph(graph_file, format = "gml")
                clust.num <- vcount(graph)
            }
            
            reference_medians <- spmm(reference_files,
            clust.num, cols = p$fold_cols, transforms = transforms,
            cluster_cols = cluster_cols, comp = comp)
            
            print("p$reference_files exist")
            
        }
        print("222")
        print(p)
        ## this part is the part that is run, this is important. The file going in has already been reduced to k cells. The FCS file, the number of nodes, a null vector called median_cols, and other less important stuff get passed in.
        #for (f in as.vector(p$panel_files)) {
        for (f in as.vector(p[[1]][1])) {
            # f <- sampled_files[match(f, basename(files))[1]]
            
            f <- paste(out_dir, f, sep="/")
            print(f)
            ## because there is only one cell in each cluster, the coefficient of variance for each cell is NA. It will be replaced with zero instead.
            
            message("Computing medians for file: ", f)
            
            if(find_global_scale) {
                
                # return(f)
                
                print("333")
                
                clust.num<-k
            } else {
                graph <- read.graph(graph_file, format = "gml")
                clust.num <- vcount(graph)
            }
            
            anno <- spmm(f, clust.num, cols = p$median_cols,
            transforms = transforms, cluster_cols = cluster_cols,
            comp = comp)
            print("444")
            if(find_global_scale) {
                attr_values<-NULL
                anno <- SPADE.flattenAnnotations(anno)
                for (c in colnames(anno)) {
                    attr_values[[c]] <- c(attr_values[[c]], anno[,c])
                }
                
                ###return(anno)
            }
            
            # anno contains the medians from here on
            
            print("555")
            if (!is.null(reference_medians)) { # i don't think this is performed, because reference_medians is null unless reference files were specified
                graph <- read.graph(graph_file, format = "gml")
                print("reference_medians are not null")
                message("Computing fold change for file: ", f)
                fold_anno <- spmm(f, vcount(graph),
                cols = p$fold_cols, transforms = transforms,
                cluster_cols = cluster_cols, comp = comp)
                fold <- fold_anno$medians - reference_medians$medians
                raw_fold <- fold_anno$raw_medians/reference_medians$raw_medians
                ratio <- log10(fold_anno$percenttotal/reference_medians$percenttotal)
                colnames(ratio) <- c("percenttotalratiolog")
                is.na(ratio) <- fold_anno$count == 0 | reference_medians$count ==
                0
                anno <- c(anno, list(percenttotalratiolog = ratio,
                fold = fold, raw_fold = raw_fold))
                
                
            }
            # this part creates a new .gml file which is the old one PLUS the medians. Probably want to make sure this is correct. The values written don't seem to correspond exactly?
            print("666")
            graph <- read.graph(graph_file, format = "gml")
            print(graph)
            SPADE.write.graph(SPADE.annotateGraph2(graph, layout = layout,
            anno = anno), paste(f, ".medians.gml", sep = ""),
            format = "gml")
            anno <- SPADE.flattenAnnotations(anno)
            for (c in colnames(anno)) {
                attr_values[[c]] <- c(attr_values[[c]], anno[,c])
            }
            print("777")
            
            # now attr_values is vector which contains the annotation data
            
            if(find_global_scale) {
                
                globout <- paste(f, "global.anno.Rsave", sep = ".")
                # save(anno, file = globout)
                ###return(globout)
            }
            
            print("888")
            save(anno, file = paste(annoout, "1.anno.Rsave", sep = "/"))
            
            # save(anno, file = paste(f, "anno.Rsave", sep = "."))
        }
    }
    
    attr_ranges <- t(sapply(attr_values, function(x) {
        #quantile(x, probs = c(0, pctile_color, 1), na.rm = TRUE)
        quantile(x, probs = c(0, 1), na.rm = TRUE)
    }))
    rownames(attr_ranges) <- sapply(rownames(attr_ranges), function(x) {
        gsub("[^A-Za-z0-9_]", "", x)
    })
    print("999")
    write.table(attr_ranges, paste(out_dir, "global_boundaries.table",
    sep = "/"), col.names = FALSE)
    ruleDir = system.file(paste("tools", "PopulationRules", "",
    sep = .Platform$file.sep), package = "spade")
    print("10000")
    if (!is.null(fcs_channel_mappings_json)) {
        library("rjson")
        library("hash")
        fcs_channel_mapping = fromJSON(fcs_channel_mappings_json)
        cat("Evaluating Population Rules....\n")
        population_rule_mappings = SPADE.createPopulationMapping(sampled_files[1],
        ruleDir, fcs_channel_mapping)
        for (filename in names(population_rule_mappings)) {
            for (sampled_file in sampled_files) {
                cat(paste("Rule:", filename, "\n", sep = " "))
                message("Evaluating population rule: ", filename,
                " for file: ", sampled_file)
                SPADE.evaluateCellTypeRule(out_dir, sampled_file,
                ruleCols = population_rule_mappings[[filename]],
                ruleDir = ruleDir, ruleFile = filename)
            }
        }
    }
    print("11111111")
    message("Producing tables...")
    print(paste(annoout, "tables", sep = "/"))
    dir.create(paste(annoout, "tables", sep = "/"), recursive = TRUE, showWarnings=FALSE)
    #  dir.create(paste(out_dir, "tables", sep = "/"), recursive = TRUE,
    #      showWarnings = FALSE)
    
    files <- dir(annoout, full.names = TRUE, pattern = glob2rx("*.anno.Rsave"))
    
    # files <- dir(out_dir, full.names = TRUE, pattern = glob2rx("*.anno.Rsave"))
    params <- unique(as.vector(sapply(files, function(f) {
        load(f)
        colnames(anno)
    })))
    print("1")
    
    dir.create(paste(annoout, "tables", "byAttribute", sep = "/"), recursive = TRUE, showWarnings = FALSE)
    
    # dir.create(paste(out_dir, "tables", "byAttribute", sep = "/"),
    #     recursive = TRUE, showWarnings = FALSE)
    
     print("2")
     print(files)
    for (p in params) {
        pivot <- c()
        names <- c()
        for (f in files) {
            load(f)
            f = basename(f)
            if (p %in% colnames(anno)) {
                #print("3")
                pivot <- cbind(pivot, anno[, p])
                names <- c(names, f)
                ##print(pivot)
            }
        }
        names <- gsub("[[:alnum:][:punct:]]+/output/([[:alnum:][:punct:]]+).fcs.density.fcs.cluster.fcs.anno.Rsave",
        "\\1", names)
        #print(names)
        pivot <- cbind(1:nrow(pivot), pivot)
        colnames(pivot) <- c("name", names)
        #print("4")
         
        if (!is.null(pivot) && ncol(pivot) > 0) {
            
            write.csv(pivot, file = paste(annoout, "tables/byAttribute/",
            p, "_table", ".csv", sep = ""), row.names = FALSE)
            
            # write.csv(pivot, file = paste(out_dir, "tables/byAttribute/",
            # p, "_table", ".csv", sep = ""), row.names = FALSE)
        }
    }
    byNodeData = list()
    
    dir.create(paste(annoout, "tables", "bySample", sep = "/"),
    recursive = TRUE, showWarnings = FALSE)
    
    
    
    #    dir.create(paste(out_dir, "tables", "bySample", sep = "/"),
    #       recursive = TRUE, showWarnings = FALSE)
    for (f in files) {
        load(f)
        f = basename(f)
        pivot <- anno
        # print("5")
        names <- colnames(pivot)
        pivot <- cbind(1:nrow(pivot), pivot)
        colnames(pivot) <- c("ID", names)
        name <- gsub("output/([[:alnum:][:punct:]]+).fcs.density.fcs.cluster.fcs.anno.Rsave",
        "\\1", f)
        write.csv(pivot, file = paste(annoout, "tables/bySample/",
        name, "_table", ".csv", sep = ""), row.names = FALSE)
        
        #        write.csv(pivot, file = paste(out_dir, "tables/bySample/",
        #            name, "_table", ".csv", sep = ""), row.names = FALSE)
    }
    
    
    dir.create(paste(annoout, "tables", "byNodeID", sep = "/"),
    recursive = TRUE, showWarnings = FALSE)
    
    
    
    #   dir.create(paste(out_dir, "tables", "byNodeID", sep = "/"),
    #        recursive = TRUE, showWarnings = FALSE)
    print("6")
    print(head(pivot[,1:3]))
    for (node in rownames(pivot)) {
        tableData = list()
        for (f in files) {
            load(f)
            f = basename(f)
            tableData[[f]] = unlist(anno[node, , drop = T])
        }
        # print("7")
        tableData = do.call("cbind", tableData)
        write.csv(tableData, file = paste(annoout, "tables/byNodeID/",
        node, "_table", ".csv", sep = ""), row.names = TRUE,
        quote = F)
    }
    #print("8")
    ###return(anno)
}



spmm<-
function (files, num.clusters, cols = NULL, arcsinh_cofactor = NULL,
transforms = flowCore::arcsinhTransform(a = 0, b = 0.2),
cluster_cols = NULL, comp = TRUE)
{
    if (!is.null(arcsinh_cofactor)) {
        warning("arcsinh_cofactor is deprecated, use transform=flowCore::arcsinhTransform(...) instead")
        transforms <- flowCore::arcsinhTransform(a = 0, b = 1/arcsinh_cofactor)
    }
    data <- c()
    
    
    
    
    files <- as.vector(files)
    for (f in files) {
        in_fcs <- spade:::SPADE.read.FCS(f, comp = comp)
        in_data <- exprs(in_fcs)
        params <- flowCore::parameters(in_fcs)
        pd <- pData(params)
        if (is.null(cols)) {
            cols <- as.vector(pd$name)
        }
        if (!"cluster" %in% cols) {
            cols <- c(cols, "cluster")
        }
        idxs <- match(cols, pd$name)
        if (any(is.na(idxs))) {
            stop("Invalid column specifier")
        }
        data <- rbind(data, in_data[, idxs, drop = FALSE])
    }
    clst <- 1:length(data[, "cluster"])
    data <- data[, colnames(data) != "cluster", drop = FALSE]
    data_t <-SPADE.transform.matrix(data, transforms)
    colnames(data) <- sapply(colnames(data), function(x) {
        if (x %in% cluster_cols)
        x <- paste(x, "clust", sep = "_")
        x
    })
    colnames(data_t) = colnames(data)
    ids <- 1:num.clusters
    
    
    if (any(is.na(match(unique(clst), ids)))) {
        stop("More clusters in FCS files than indicated")
    }
    count <- matrix(0, nrow = num.clusters, ncol = 1, dimnames = list(ids,
    "count"))
    medians <- matrix(NA, nrow = num.clusters, ncol = ncol(data),
    dimnames = list(ids, colnames(data)))
    raw_medians <- matrix(NA, nrow = num.clusters, ncol = ncol(data),
    dimnames = list(ids, colnames(data)))
    cvs <- matrix(NA, nrow = num.clusters, ncol = ncol(data),
    dimnames = list(ids, colnames(data)))
    for (i in ids) {
        
        ## check on data types(are the data subsets still matrices?)
        data_s <- subset(data, clst == i)
        data_s_t <- subset(data_t, clst == i)
        count[i, 1] <- nrow(data_s_t)
       # medians[i, ] <- apply(data_s_t, 2, function(x) median(x,na.rm = TRUE))
        #raw_medians[i, ] <- apply(data_s, 2, function(x) median(x,na.rm = TRUE))
        medians[i, ] <- apply(data_s_t, 2, function(x) medianna0(x))
        raw_medians[i, ] <- apply(data_s, 2, function(x) medianna0(x))
        
        #   cvs[i, ] <- apply(data_s_t, 2, function(d) {
        #       100 * sd(d)/abs(mean(d))
        #   })
        
        cvs[i, ] <- 0
        
    }
    
    # change the Percenttotal to be the percentage of total cells as determine dy the "cell_count" column
    
    # print(data)
    totalcells <- 0
    for(i in 1:nrow(data)) {
        totalcells <- (totalcells+data[i,"cell_count"])
        
    }
    cellcountmatrix <- NULL
    
    for(i in 1:nrow(data)) {
        cellcountmatrix<-c(cellcountmatrix, as.numeric(data[i,"cell_count"])/totalcells)
    }
    
    if(nrow(data)<100) {
        cellcountmatrix<-cellcountmatrix*1.5
    } else if(nrow(data)<50) {
        cellcountmatrix<-cellcountmatrix*2
    } else if(nrow(data)<25) {
        cellcountmatrix<-cellcountmatrix*3
    } else if(nrow(data)<12) {
        cellcountmatrix<-cellcountmatrix*4
    }
    
    percenttotal <- matrix(cellcountmatrix * 100, nrow = num.clusters,
    ncol = 1, dimnames = list(ids, "percenttotal"))
    list(count = count, medians = medians, raw_medians = raw_medians,
    cvs = cvs, percenttotal = percenttotal)
}


dotahc<-
function(filein,cols, graph_cols, dirout, gated,leafsort = TRUE) {
    ## oldwd = getwd()
    ### setwd(paste(oldwd,"spadeforest",sep="/"))
    
    ###embed <- paste(oldwd,"/ACCENSE_embed.", costtol, ".r", sep="")
    
    mst_graph <- igraph:::read.graph(paste(filein,"mst.gml",sep=.Platform$file.sep),format="gml")
    g=mst_graph
    C=distances(g)
    ###C=distances(g,algorithm ="unweighted")
    cluster=1:nrow(C)
    C1=cor(C, use="complete.obs", method="spearman" )
    Dsp=1-C1
    save(Dsp,file="Dsp.rdata")
    if(gated) {
    ###load("marrownodeid.rdata")
    nodeid=V(g)$name
     groupCodes <- nodeid
     rownames(Dsp) <- groupCodes
     colorCodes=c( "HSC"="#5E4FA2","MPP"="#4A67AD","NK"="#3780B9","Immature B"="#4199B5", "Pre-B I"="#58B2AB","Pre-B II"="#71C6A4", "Mature CD38lo B"="#8FD2A4"  , "Mature CD38mid B"="#ADDEA3","Plasma cell"="#C7E89E"  , "Naive CD4+ T"="#E0F299" ,"Naive CD8+ T"="#EEF8A5", "Mature CD4+ T"="#F9FCB6"    ,"Mature CD8+ T"="#FEF8B3"   ,  "Plasmacytoid DC"="#FEEA9D"  ,"CMP"="#FDDB87" , "Myelocyte"="#FDC575" ,"CD11bhi Monocyte"="#FDB062"  ,"CD11bmid Monocyte"="#F99455","CD11b- Monocyte"="#F57848", "Erythroblast"="#EB6046","GMP"="#DE4C4B",  "MEP"="#CD364D","Platelet"="#B51B47", "Megakaryocyte"="#9E0142")
     } else {
       groupCodes <- paste("cellid",1:nrow(Dsp))
        colorCodes <- jet.col (n = nrow(Dsp), alpha = 1)
     }
    
    h.avg  <- agnes(Dsp, method = "average")
    (d2 <- as.dendrogram(h.avg)) # two main branches
    labels_colors(d2) <- colorCodes[groupCodes][order.dendrogram(d2)]
    
    
    subDir <- "dendogramplots"
    
    if (! file.exists(paste(dirout, subDir, sep=""))){
        dir.create(paste(dirout, subDir, sep=""))
    }
    
    setwd(paste(dirout, subDir, sep=""))
    
    if(leafsort){
       t=dendsort(d2)
    } else {
        t=d2
    }
    ##### plot pdf ### d2[[2]][[1]][[2]]%>% labels
    pdf("MSTclustering.pdf",width=25, height=15)
    plot(t)
    plot(t[[1]][[1]])
    plot(t[[1]][[2]])
    plot(t[[2]][[1]])
    plot(t[[2]][[2]])
    #plot(d2[[1]][[1]][[1]])
    #plot(d2[[1]][[1]][[2]])
    #plot(d2[[1]][[2]][[1]])
    #plot(d2[[1]][[2]][[2]])
    #plot(d2[[2]][[1]][[1]])
    #plot(d2[[2]][[1]][[2]])
    #plot(d2[[2]][[2]][[1]])
    #plot(d2[[2]][[2]][[2]])
    dev.off();
    
 setwd("..")
 
 return(list(dendo=t,dist=Dsp))
 
}

dotahc2<-
function(filein,cols, graph_cols, dirout,dcut, gated) {
    require(sp)
    require(igraph)
    require(cluster)
    require(dendextend)
    require(ggplot2)
    require(dynamicTreeCut)
    require(colorspace)
    ## oldwd = getwd()
    ### setwd(paste(oldwd,"spadeforest",sep="/"))
    
    
    ###embed <- paste(oldwd,"/ACCENSE_embed.", costtol, ".r", sep="")
    
    mst_graph <- igraph:::read.graph(paste(filein,"mst.gml",sep=.Platform$file.sep),format="gml")
    g=mst_graph
    C=distances(g)
    ###C=distances(g,algorithm ="unweighted")
    cluster=1:nrow(C)
    C1=cor(C, use="complete.obs", method="spearman" )
    Dsp=1-C1
    ### if(gated) {
    ###load("marrownodeid.rdata")
    groupCodes <- nodeid
    rownames(Dsp) <- groupCodes
    colorCodes=c( "HSC"="#5E4FA2","MPP"="#4A67AD","NK"="#3780B9","Immature B"="#4199B5", "Pre-B I"="#58B2AB","Pre-B II"="#71C6A4", "Mature CD38lo B"="#8FD2A4"  , "Mature CD38mid B"="#ADDEA3","Plasma cell"="#C7E89E"  , "Naive CD4+ T"="#E0F299" ,"Naive CD8+ T"="#EEF8A5", "Mature CD4+ T"="#F9FCB6"    ,"Mature CD8+ T"="#FEF8B3"   ,  "Plasmacytoid DC"="#FEEA9D"  ,"CMP"="#FDDB87" , "Myelocyte"="#FDC575" ,"CD11bhi Monocyte"="#FDB062"  ,"CD11bmid Monocyte"="#F99455","CD11b- Monocyte"="#F57848", "Erythroblast"="#EB6046","GMP"="#DE4C4B",  "MEP"="#CD364D","Platelet"="#B51B47", "Megakaryocyte"="#9E0142")
    ##  } else {
    ###   groupCodes <- paste("celltype",1:nrow(Dsp))
    ###     colorCodes <- rep(c("red", "green", "blue", "yellow","black"),40)
    ### }
    
    h.avg  <- agnes(Dsp, method = "average")
    (d2 <- as.dendrogram(h.avg)) # two main branches
    labels_colors(d2) <- colorCodes[groupCodes][order.dendrogram(d2)]
    clusters <- cutreeDynamic(havg, distM = Dsp, method = "tree")
    clusters <- clusters[order.dendrogram(d2)]
    clusters_numbers <- unique(clusters) - (0 %in% clusters)
    n_clusters <- length(clusters_numbers)
    dend2 <- d2 %>%
    branches_attr_by_clusters(clusters, values = cols) %>%
    color_labels(col =   true_species_cols)
    
    subDir <- "dendogramplots"
    
    if (! file.exists(paste(dirout, subDir, sep=""))){
        dir.create(paste(dirout, subDir, sep=""))
    }
    
    setwd(paste(dirout, subDir, sep=""))
    
    ##### plot pdf ###
    pdf("MSTclustering.pdf",width=25, height=15)
    
    plot(d2)
    for (k in 1:dcut)
    plot(d2[[1]][[1]])
    plot(d2[[1]][[2]])
    plot(d2[[2]][[1]])
    plot(d2[[2]][[2]])
    
    dev.off();
    
    setwd("..")
    
    return(d2)
    
}

####Calinski-Harabasz Index and Boostrap Evaluation with Clustering Methods
# parameters : kmax = 10 ; clustermethod = "hclust"
# method = "ward.D" is one of parameters you can use to specify the
# the agglomeration method to be used for the hclust function in R

#Calinski-Harabasz Index : The math formula to the measure.
#  (SSB/SSW) * (N−k/k−1)  k=cluster N

#boot_clust <- ClusterBootstrap(data = mtcars_scaled,k = 4,bootstrap = 150,clustermethod = "kmeanspp",nstart = 10,iter.max = 100)

# print the returned list for clarity
#boot_clust
#criteria <- CHCriterion(data = mtcars_scaled, kmax = 10, clustermethod = "hclust", method = "ward.D")
# result
#head(criteria$data)
convertdata<-
function(infilename, k) {
    
    indata = spade:::SPADE.read.FCS(infilename)
    indata = exprs(indata)
    outdata = NULL
    lengthvec = NULL
    
    realgroups <- indata[,"cluster"]
    types <- indata[,"type"]
    
    
    indata = indata[,-which(colnames(indata)=="cluster")]
    coln = c(colnames(indata), "cell_count")
    
    
    # print("the following clusters were eliminated")
    for(i in 1:k) {
        rlist = which(realgroups==i)
        if(length(rlist) >= 1) {
            temp = indata[rlist,,drop=FALSE]
            
            lengthvec = c(lengthvec, length(rlist))
            
            if (length(rlist) > 1) {
                #outdata = rbind(outdata, apply(temp, 2, function(x) median(x,na.rm = TRUE)))
                outdata = rbind(outdata, apply(temp, 2, function(x) medianna0(x)))
                
            } else {
                outdata = rbind(outdata, temp)
                
            }
        } else {
            print(i)
        }
    }
    
    outdata = cbind(outdata, lengthvec)
    colnames(outdata) = coln
    ###print(head(outdata))
    return(outdata)
}


#####
##### SPADE.FCSToTree modified to return the cluster object that was generated

spadefcstotree2 <-
function (infilenames, outfilename, graphfilename, clusterfilename,
cols = NULL, k = 200, arcsinh_cofactor = NULL, transforms, desired_samples = NULL, comp = FALSE, remove_single_cell_clusters)
{
    if (!is.null(arcsinh_cofactor)) {
        warning("arcsinh_cofactor is deprecated, use transform=flowCore::arcsinhTransform(...) instead")
        transforms <- flowCore::arcsinhTransform(a = 0, b = 1/arcsinh_cofactor)
    }
    data = c()
    for (f in infilenames) {
        in_fcs <- spade:::SPADE.read.FCS(f, comp = comp)
        in_data <- exprs(in_fcs)
        params <- flowCore::parameters(in_fcs)
        pd <- pData(params)
        if (is.null(cols)) {
            cols <- as.vector(pd$name)
        }
        idxs <- match(cols, pd$name)
        if (any(is.na(idxs))) {
            stop("Invalid column specifier")
        }
        data <- rbind(data, in_data[, idxs, drop = FALSE])
        colnames(data) <- pd$name[idxs]
    }
    if(is.null(desired_samples)) {
        data <- data
    }else{
    if (nrow(data) > desired_samples) {
        data <- data[sample(1:nrow(data), desired_samples), ]
    } else if (nrow(data) == 0) {
        stop("Number of observations to cluster is zero. Did the density/downsampling warn about data similarity?")
    }
    }
    
    # return(list(data, transforms, k, remove_single_cell_clusters))
    print(" clustering code starts")
    
    clust <- spadecluster2(SPADE.transform.matrix(data, transforms),
    k, remove_single_cell_clusters=remove_single_cell_clusters)
    
    print(" clustering code ok")
    # re-assign k to the total number of clusters produced)
    k = max(clust$assign)
    
    mergeOrderPath = paste(dirname(outfilename), "/", "merge_order.txt",
    sep = "")
    write.table(clust$hclust$merge, file = mergeOrderPath, sep = "\t",
    quote = F, row.names = F, col.names = F)
    ff <- spade:::SPADE.build.flowFrame(subset(cbind(data, cluster = clust$assign),
    !is.na(clust$assign)))
    write.FCS(ff, outfilename)
    
    spade:::SPADE.writeGraph(spade:::SPADE.clustersToMST(clust$centers), graphfilename)
    write.table(clust$centers, file = clusterfilename, row.names = FALSE,
    col.names = colnames(data))
    
    
    return(list(clust, k))
    
}



spadecluster2<-
function (tbl, k, remove_single_cell_clusters)
{
    if (nrow(tbl) > 100000) {
        warning("Potentially too many observations for the clustering step",
        immediate = TRUE)
    }
    if (nrow(tbl) < k) {
        stop("Number of requested clusters exceeds number of events")
    }
    
    #d = dist(tbl)
    
    ###cluster = Rclusterpp.hclust(tbl)
    ##cluster = fastcluster:::hclust.vector(tbl, method='ward')
    cluster <- hclust2(objects=as.matrix(tbl), thresholdGini=0.2)
    #cluster <- hclust2(d, thresholdGini=0.2)
    # cluster = hclust(d)
    rawclust = cutree(cluster, k = k)
    
    clust = list(assign = rawclust)
    centers = c()
    is.na(clust$assign) <- which(clust$assign == 0)
    for (i in c(1:max(clust$assign, na.rm = TRUE))) {
        obs <- which(clust$assign == i)
        if (length(obs) > 1) {
            centers <- rbind(centers, colMeans(tbl[obs, , drop = FALSE]))
            clust$assign[obs] <- nrow(centers)
        } else {
            if(remove_single_cell_clusters) {
                is.na(clust$assign) <- obs
            } else {
                centers <- rbind(centers, tbl[obs, , drop = FALSE])
                clust$assign[obs] <- nrow(centers)
                
            }
        }
    }
    ###return(list(centers = centers, assign = clust$assign, hclust = cluster, distance = d))
    return(list(centers = centers, assign = clust$assign, hclust = cluster))
}


#####

##### from ccast_silhouette, should reassign cells to their nearest neighbor cluster if their silhouette width is <0. Input: distance matrix, outputfile for silhouette PDF, number of clusters (k)

ccast_silhouette2<-
function (d, fileout, nCluster.init = k, assign, iter.max = 5,
diagnosis = TRUE) {
    
    print("ccast_silhouette2 entered")
    
    groups = assign	 ### a tree output from cutree
    suppressMessages(require(ggplot2))
    sk0 = silhouette(groups, d)
    sk = sk0
    i = 1
    nc = vector()
    sw = vector()
    ssw = vector()
    while (any(sk[, 3] < 0) && i <= iter.max) {
        
        nc[i] = length(unique(sk[, "cluster"]))
        ssw[i] = sum(sk[, "sil_width"][sk[, "sil_width"] < 0])
        assign2 = ifelse(sk[, "sil_width"] < 0, sk[, "neighbor"],
        sk[, "cluster"])
        sk <- silhouette(assign2, d)
        i = i + 1
    }
    #  print(ssw)
    ssw.opt = max(ssw)
    p1 = ggplot(data.frame(numberOfClusters = nc, Iteration = seq_along(nc)),
    aes(x = Iteration, y = numberOfClusters)) + geom_point() +
    labs(x = "Iteration", y = "Number of clusters", title = "Profiling")
    dmat = data.frame(ssw = ssw, Iteration = seq_along(ssw))
    p2 = ggplot(dmat, aes(x = Iteration, y = ssw)) + geom_point() +
    geom_point(data = subset(dmat, ssw == ssw.opt), size = 5,
    colour = alpha("red", 1/2)) + labs(x = "Iteration",
    y = "Sum of negative Silhouette width", title = "Profiling")
    sk = sk0
    i = 1
    sw = -1
    ssw = -1e+06
    nc = vector()
    ssw.vec = vector()
    while (ssw < ssw.opt && i <= iter.max) {
        nc.opt = length(unique(sk[, "cluster"]))
        nc[i] = nc.opt
        ssw = sum(sk[, "sil_width"][sk[, "sil_width"] < 0])
        ssw.vec[i] = ssw
        assign2 = ifelse(sk[, "sil_width"] < 0, sk[, "neighbor"],
        sk[, "cluster"])
        sk <- silhouette(assign2, d)
        i = i + 1
    }
    p3 = ggplot(data.frame(numberOfClusters = nc, Iteration = seq_along(nc)),
    aes(x = Iteration, y = numberOfClusters)) + geom_point() +
    labs(x = "Iteration", y = "Number of clusters", title = "Iteration to optimal fit")
    dmat = data.frame(ssw = ssw.vec, Iteration = seq_along(ssw.vec))
    p4 = ggplot(dmat, aes(x = Iteration, y = ssw.opt)) + geom_point() +
    geom_point(data = subset(dmat, ssw == ssw.opt), size = 5,
    colour = alpha("red", 1/2)) + geom_text(data = subset(dmat,
    ssw == ssw.opt), aes(x = Iteration, y = ssw, label = paste("(",
    Iteration, ",", round(ssw.opt, 2), ")", sep = "")), hjust = 1,
    vjust = -1) + labs(x = "Iteration", y = "Sum of negative silhouette width",
    title = "Iterate to optimal fit") + ylim(c(min(ssw.vec),
    ifelse(max(ssw.vec) > 0, ssw.vec, 0)))
    if (diagnosis) {
        pdf(fileout)
        plot(sk0, main = "Silhouette plot: Original clusters")
        ###  print(p1)
        ###  print(p2)
        plot(sk, main = "Silhouette plot: Final clusters")
        # print(p3)
        dev.off()
    }
    gc()
    return(list(orgFit = sk0, optFit = sk, optFitCluster = data.frame(cluster = sk[,
    "cluster"], stringsAsFactors = F), optNumCluster = nc.opt,
    stopAtIteration = i, optSumOfNegSW = ssw.opt))
}


original_spade<-
function(files, file_pattern, out_dir, cluster_cols, panels, comp, arcsinh_cofactor, transforms, downsampling_target_number, downsampling_target_pctile, downsampling_target_percent, downsampling_exclude_pctile, k, clustering_samples, layout, pctile_color, fcs_channel_mappings_json, groups, cells_file, graph_file, clust_file, density_files, sampled_files) {
    
    
    
    sampled_files <- c()
    for (f in density_files) {
        message("Upsampling file: ", f)
        f_sampled <- paste(f, ".cluster.fcs", sep = "")
        ##spade:::SPADE.addClusterToFCS(f, f_sampled, cells_file, cols = cluster_cols,
        #transforms = transforms, comp = comp)
        SPADE.addClusterToFCS2(f, f_sampled, cells_file, cols = cluster_cols,
               transforms = transforms, comp = comp)
        sampled_files <- c(sampled_files, f_sampled)
    }
    # return(list(density_files, sampled_files))
    
    # return(sampled_files)
    graph <- read.graph(graph_file, format = "gml")
    #layout_table <- layout(graph)
    layout_table <- layout_with_kk(graph)
    if (deparse(substitute(layout)) != "SPADE.layout.arch")
    layout_table = layout_table * 50
    write.table(layout_table, paste(out_dir, file = "layout.table",
    sep = ""), row.names = FALSE, col.names = FALSE)
    attr_values <- list()
    if (is.null(panels)) {
        panels <- list(list(panel_files = basename(files), median_cols = NULL))
    }
    for (p in panels) {
        reference_medians <- NULL
        if (!is.null(p$reference_files)) {
            reference_files <- sapply(as.vector(p$reference_files),
            function(f) {
                sampled_files[match(f, basename(files))[1]]
            })
            reference_medians <- SPADE.markerMedians2(reference_files,
            vcount(graph), cols = p$fold_cols, transforms = transforms,
            cluster_cols = cluster_cols, comp = comp)
        }
        for (f in as.vector(p$panel_files)) {
            f <- sampled_files[match(f, basename(files))[1]]
            message("Computing medians for file: ", f)
            anno <- SPADE.markerMedians2(f, vcount(graph), cols = p$median_cols,
            transforms = transforms, cluster_cols = cluster_cols,
            comp = comp)
            
            print("1")
            if (!is.null(reference_medians)) {
                message("Computing fold change for file: ", f)
                fold_anno <- SPADE.markerMedians2(f, vcount(graph),
                cols = p$fold_cols, transforms = transforms,
                cluster_cols = cluster_cols, comp = comp)
                fold <- fold_anno$medians - reference_medians$medians
                raw_fold <- fold_anno$raw_medians/reference_medians$raw_medians
                ratio <- log10(fold_anno$percenttotal/reference_medians$percenttotal)
                colnames(ratio) <- c("percenttotalratiolog")
                is.na(ratio) <- fold_anno$count == 0 | reference_medians$count == 
                0
                anno <- c(anno, list(percenttotalratiolog = ratio, 
                fold = fold, raw_fold = raw_fold))
            }
            print("2")
            SPADE.write.graph(SPADE.annotateGraph(graph, layout = layout_table, 
            anno = anno), paste(f, ".medians.gml", sep = ""), 
            format = "gml")
            print("3")
            anno <- SPADE.flattenAnnotations(anno)
            for (c in colnames(anno)) {
                attr_values[[c]] <- c(attr_values[[c]], anno[, 
                c])
            }
            print("4")
            # save(anno, file = paste(annoout, "anno.Rsave", sep = "."))
            print("5")
            save(anno, file = paste(f, "anno.Rsave", sep = "."))
        }
    }
    attr_ranges <- t(sapply(attr_values, function(x) {
        #quantile(x, probs = c(0, pctile_color, 1), na.rm = TRUE)
        quantile(x, probs = c(0, 1), na.rm = TRUE)
    }))
    rownames(attr_ranges) <- sapply(rownames(attr_ranges), function(x) {
        gsub("[^A-Za-z0-9_]", "", x)
    })
    write.table(attr_ranges, paste(out_dir, "global_boundaries.table", 
    sep = ""), col.names = FALSE)
    ruleDir = system.file(paste("tools", "PopulationRules", "", 
    sep = .Platform$file.sep), package = "spade")
    if (!is.null(fcs_channel_mappings_json)) {
        library("rjson")
        library("hash")
        fcs_channel_mapping = fromJSON(fcs_channel_mappings_json)
        cat("Evaluating Population Rules....\n")
        population_rule_mappings = SPADE.createPopulationMapping(sampled_files[1], 
        ruleDir, fcs_channel_mapping)
        for (filename in names(population_rule_mappings)) {
            for (sampled_file in sampled_files) {
                cat(paste("Rule:", filename, "\n", sep = " "))
                message("Evaluating population rule: ", filename, 
                " for file: ", sampled_file)
                SPADE.evaluateCellTypeRule(out_dir, sampled_file, 
                ruleCols = population_rule_mappings[[filename]], 
                ruleDir = ruleDir, ruleFile = filename)
            }
        }
    }
    message("Producing tables...")
    dir.create(paste(out_dir, "tables", sep = "/"), recursive = TRUE, 
    showWarnings = FALSE)
    files <- dir(out_dir, full.names = TRUE, pattern = glob2rx("*.anno.Rsave"))
    params <- unique(as.vector(sapply(files, function(f) {
        load(f)
        colnames(anno)
    })))
    
    
    dir.create(paste(out_dir, "tables", "byAttribute", sep = "/"), 
    recursive = TRUE, showWarnings = FALSE)
    for (p in params) {
        pivot <- c()
        names <- c()
        for (f in files) {
            load(f)
            f = basename(f)
            if (p %in% colnames(anno)) {
                pivot <- cbind(pivot, anno[, p])
                names <- c(names, f)
            }
        }
        names <- gsub("[[:alnum:][:punct:]]+/output/([[:alnum:][:punct:]]+).fcs.density.fcs.cluster.fcs.anno.Rsave", 
        "\\1", names)
        pivot <- cbind(1:nrow(pivot), pivot)
        colnames(pivot) <- c("name", names)
        if (!is.null(pivot) && ncol(pivot) > 0) {
            write.csv(pivot, file = paste(out_dir, "tables/byAttribute/", 
            p, "_table", ".csv", sep = ""), row.names = FALSE)
        }
    }
    byNodeData = list()
    dir.create(paste(out_dir, "tables", "bySample", sep = "/"), 
    recursive = TRUE, showWarnings = FALSE)
    for (f in files) {
        load(f)
        f = basename(f)
        pivot <- anno
        names <- colnames(pivot)
        pivot <- cbind(1:nrow(pivot), pivot)
        colnames(pivot) <- c("ID", names)
        name <- gsub("output/([[:alnum:][:punct:]]+).fcs.density.fcs.cluster.fcs.anno.Rsave", 
        "\\1", f)
        write.csv(pivot, file = paste(out_dir, "tables/bySample/", 
        name, "_table", ".csv", sep = ""), row.names = FALSE)
    }
    dir.create(paste(out_dir, "tables", "byNodeID", sep = "/"), 
    recursive = TRUE, showWarnings = FALSE)
    for (node in rownames(pivot)) {
        tableData = list()
        for (f in files) {
            load(f)
            f = basename(f)
            tableData[[f]] = unlist(anno[node, , drop = T])
        }
        tableData = do.call("cbind", tableData)
        write.csv(tableData, file = paste(out_dir, "tables/byNodeID/", 
        node, "_table", ".csv", sep = ""), row.names = TRUE, 
        quote = F)
    }
    
    dotrees(graph_file, output_dir = out_dir)
    invisible(NULL)
    
    
    
}


spade_forest<-
function(files, file_pattern, out_dir, cluster_cols, panels, comp, arcsinh_cofactor, transforms, downsampling_target_number, downsampling_target_pctile, downsampling_target_percent, downsampling_exclude_pctile, k, clustering_samples, layout, pctile_color, fcs_channel_mappings_json, groups, cells_file, graph_file, clust_file, density_files, sampled_files,refgraph) {
    
    
    
    #sampled_files <- c()
    # for (f in density_files) {
    #   message("Upsampling file: ", f)
    #   f_sampled <- paste(f, ".cluster.fcs", sep = "")
    #   SPADE.addClusterToFCS(f, f_sampled, cells_file, cols = cluster_cols,
    #    transforms = transforms, comp = comp)
    #   sampled_files <- c(sampled_files, f_sampled)
    #}
    # return(list(density_files, sampled_files))
    
    # return(sampled_files)
    graph <- read.graph(graph_file, format = "gml")
    layout_table <- layout(graph)
    if (deparse(substitute(layout)) != "SPADE.layout.arch")
    layout_table = layout_table * 50
    write.table(layout_table, paste(out_dir, file = "layout.table",
    sep = ""), row.names = FALSE, col.names = FALSE)
    attr_values <- list()
    if (is.null(panels)) {
        panels <- list(list(panel_files = basename(files), median_cols = NULL))
    }
    for (p in panels) {
        reference_medians <- NULL
        if (!is.null(p$reference_files)) {
            reference_files <- sapply(as.vector(p$reference_files),
            function(f) {
                sampled_files[match(f, basename(files))[1]]
            })
            reference_medians <- SPADE.markerMedians2(reference_files,
            vcount(graph), cols = p$fold_cols, transforms = transforms,
            cluster_cols = cluster_cols, comp = comp)
        }
        for (f in as.vector(p$panel_files)) {
            f <- sampled_files[match(f, basename(files))[1]]
            message("Computing medians for file: ", f)
            anno <- SPADE.markerMedians2(f, vcount(graph), cols = p$median_cols,
            transforms = transforms, cluster_cols = cluster_cols,
            comp = comp)
            
            print("1")
            if (!is.null(reference_medians)) {
                message("Computing fold change for file: ", f)
                fold_anno <- SPADE.markerMedians2(f, vcount(graph),
                cols = p$fold_cols, transforms = transforms,
                cluster_cols = cluster_cols, comp = comp)
                fold <- fold_anno$medians - reference_medians$medians
                raw_fold <- fold_anno$raw_medians/reference_medians$raw_medians
                ratio <- log10(fold_anno$percenttotal/reference_medians$percenttotal)
                colnames(ratio) <- c("percenttotalratiolog")
                is.na(ratio) <- fold_anno$count == 0 | reference_medians$count ==
                0
                anno <- c(anno, list(percenttotalratiolog = ratio,
                fold = fold, raw_fold = raw_fold))
            }
            print("2")
            SPADE.write.graph(SPADE.annotateGraph2(graph, layout = layout_table,
            anno = anno), paste(f, ".medians.gml", sep = ""),
            format = "gml")
            print("3")
            anno <- SPADE.flattenAnnotations(anno)
            for (c in colnames(anno)) {
                attr_values[[c]] <- c(attr_values[[c]], anno[,
                c])
            }
            print("4")
            # save(anno, file = paste(annoout, "anno.Rsave", sep = "."))
            print("5")
            save(anno, file = paste(f, "anno.Rsave", sep = "."))
        }
    }
    attr_ranges <- t(sapply(attr_values, function(x) {
        ##quantile(x, probs = c(0, pctile_color, 1), na.rm = TRUE)
        quantile(x, probs = c(0, 1), na.rm = TRUE)
    }))
    rownames(attr_ranges) <- sapply(rownames(attr_ranges), function(x) {
        gsub("[^A-Za-z0-9_]", "", x)
    })
    write.table(attr_ranges, paste(out_dir, "global_boundaries.table",
    sep = ""), col.names = FALSE)
    ruleDir = system.file(paste("tools", "PopulationRules", "",
    sep = .Platform$file.sep), package = "spade")
    if (!is.null(fcs_channel_mappings_json)) {
        library("rjson")
        library("hash")
        fcs_channel_mapping = fromJSON(fcs_channel_mappings_json)
        cat("Evaluating Population Rules....\n")
        population_rule_mappings = SPADE.createPopulationMapping(sampled_files[1],
        ruleDir, fcs_channel_mapping)
        for (filename in names(population_rule_mappings)) {
            for (sampled_file in sampled_files) {
                cat(paste("Rule:", filename, "\n", sep = " "))
                message("Evaluating population rule: ", filename,
                " for file: ", sampled_file)
                SPADE.evaluateCellTypeRule(out_dir, sampled_file,
                ruleCols = population_rule_mappings[[filename]],
                ruleDir = ruleDir, ruleFile = filename)
            }
        }
    }
    message("Producing tables...")
    dir.create(paste(out_dir, "tables", sep = "/"), recursive = TRUE,
    showWarnings = FALSE)
    files <- dir(out_dir, full.names = TRUE, pattern = glob2rx("*.anno.Rsave"))
    params <- unique(as.vector(sapply(files, function(f) {
        load(f)
        colnames(anno)
    })))
    
    
    dir.create(paste(out_dir, "tables", "byAttribute", sep = "/"),
    recursive = TRUE, showWarnings = FALSE)
    for (p in params) {
        pivot <- c()
        names <- c()
        for (f in files) {
            load(f)
            f = basename(f)
            if (p %in% colnames(anno)) {
                pivot <- cbind(pivot, anno[, p])
                names <- c(names, f)
            }
        }
        names <- gsub("[[:alnum:][:punct:]]+/output/([[:alnum:][:punct:]]+).fcs.density.fcs.cluster.fcs.anno.Rsave",
        "\\1", names)
        pivot <- cbind(1:nrow(pivot), pivot)
        colnames(pivot) <- c("name", names)
        if (!is.null(pivot) && ncol(pivot) > 0) {
            write.csv(pivot, file = paste(out_dir, "tables/byAttribute/",
            p, "_table", ".csv", sep = ""), row.names = FALSE)
        }
    }
    byNodeData = list()
    dir.create(paste(out_dir, "tables", "bySample", sep = "/"),
    recursive = TRUE, showWarnings = FALSE)
    for (f in files) {
        load(f)
        f = basename(f)
        pivot <- anno
        names <- colnames(pivot)
        pivot <- cbind(1:nrow(pivot), pivot)
        colnames(pivot) <- c("ID", names)
        name <- gsub("output/([[:alnum:][:punct:]]+).fcs.density.fcs.cluster.fcs.anno.Rsave",
        "\\1", f)
        write.csv(pivot, file = paste(out_dir, "tables/bySample/",
        name, "_table", ".csv", sep = ""), row.names = FALSE)
    }
    dir.create(paste(out_dir, "tables", "byNodeID", sep = "/"),
    recursive = TRUE, showWarnings = FALSE)
    for (node in rownames(pivot)) {
        tableData = list()
        for (f in files) {
            load(f)
            f = basename(f)
            tableData[[f]] = unlist(anno[node, , drop = T])
        }
        tableData = do.call("cbind", tableData)
        write.csv(tableData, file = paste(out_dir, "tables/byNodeID/",
        node, "_table", ".csv", sep = ""), row.names = TRUE,
        quote = F)
    }
    
    dotrees2(graph_file, output_dir = out_dir,refgraph=refgraph)
    invisible(NULL)
    
    
}

separate.regions<-function(dendo,dcut,minsize) {
    require(dendextend)
    require(dynamicTreeCut)
    require(colorspace)
    # given data frame of 2D coordinates, find distinct regions.
    # calculate pairwise distance?
    # then set distance-cutoff which will determine if neighbor or not
    # there may be better alternative methods to determine if nodes are neighbors. probability based may work.
    # then add second neighbors to each node's list. 
    # Iterate a certain number of times until we are confident all neighbors have been added.
    # Then identify non-overlapping sets of neighborhoods 
    # calculate pairwise distance?
    ##n=dendo %>% nleaves
    print("generating separate clusters of cells from spade branches")
    if(!is.null(dcut)) {
     reg <- cutree(dendo[[1]], dcut)
     clusterlist<-list()
     for (i in 1:dcut){
      clusterlist[[i]] = c(which(reg==i))
     }
    } else {
    # Find special clusters:
    dist=dendo[[2]]
    h.avg  <- agnes(dist, method = "average")
    dendo1=as.hclust(h.avg)
    maxh=round(max(dendo1$height),0)
    #clusters <- cutreeDynamic(dendo1, method = "tree",minClusterSize = minsize)
    clusters <-cutreeDynamicTree(dendo1, maxTreeHeight = maxh, deepSplit = TRUE, minModuleSize = minsize)
    #cutreeDynamicTree(dendo1, distM = as.matrix(dist), method = "tree")
    #clusters <- cutreeDynamicTree(dendo1, maxTreeHeight = 1, deepSplit = TRUE, minModuleSize = 5)
   
    ###h.avg  <- agnes(Dsp, method = "average")
     if(0 %in% clusters) {
    clusters1 <- clusters[order.dendrogram(dendo[[1]])]+1
     }else{
    clusters1 <- clusters[order.dendrogram(dendo[[1]])]
     }

    #clusters <- clusters[dendo]
    #clusters_numbers <- unique(clusters1) - (0 %in% clusters1)
    clusters_numbers <- unique(clusters1)
    n_clusters <- length(clusters_numbers)
    colorCodes <- jet.col (n =  n_clusters , alpha = 1)
    labels_colors(dendo[[1]]) <- colorCodes[clusters1]
    cols=colorCodes[clusters1]
    
    pdf("cutreeDynamiclusters.pdf",width=25, height=15)
    #true_species_cols <- rainbow_hcl(max(clusters))[as.factor(clusters)]
    dend2 <- dendo[[1]] %>%
    branches_attr_by_clusters(clusters, values =  colorCodes[clusters_numbers]) %>%
    color_labels(col =  cols)
   
    plot(dend2)
    dev.off();
    
    clusterlist<-list()
    for (i in 1:n_clusters){
     clusterlist[[i]] = c(which(clusters1==i))
    }
    }
    return(clusterlist)
}

SPADE.annotateGraph2 <- function (graph, layout = NULL, anno = NULL)
{
    if (!is.igraph(graph)) {
        stop("Not a graph object")
    }
    if (!is.null(layout) && is.matrix(layout)) {
        if (nrow(layout) != vcount(graph) || ncol(layout) !=
        2) {
            stop("Ill-formated layout matrix, must 2 columns (x,y) and as many rows as vertices")
        }
        v_attr <- list.vertex.attributes(graph)
        for (i in grep("^graphics$", v_attr)) graph <- remove.vertex.attribute(graph,
        v_attr[i])
        graph <- set.vertex.attribute(graph, "graphics.x", value = layout[,
        1])
        graph <- set.vertex.attribute(graph, "graphics.y", value = layout[,
        2])
    }
    if (is.null(anno)) {
        return(graph)
    }
    else if (!is.list(anno)) {
        stop("anno must be a list with named entries")
    }
    for (i in seq_len(length(anno))) {
        l <- as.matrix(anno[[i]]) #### added as matrix
        if (!is.matrix(l)) {
            stop(paste("Argument:", quote(l), "must be a matrix"))
        }
        vt <- rownames(l)[rownames(l)%in%V(graph)$name]
        for (c in colnames(l)) {
            graph <- set.vertex.attribute(graph, ifelse(names(anno)[i] ==
            c, c, paste(names(anno)[i], c, sep = "")), index = vt,
            value = l[vt, c])
        }
    }
    graph
}

#### Beta incomplete version #####
matrixtofcs <-function(file1,file2,filename) {
    d1=read.FCS(file1,transformation=FALSE)
    
    in_data <- exprs(d1)
    cc=colnames(file2)
    d1@parameters@data$name[1:dim(d1)[2]]<-cc[1:dim(d1)[2]]
    params <- flowCore::parameters(d1)
    desc <- description(d1)
    pd <- pData(params)
    pd$desc<-d1@parameters@data$name
    ## in_dat=in_data[1:dim(file2)[1],1:dim(in_data)[2]]
   
    
    if (dim(file2)[2]<=dim(in_data)[2]) {
        in_dat=file2
        out_frame<- flowFrame(in_dat, params, description = desc)
         write.FCS(out_frame,filename)
    }else{
   
    in_dat=file2[,1:dim(in_data)[2]]
    counter=dim(in_dat)[2]
    
    for (k in (dim(in_dat)[2]+1):dim(file2)[2]) {
        ###typevec1 <- file2[,k]
        channel_number <- ncol(in_dat) + 1
        channel_id <- paste("$P", channel_number, sep = "")
        channel_name <- cc[k]
        channel_range <- 46
        plist <- matrix(c(channel_name, channel_name, channel_range,
        0, channel_range - 1))
        rownames(plist) <- c("name", "desc", "range", "minRange",
        "maxRange")
        colnames(plist) <- c(channel_id)
        pd <- rbind(pd, t(plist))
        pData(params) <- pd
        #channel_ids <- paste("$P", channel_number,"S", sep = "")
        #desc[[channel_ids]]<-channel_name
        
        ##out_data <- cbind(in_data[,i], noquote(paste(colnames(file2)[1],"=","typevec")))
      in_dat <- cbind(in_dat, file2[,k])
        colnames(in_dat)=cc[1:k]
        
        out_frame<- flowFrame(in_dat, params, description = desc)
        keyval <- list()
        keyval[[paste("$P", channel_number, "B", sep = "")]] <- "32"
        keyval[[paste("$P", channel_number, "R", sep = "")]] <- toString(channel_range)
        keyval[[paste("$P", channel_number, "E", sep = "")]] <- "0,0"
        keyval[[paste("$P", channel_number, "N", sep = "")]] <- channel_name
        keyval[[paste("$P", channel_number, "S", sep = "")]] <- channel_name
        keyword(out_frame) <- keyval
        counter=counter+1
        }
    
      print(paste("counter=",counter))
      write.FCS(out_frame,filename)
    }
}


matrixtofcs2 <-function(file1,file2,filename) {
    
    d1=read.FCS(file1,transformation=FALSE)
    
    Dall=as.matrix(file2)
    d1@parameters@data$name[1:dim(Dall)[2]]=colnames(Dall)
    
    param=flowCore::parameters(d1)
    out_frame <- flowFrame(Dall, param, description = description(d1))
    write.FCS(out_frame, filename)
    return(Dall)
    
}

SPADE.addClusterToFCS2<-function (infilename, outfilename, clusterfilename, cols = NULL,
    arcsinh_cofactor = NULL, transforms = flowCore::arcsinhTransform(a = 0,
        b = 0.2), comp = TRUE)
{
    if (!is.null(arcsinh_cofactor)) {
        warning("arcsinh_cofactor is deprecated, use transform=flowCore::arcsinhTransform(...) instead")
        transforms <- flowCore::arcsinhTransform(a = 0, b = 1/arcsinh_cofactor)
    }
    in_fcs <- spade:::SPADE.read.FCS(infilename, comp = comp)
    in_data <- exprs(in_fcs)
    params <- flowCore::parameters(in_fcs)
    pd <- pData(params)
    if (is.null(cols)) {
        cols <- as.vector(pd$name)
    }
    idxs <- match(cols, pd$name)
    if (any(is.na(idxs))) {
        stop("Invalid column specifier")
    }
    cluster_fcs <- spade:::SPADE.read.FCS(clusterfilename, comp = comp)
    cluster_data <- exprs(cluster_fcs)
    cluster_params <- flowCore::parameters(cluster_fcs)
    cluster_pd <- pData(cluster_params)
    c_idxs <- match(cols, cluster_pd$name)
    na.fail(c_idxs)
    assign<-cluster_data[, "cluster"]
    
    #assign <- SPADE.assignToCluster(SPADE.transform.matrix(in_data[,
   #     idxs], transforms), SPADE.transform.matrix(cluster_data[,
    #    c_idxs], transforms), cluster_data[, "cluster"])
    in_fcs <- spade:::SPADE.read.FCS(infilename, comp = FALSE, transform = FALSE)
    in_data <- exprs(in_fcs)
    channel_number <- ncol(in_fcs) + 1
    channel_id <- paste("$P", channel_number, sep = "")
    channel_name <- "cluster"
    channel_range <- max(assign) + 1
    plist <- matrix(c(channel_name, channel_name, channel_range,
        0, channel_range - 1))
    rownames(plist) <- c("name", "desc", "range", "minRange",
        "maxRange")
    colnames(plist) <- c(channel_id)
    pd <- rbind(pd, t(plist))
    pData(params) <- pd
    out_data <- cbind(in_data, cluster = assign)
    out_frame <- flowFrame(out_data, params, description = description(in_fcs))
    keyval <- list()
    keyval[[paste("$P", channel_number, "B", sep = "")]] <- "32"
    keyval[[paste("$P", channel_number, "R", sep = "")]] <- toString(channel_range)
    keyval[[paste("$P", channel_number, "E", sep = "")]] <- "0,0"
    keyval[[paste("$P", channel_number, "N", sep = "")]] <- channel_name
    keyval[[paste("$P", channel_number, "S", sep = "")]] <- channel_name
    keyword(out_frame) <- keyval
    write.FCS(out_frame, outfilename)
}
# Transpose data before call to put in row major order
SPADE.assignToCluster <- function(tbl, cluster_data, cluster_assign)
#require(spade)
   # .Call("SPADE_assign",t(tbl),t(cluster_data),as.integer(cluster_assign))

SPADE.addClusterToFCS <- function(
    infilename,
    outfilename,
    clusterfilename,
    cols=NULL,
    arcsinh_cofactor=NULL,
    transforms=flowCore::arcsinhTransform(a=0, b=0.2),
    comp=TRUE
) {


    if (!is.null(arcsinh_cofactor)) {
        warning("arcsinh_cofactor is deprecated, use transform=flowCore::arcsinhTransform(...) instead")
        transforms <- flowCore::arcsinhTransform(a=0, b=1/arcsinh_cofactor)
    }
    

    # Load in FCS file
    in_fcs  <- spade:::SPADE.read.FCS(infilename,comp=comp);
    in_data <- exprs(in_fcs);
    
    params <- flowCore::parameters(in_fcs);
    pd     <- pData(params);
    
    # Select out the desired columns
    if (is.null(cols)) {
        cols <- as.vector(pd$name)
    }
    idxs <- match(cols,pd$name)
    if (any(is.na(idxs))) {
        stop("Invalid column specifier")
    }
    
    # Load in clustered FCS file
    cluster_fcs  <- spade:::SPADE.read.FCS(clusterfilename, comp=comp)
    cluster_data <- exprs(cluster_fcs)
    
    cluster_params <- flowCore::parameters(cluster_fcs)
    cluster_pd     <- pData(cluster_params)
    
    c_idxs <- match(cols, cluster_pd$name)
    na.fail(c_idxs)
    
    # Assign observations to clusters
    assign <- SPADE.assignToCluster(
        SPADE.transform.matrix(in_data[,idxs], transforms),
        SPADE.transform.matrix(cluster_data[,c_idxs], transforms),
        cluster_data[,"cluster"]
    )
    
    # Reload FCS file without transformation, so it can be compactly rewritten...
    in_fcs <- spade:::SPADE.read.FCS(infilename,comp=FALSE,transform=FALSE)
    in_data <- exprs(in_fcs);
    
    # Add column named "cluster" to the FCS file
    channel_number <- ncol(in_fcs)+1;
    channel_id     <- paste("$P",channel_number,sep="");
    channel_name   <- "cluster";
    channel_range  <- max(assign)+1;
    
    plist <- matrix(c(channel_name,channel_name,channel_range,0,channel_range-1));
    rownames(plist) <- c("name","desc","range","minRange","maxRange");
    colnames(plist) <- c(channel_id);
    
    pd <- rbind(pd,t(plist));
    pData(params) <- pd;
    
    out_data <- cbind(in_data,"cluster"=assign);
    out_frame <- flowFrame(out_data,params,description=description(in_fcs));
    
    keyval <- list();
    keyval[[paste("$P",channel_number,"B",sep="")]] <- "32";            # Number of bits
    keyval[[paste("$P",channel_number,"R",sep="")]] <- toString(channel_range); # Range
    keyval[[paste("$P",channel_number,"E",sep="")]] <- "0,0";            # Exponent
    keyval[[paste("$P",channel_number,"N",sep="")]] <- channel_name;        # Name
    keyval[[paste("$P",channel_number,"S",sep="")]] <- channel_name;        # Desc
    keyword(out_frame) <- keyval;
    
    write.FCS(out_frame,outfilename);
}

SPADE.downsampleFCS3<-
function (infilename, outfilename, exclude_pctile = 0.01, target_pctile = NULL,
    target_number = NULL, target_percent = 0.1)
{
    in_fcs <- spade:::SPADE.read.FCS(infilename, comp = FALSE, transform = FALSE)
    in_data <- exprs(in_fcs)
    params <- flowCore::parameters(in_fcs)
    pd <- pData(params)
    d_idx <- match("density", pd$name)
    if (is.na(d_idx)) {
        stop("No density parameter in FCS file")
    }
    boundary <- quantile(in_data[, d_idx], c(exclude_pctile,
        target_pctile), names = FALSE)
    out_data <- subset(in_data, in_data[, d_idx] > boundary[1])
    if (!is.null(target_percent)) {
        target_number = round(target_percent * nrow(out_data))
        message("Targeting ", target_number, " events for ",
            infilename)
    }
    density <- out_data[, d_idx]
    if (is.null(target_number)) {
        boundary <- boundary[2]
        out_data <- subset(out_data, boundary/density > runif(nrow(out_data)))
    }
    else if (target_number < nrow(out_data)) {
        density_s <- sort(density)
        cdf <- rev(cumsum(1/rev(density_s)))
        boundary <- target_number/cdf[1]
        if (boundary > density_s[1]) {
            targets <- (target_number - 1:length(density_s))/cdf
            boundary <- targets[which.min(targets - density_s >
                0)]
        }
        out_data <- subset(out_data, boundary/density > runif(length(density)))
    }
    else if (target_number > nrow(out_data)) {
        stop("More events requested than present in file")
    }
    out_frame <- flowFrame(out_data, params, description = description(in_fcs))
    write.FCS(out_frame, outfilename)
}
#####
Modes2 <- function(x) {
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  r=ux[tab == max(tab)]
  return(r)
}

Mode <- function(x) {
  ux <- unique(x)
  r=ux[which.max(tabulate(match(x, ux)))]
  return(r)
}

SPADE.markerMedians2<-function (files, num.clusters, cols = NULL, arcsinh_cofactor = NULL,
    transforms = flowCore::arcsinhTransform(a = 0, b = 0.2),
    cluster_cols = NULL, comp = TRUE)
{
    if (!is.null(arcsinh_cofactor)) {
        warning("arcsinh_cofactor is deprecated, use transform=flowCore::arcsinhTransform(...) instead")
        transforms <- flowCore::arcsinhTransform(a = 0, b = 1/arcsinh_cofactor)
    }
    data <- c()
    files <- as.vector(files)
    for (f in files) {
        in_fcs <- spade:::SPADE.read.FCS(f, comp = comp)
        in_data <- exprs(in_fcs)
        params <- flowCore::parameters(in_fcs)
        pd <- pData(params)
        if (is.null(cols)) {
            cols <- as.vector(pd$name)
        }
        if (!"cluster" %in% cols) {
            cols <- c(cols, "cluster")
        }
        idxs <- match(cols, pd$name)
        if (any(is.na(idxs))) {
            stop("Invalid column specifier")
        }
        data <- rbind(data, in_data[, idxs, drop = FALSE])
    }
    clst <- data[, "cluster"]
    data <- data[, colnames(data) != "cluster", drop = FALSE]
    data_t <- SPADE.transform.matrix(data, transforms)
    colnames(data) <- sapply(colnames(data), function(x) {
        if (x %in% cluster_cols)
            x <- paste(x, "clust", sep = "_")
        x
    })
    colnames(data_t) = colnames(data)
    ids <-1:num.clusters
    if (any(is.na(match(unique(clst), ids)))) {
        print("Cluster ids in FCS files not in graph")
        print(length(ids))
        ids<-1:max(max(unique(clst)),num.clusters)
        print(length(ids))
    }
    count <- matrix(0, nrow = max(ids), ncol = 1, dimnames = list(ids,
        "count"))
        rownames(count)=ids
    medians <- matrix(NA, nrow = max(ids), ncol = ncol(data),
        dimnames = list(ids, colnames(data)))
    raw_medians <- matrix(NA, nrow = max(ids), ncol = ncol(data),
        dimnames = list(ids, colnames(data)))
    cvs <- matrix(NA, nrow = max(ids), ncol = ncol(data),
        dimnames = list(ids, colnames(data)))
    for (i in ids) {
        data_s <- subset(data, clst == i)
        data_s_t <- subset(data_t, clst == i)
        count[i, 1] <- nrow(data_s_t)
        #medians[i, ] <- apply(data_s_t, 2, function(x) median(x,na.rm = TRUE))
        #raw_medians[i, ] <- apply(data_s, 2, function(x) median(x,na.rm = TRUE))
        medians[i, ] <- apply(data_s_t, 2, function(x) medianna0(x))
        raw_medians[i, ] <- apply(data_s, 2, function(x) medianna0(x))
        #cvs[i, ] <- apply(data_s_t, 2, function(d) {
         #   100 * sd(d)/abs(mean(d))
        #})
        cvs[i,]<-0
    }
    percenttotal <- matrix((count/sum(count)) * 100, nrow = max(max(ids),num.clusters),
        ncol = 1, dimnames = list(ids, "percenttotal"))
    list(count = count, medians = medians, raw_medians = raw_medians,
        cvs = cvs, percenttotal = percenttotal)
}

SPADE.write.graph1111<-function (graph, file = "", format = c("gml"))
{
    if (!is.igraph(graph)) {
        stop("Not a graph object")
    }
    if (file == "") {
        file <- stdout()
    }
    else if (is.character(file)) {
        file <- file(file, "w")
        on.exit(close(file))
    }
    else if (!isOpen(file, "w")) {
        open(file, "w")
        on.exit(close(file))
    }
    if (!inherits(file, "connection")) {
        stop("'file' must be a character string or a connection")
    }
    write.gml <- function(graph, file) {
        write.attr <- function(name, attr) {
            name <- gsub("[^A-Za-z0-9_]", "", name)
            if (length(grep("^[0-9]", name))) {
                name <- paste("spade", name, sep = "")
            }
            if (is.na(attr) || is.nan(attr))
                stop("Unexpected NA or NaN attribute")
            else if (is.character(attr) && nchar(attr) > 0)
                paste(name, " \"", attr, "\"", sep = "")
            else if (is.integer(attr))
                paste(name, attr)
            else paste(name, formatC(attr, format = "f"))
        }
        writeLines(c("graph [", paste("directed", ifelse(is.directed(graph),
            1, 0))), con = file)
        v_attr <- list.vertex.attributes(graph)
        v_attr_g <- v_attr[grep("graphics[.]", v_attr)]
        v_attr <- setdiff(v_attr, c(v_attr_g, "id"))
        for (v in V(graph)) {
            writeLines("node [", con = file)
            writeLines(paste("id", v), con = file)
            for (a in v_attr) {
                val <- get.vertex.attribute(graph, a, index = v)
                if (!is.na(val) && !is.nan(val))
                  writeLines(write.attr(a, val), con = file)
            }
            if (length(v_attr_g) > 0) {
                writeLines("graphics [", con = file)
                for (a in v_attr_g) {
                  parts <- unlist(strsplit(a, "[.]"))
                  val <- get.vertex.attribute(graph, a, index = v)
                  if (!is.na(val) && !is.nan(val))
                    writeLines(write.attr(parts[2], val), con = file)
                }
                writeLines("]", con = file)
            }
            writeLines("]", con = file)
        }
        e_attr <- list.edge.attributes(graph)
        if (length(grep("[.]", e_attr)) > 0) {
            stop("Unsupported struct in edge attributes")
        }
        for (e in E(graph)) {
            writeLines("edge [", con = file)
            pts <- get.edges(graph, e)
            writeLines(c(paste("source", pts[1]), paste("target",
                pts[2])), con = file)
            for (a in e_attr) {
                val <- get.edge.attribute(graph, a, index = e)
                if (!is.na(val) && !is.nan(val))
                  writeLines(write.attr(a, val), con = file)
            }
            writeLines("]", con = file)
        }
        writeLines("]", con = file)
    }
    res <- switch(format, gml = write.gml(graph, file), stop(paste("Unsupported output format:",
        format)))
    invisible(res)
}

medianna0<-function(x) {
    if(all(x==0)) {
    c=0
    } else {
        r1=which(x==0)
        if(length(r1)>0) {
        lr1=length(x[-r1])
        lr2=length(x[r1])
         x1=x[-r1]
         if(lr1>lr2) {
            c=median(x,na.rm = TRUE)
         }else{
          c=min(x1,na.rm = TRUE)
         }
       }else{
         c=median(x,na.rm = TRUE)
        }
    }
    #print(paste("median=",c))
    return(c)
}

SPADE.read.FCS <- function(file, comp=TRUE, verbose=FALSE, ...) {
    if (verbose)
        fcs <- read.FCS(file, ...)
    else
        fcs <- suppressWarnings(read.FCS(file, ...))
    
    params <- parameters(fcs)
    pd     <- pData(params)

    # Replace any null descs with names (for FSC-A, FSC-W, SSC-A)
    bad_col <- grep("^[a-zA-Z0-9]+",pd$desc,invert=TRUE)
    if (length(bad_col) > 0) {
        keyval <- keyword(fcs)
        for (i in bad_col) {
            pd$desc[i] <- pd$name[i]
            keyval[[paste("$P",i,"S",sep="")]] <- pd$name[i]
        }
        pData(params) <- pd;
        fcs <- flowFrame(exprs(fcs),params,description=description(fcs));
        keyword(fcs) <- keyval
    }

    # Compensate data if SPILL or SPILLOVER present, stripping compensation matrix
    # out of the flowFrame, i.e we should only have to do this once
    apply.comp <- function(in_fcs, keyword) {
        comp_fcs <- compensate(in_fcs, description(in_fcs)[[keyword]])
        flowFrame(exprs(comp_fcs), parameters(comp_fcs), description(comp_fcs)[grep("SPILL",names(description(comp_fcs)),invert=TRUE)])
    }
    
    if (comp && !is.null(description(fcs)$SPILL)) {
        fcs <- apply.comp(fcs, "SPILL")
    } else if (comp && !is.null(description(fcs)$SPILLOVER)) {
        fcs <- apply.comp(fcs, "SPILLOVER")
    }
    
    fcs
}

SPADE.transform.matrix <- function(mat, tform=NULL) {
    if (is.null(tform)) {
        mat  # No-op
    } else {
        if (class(tform) == "transform") {
            apply(mat, 2, tform)
        } else {
            for (name in intersect(colnames(mat), names(tform))) {
                mat[,name] <- tform[[name]](mat[,name])
            }
            mat
        }
    }
}

SPADE.transform.FCS <- function(ff, tform=NULL) {
    if (is.null(tform)) {
        ff  # No-op
    } else {

        if (class(tform) == "transform") {
            new_exprs <- apply(exprs(ff), 2, tform)
            new_range <- apply(range(ff), 2, tform)
        } else {
            new_exprs <- exprs(ff)
            new_range <- range(ff)
            for (name in intersect(colnames(ff), names(tform))) {
                new_exprs[,name] <- tform[[name]](new_exprs[,name])
                new_range[,name] <- tform[[name]](new_range[,name])
            }
            new_range <- as.matrix(new_range)
        }
                
        new_par <- parameters(ff)
        new_par$minRange <- new_range[1,]
        new_par$maxRange <- new_range[2,]
        new_par$range    <- (new_range[2,] - new_range[1,]) + 1

        flowFrame(new_exprs, new_par, description(ff))
    }
}

SPADE.build.flowFrame <- function(x) {
    if (!is.matrix(x)) {
        stop("Input must be matrix")
  }

  # Build metadata for FCS file
  pd <- c()  # 'params' phenoData
  dl <- list()  # 'description' list
    
  dl[["$DATATYPE"]] <- "F"
  for (c in 1:ncol(x)) {
        c_name <- colnames(x)[c]
            
        c_min <- min(x[,c])
        c_max <- max(x[,c])
        c_rng <- c_max - c_min + 1

        pl <- matrix(c(c_name, c_name, c_rng, c_min, c_max),nrow=1)
        colnames(pl) <- c("name", "desc", "range", "minRange", "maxRange")
        rownames(pl) <- paste("$P",c,sep="")
        pd <- rbind(pd, pl)
        
        dl[[paste("$P",c,"B",sep="")]] <- "32";        # Number of bits
        dl[[paste("$P",c,"R",sep="")]] <- toString(c_rng); # Range
        dl[[paste("$P",c,"E",sep="")]] <- "0,0";        # Exponent
        dl[[paste("$P",c,"N",sep="")]] <- c_name;        # Name
        dl[[paste("$P",c,"S",sep="")]] <- c_name;        # Desc
    }
   
  flowFrame(x, as(data.frame(pd), "AnnotatedDataFrame"), description=dl)
}

SPADE.removeExistingDensityAndClusterColumns <- function(file) {
    # Do not comp or transform ... make this step invisible.
    input_file <- suppressWarnings(read.FCS(file))

    input_file_names <- names(input_file)

    if ("<cluster> cluster" %in% input_file_names ||
        "<density> density" %in% input_file_names) {
        # Drop those columns
        cleaned <- input_file[,!(input_file_names %in% c("<cluster> cluster", "<density> density"))]

        # Rename the original file. Increment the suffix if it already exists.
        suffix <- as.integer(regmatches(file, gregexpr("(\\d+)$", file, perl=TRUE)))

        if (is.na(suffix)) {
            suffix <- ".orig1"
            new_file_name <- paste(file, suffix, sep = "")
            # NB: file.rename is a namespaced function, not a method of the argument.
            file.rename(file, new_file_name)
        } else {
            suffix <- paste(".orig", suffix + 1, sep = "")
            new_file_name <- sub("(.orig\\d+)$", suffix, file);
            file.rename(file, new_file_name)
        }

        # Save with the original file name
        write.FCS(cleaned, file)
    }
}
