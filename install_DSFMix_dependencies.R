get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}
########
r_version_int <- 0
if (as.numeric(R.Version()$minor) < 1.0){
  r_version_int <- as.numeric(R.Version()$minor )+ as.numeric(R.Version()$major)
}else{
  minor <- as.numeric(R.Version()$minor)
  if (minor >=1.0)
    minor=minor*0.1
  r_version_int <- minor + as.numeric(R.Version()$major)
}

isBiocManager <- FALSE

# R version >= 3.5 we should use BiocManager::install instead of biocLite
# function
if (r_version_int < 3.5){
  if (!require(biocLite, quietly = TRUE)){
    source("http://bioconductor.org/biocLite.R")
  }
}else{
  isBiocManager <- TRUE
  if (!require(BiocManager)){
    install.packages("BiocManager", repos = "https://cran.rstudio.com/", dependencies = TRUE)
  }
}

# We need also the package remotes to be installed
if(!require(remotes)){
   install.packages("remotes",repos = "https://cran.rstudio.com/", dependencies = TRUE)
}

### List of needed packages

# We do not include plot3D because it does not load, and halts the entire session
needed_packages <- c("spade","Seurat", "xgboost", "car", "scales", "flowCore", "RColorBrewer", 
                     "cluster", "fastcluster", "sp", "igraph", "dendextend", "ggplot2", "plot3D",
                     "diptest", "mutoss", "lawstat","robustbase", "uwot", 
                     "Rtsne", "phateR", "party", "genie", "dendsort", "networkD3", 
                     "magrittr", "htmlwidgets", "plyr", "tidyr", "dynamicTreeCut", 
                     "gTests", "parallel","modeltools", "sandwich", "strucchange", "genieclust","Matrix")

# Before we continue further, if the OS is mac, we need to install odl misc3d version
# Since the current one DOES NOT load on macos and makes plot3D not to load either
os <- get_os()
if(os=="osx"){
   # delete curently installe misc3d package
  if(require(misc3d))
      remove.packages("misc3d")
  
  # install old version of misc3d provided in the repository
  install.packages("misc3d_0.8-4.tar.gz",repos = NULL, type = "source")
}

# List of Bioconductor packages, including the installed ones
if (isBiocManager){
  avail_list <- BiocManager::available(include_installed = TRUE)
}else{
  avail_list <- all_group()
}
# Distinguish those from packages from Bioconductor from those from CRAN
needed_bioc_pkgs <- needed_packages[needed_packages %in% avail_list]
needed_CRAN_pkgs <- needed_packages[! needed_packages %in% avail_list]

# xgboost is not Bioconductor although is listed as such
needed_bioc_pkgs <- needed_bioc_pkgs[-grep("xgboost",needed_bioc_pkgs)]
needed_CRAN_pkgs <- c("xgboost",needed_CRAN_pkgs)

# Install those needed BIOConductor packages
new.packages <- needed_bioc_pkgs[!(needed_bioc_pkgs %in%  installed.packages()[,"Package"])]
#new.packages <- needed_bioc_pkgs
if(length(new.packages)){
  for (i in 1:length(new.packages)){
    pck <- new.packages[i]
    
    M <- NULL
    dotR <- file.path(Sys.getenv("HOME"), ".R")
    if (!file.exists(dotR))
      dir.create(dotR)
    M <- file.path(dotR, "Makevars")
    if (!file.exists(M))
      file.create(M)
    
    cat("\nLIB_DIR=$(LD_LIBRARY_PATH)",
        "INCLUDE_DIR=$(CPATH)",
        file = M, sep = "\n", append = F)
    
    
    print(paste0("Installing package ",pck,""))
    
    if(pck=="Seurat"){
      print(paste0("Installing package ",pck,".."))
      remotes::install_github("satijalab/seurat", dependencies=TRUE, upgrade="never",args="--configure.vars='INCLUDE_DIR=$(CPATH) LIB_DIR=$(LD_LIBRARY_PATH)'")
      
    }else if(pck=="flowCore" && isBiocManager){
      # flowCore dependencies have to be installed directly from github
      # 1. First install latest version of RProtoBufLib 2.3.1
      remotes::install_github("RGLab/RProtoBufLib", args="--configure.vars='INCLUDE_DIR=$(CPATH) LIB_DIR=$(LD_LIBRARY_PATH)'")
      
      # 2. Install cytolib without using checking the dependencies, otherwise it always tries to install RProtoBufLib 2.2.0 which the 
      # default version at Bioconductor repository
      remotes::install_github("RGLab/cytolib", upgrade="never", dependencies=FALSE, args="--configure.vars='INCLUDE_DIR=$(CPATH) LIB_DIR=$(LD_LIBRARY_PATH)'")
      
      # 3. Install flowCore package itself:
      BiocManager::install(pck, args="--configure.vars='INCLUDE_DIR=$(CPATH) LIB_DIR=$(LD_LIBRARY_PATH)'")
      
    }else if(pck=="party"){
      # party is from CRAN, and for some reason the dependencies cannot be handled by install.packages
      # so we installed them first
      pack <- available.packages(repos="https://cran.rstudio.com/")
      dependencies <- gsub("\\s+|\\n","",gsub("\\s+\\(.+","",unlist(strsplit(split=",",x=pack[pck,"Depends"])),perl=TRUE),perl=TRUE)
      dependencies <- dependencies[ ! dependencies %in% "R"]
      
      # installing dependeencies
      install.packages(dependencies,repos="https://cran.rstudio.com/", dependencies = TRUE,force = TRUE, args="--configure.vars='INCLUDE_DIR=$(CPATH) LIB_DIR=$(LD_LIBRARY_PATH)'")
      
      # And then we install the party package itself
      BiocManager::install(pck, args="--configure.vars='INCLUDE_DIR=$(CPATH) LIB_DIR=$(LD_LIBRARY_PATH)'")
    }else{
      
      if(isBiocManager){
        BiocManager::install(pck, dependencies=TRUE, args="--configure.vars='INCLUDE_DIR=$(CPATH) LIB_DIR=$(LD_LIBRARY_PATH)'")
      }else{
        biocLite(pck, dependencies = TRUE,args="--configure.vars='INCLUDE_DIR=$(CPATH) LIB_DIR=$(LD_LIBRARY_PATH)'")
      }
    }
    
  }
}
# Install those needed CRAN packages
new.packages <- needed_CRAN_pkgs[!(needed_CRAN_pkgs %in%  installed.packages()[,"Package"])]
new.packages <- new.packages[! new.packages %in% "parallel"] # parallel is now a base package, no installation is needed
if(length(new.packages)){
  for (i in 1:length(new.packages)){
    pck <- new.packages[i]
    print(paste0("Installing package ",pck,""))
    
    M <- NULL
    dotR <- file.path(Sys.getenv("HOME"), ".R")
    
    if (!file.exists(dotR))
      dir.create(dotR)
    M <- file.path(dotR, "Makevars")
    if (!file.exists(M))
      file.create(M)
    if(pck=="spade"){
      if(grepl("apple",R.Version()$platform,perl = TRUE, ignore.case = TRUE)){
        # if it is MacOSX it requires CLANG
        cat("\nC=/usr/bin/clang",
            "CXX=/usr/bin/clang++",
            "LIB_DIR=$(LD_LIBRARY_PATH)",
            "INCLUDE_DIR=$(CPATH)",
            file = M, sep = "\n", append = F)
      }else{
        cat("\nLIB_DIR=$(LD_LIBRARY_PATH)",
            "INCLUDE_DIR=$(CPATH)",
            file = M, sep = "\n", append = F)
        
      }
      # install spade from source file in gitfront
      #install.packages("https://gitfront.io/r/anchangslab/431118cea12808ca2a76813342377dc01eb26ebe/DSFMix/blob/spade_1.0.0.tar", repos = NULL, type = "source")
      install.packages("spade_1.0.0.tar", repos = NULL, type = "source")
      
    }else if(pck=="xgboost"){
      cat("\nCXX14FLAGS=-fPIC",
          "CXX14 = $(HOME)/Work/Software/bin/g++ -std=c++1y",
          "CXX11FLAGS=-fPIC",
          "LIB_DIR=$(LD_LIBRARY_PATH)",
          "INCLUDE_DIR=$(CPATH)",
          file = M, sep = "\n", append = F)
      install.packages(pck,repos="https://cran.rstudio.com/", dependencies = TRUE,args="--configure.vars='INCLUDE_DIR=$(CPATH) LIB_DIR=$(LD_LIBRARY_PATH) CXX14 = $(HOME)/Work/Software/bin/g++ -std=c++1y'")
      
    }else{
      cat("\nLIB_DIR=$(LD_LIBRARY_PATH)",
          "INCLUDE_DIR=$(CPATH)",
          file = M, sep = "\n", append = F)
      
      install.packages(pck,repos="https://cran.rstudio.com/", dependencies = TRUE,args="--configure.vars='INCLUDE_DIR=$(CPATH) LIB_DIR=$(LD_LIBRARY_PATH) CXX14 = $(HOME)/Work/Software/bin/g++ -std=c++1y'")
      
    }
    
    
    if(!is.null(M) && file.exists(M) && pck=="xgboost"){
      file.remove(M)
    }
  }
}

# Loading those libraries (remember spade has been installed externally with )
load_libraries <- needed_packages

for (pck in load_libraries){
  print(paste0("Loading package ",pck,"..."))
  if (pck == "party"){
    # Load dependencies first
    library("modeltools",character.only = TRUE)
    library("sandwich",character.only = TRUE)
    library("strucchange",character.only = TRUE)
    
  }else if(pck == "genie"){
    library("genieclust",character.only = TRUE)
  }else if(pck == "uwot"){
    library("Matrix",character.only = TRUE)
  }
  
  library(pck,character.only = TRUE)
  
}
