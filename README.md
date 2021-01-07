# DSFMix

RUNNING DSFMix in R
1)Open R
2)Install old version spade (spade_1.0.0.tar) included in the folder
###COMMENT on installing spade (spade_1.0.0.tar)##
#installing spade algorithm (spade_1.0.0.tar) might produce C++ clang related errors. I recommend open terminal and run the following before installation:
#cat <<- EOF > ~/.R/Makevars
#C=/usr/bin/clang
#CXX=/usr/bin/clang++
#EOF
#R CMD INSTALL spade_1.0.0.tar
2)Install all other dsfmix dependent R libraries highlighted in the file "feature dsfmix prerequisite.R"

3)For each dsfmix application create an input folder with the following input data:
a)dat1.rdata #scRNA-seq normalized data matrix (dat1_****.rdata) with rows genes and columns cells
b)seuratclust # vector of cluster labels associated with each cell (seuratclust_****.rdata)
c)timeclust # vector of time point labels associated with each cell (timeclust_****.rdata)
d)feature dsfmix prerequisite.R
e)feature dsfmix running code.R
f)feature dsfmix spadefunctions.R
g)fcstospadedatamatrix.fcs

4)Set the input folder as working directory to run dsfmix
5)In line 16 of the file "feature dsfmix running code.R" set file name for spade input FCS data as spadefilename="dsfmix***.fcs" E.g. for hormonedata set spadefilename="dsfmixhormone.fcs"
6)In line 428 of the file "feature dsfmix running code.R" set number of spade initial clusters to k=300 for spermdata, k=200 for ipscdata and emtdata and k=100 for hormonedata. Also due to small sample size for dsfmix analysis on hormonedata set the parameter do_real_filtering = FALSE in the dsfmixdriverc function.
7)In R run dsfmix:
   source("feature dsfmix running code.R")
8)After a couple of minutes to hours depending on the size of data, in the input folder you can view several output figures and r objects corresponding to several output data associated with different stages in  dsfmix algorithm.
9)Within the input folder there is an output folder named "output", which contains a) the dsfmix input data stored as FCS file and other spade related files
    b) the results of  spade found within the output folder called "orig_spade" called "pdf" showing the profile of each gene on a spade tree
    c) the results of dsfmix forest found within the orig_spade folder called "spadeforest" showing the profile of each gene on a dynamic spanning forest.
