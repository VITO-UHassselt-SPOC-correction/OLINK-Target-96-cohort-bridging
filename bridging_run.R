# Just put your raw NPX data already in 'data' folder inside the repository so that the bridging_fun.R script can read the files.
source("bridging_fun.R")
#filename1 <- "./data/yourdata_NPX.xlsx"
#filename2 <- "./data/yourdata_NPX.xlsx"
#Additional batches can be included completely analogue to filename 1 and filename 2

#load the reference pool that includes all proteins that OLINK measures with target 96.
filename4 <- "./data/output/pooleffect_scriptinput.xlsx"

#additional filenames should be included below
run_fun(filename1, filename2)




