# or this if the data already in 'data' folder inside the repository
source("bridging_fun1.R")
filename1 <- "./data/yourdata_NPX.xlsx"
filename2 <- "./data/yourdata_NPX.xlsx"
#Additional batches can be included completely analogue to filename 1 and filename 2

#load the reference pool that includes all proteins that OLINK measures with target 96.
filename3 <- "./data/reference pool/pooleffect_prep_correct.xlsx"
filename4 <- "./data/output/pooleffect_scriptinput.xlsx"

#additional filenames should be included below
run_fun(filename1, filename2, filename3)




