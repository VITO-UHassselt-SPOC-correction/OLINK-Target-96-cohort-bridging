#A VITO/DSI workflow to apply a Synthetic Plasma pOol Cohort (SPOC) correction for OLINK proteomics data
#required packages
install.packages("readxl")
library(dplyr)
library(stringr)
library(lsa)
library(ggplot2)
library(data.table)
library(dplyr)
library(tidyr)
library(tidyverse)
library(reshape2)
library(pastecs)
library(readxl)

#### FIRST STEP: Function to read NPX data into long format
#'
#' Imports an NPX file exported from NPX Manager. 
#' No alterations to the output NPX Manager format are allowed.
read_NPX <- function(filename){
  
  # If the file is csv or txt, read_NPX assumes that the file is explore data in long format
  if (tools::file_ext(filename) %in% c("csv","txt")) {
    #read file using ; as delimiter
    out <- read.table(filename, header = T, sep=";", stringsAsFactors = F,
                      na.strings = c("NA",""))
    if (is.data.frame(out) & ncol(out) == 1) {
      #if only one column in the data, wrong delimiter. use , as delimiter
      out <- read.table(filename, header = T, sep=",", stringsAsFactors = F,
                        na.strings = c("NA",""))
    }
    #check that all colnames are present
    match_old_header <- all(c("SampleID", "Index", "OlinkID", "UniProt", "Assay",
                              "MissingFreq", "Panel", "Panel_Version", "PlateID",
                              "QC_Warning", "LOD", "NPX") %in% colnames(out))
    match_new_header <- all(c("SampleID", "Index", "OlinkID", "UniProt", "Assay",
                              "MissingFreq", "Panel", "Panel_Lot_Nr", "PlateID",
                              "QC_Warning", "LOD", "NPX") %in% colnames(out))
    if (match_old_header || match_new_header) {
      out$NPX <- as.numeric(out$NPX)
      out$LOD <- as.numeric(out$LOD)
      out$MissingFreq <- as.numeric(out$MissingFreq)
      out$SampleID <- as.character(out$SampleID)
      return(as_tibble(out))
    } else {
      missing.cols <- setdiff(c("SampleID", "Index", "OlinkID", "UniProt", "Assay", "MissingFreq", "Panel", "Panel_Lot_Nr", "PlateID", "QC_Warning", "LOD", "NPX"),
                              colnames(out))
      #If columns are missing, stop and print out which are missing
      stop(paste0("Cannot find columns ", paste(missing.cols,collapse=",")))
    }
  }
  NORM_FLAG <-  F
  
  # Check if the data is npx or concentration as well as if it is tg48 or tg96
  
  data_type <- readxl::read_excel(filename, range='A2',
                                  col_names = F, .name_repair="minimal")
  if (grepl('NPX', data_type, fixed=TRUE)) {
    is_npx_data <- TRUE
    n_max_meta_data <- 4
    
    # Check whether it is target 48 or 96
    panel_name <- readxl::read_excel(filename, range='B3',
                                     col_names = F, .name_repair="minimal")
    if (grepl('Target 48', panel_name, fixed=TRUE)) {
      target_type <- '48'
      BASE_INDEX <- 45
    } else {
      target_type <- '96'
      BASE_INDEX <- 92
    }
  } else if (grepl('Quantified', data_type, fixed=TRUE)) {
    # Quant data given, which also means it is target 48
    is_npx_data <- FALSE
    n_max_meta_data <- 5
    target_type <- '48'
    BASE_INDEX <- 45
  } else {
    stop("Cannot find whether the given data is NPX or concentration")
  }
  
  # Load initial meta data (the first rows of the wide file)
  meta_dat <-  readxl::read_excel(filename, skip = 2, n_max = n_max_meta_data,
                                  col_names = F, .name_repair="minimal")
  meta_dat[4,1] <- 'SampleID'
  NR_DEVIATIONS <- sum(stringr::str_detect(meta_dat[2,],
                                           'QC Deviation from median'))
  control_index <- (stringr::str_detect(meta_dat[2,], 'Det Ctrl') |
                      stringr::str_detect(meta_dat[2,], 'Inc Ctrl 2') |
                      stringr::str_detect(meta_dat[2,], 'Inc Ctrl 1') |
                      stringr::str_detect(meta_dat[2,], 'Ext Ctrl'))
  meta_dat[4, control_index] <- meta_dat[2, control_index]
  meta_dat[3, control_index] <- '-'
  NR_CONTROLS <- sum(control_index)
  nr_panel<-(ncol(meta_dat)-1-NR_DEVIATIONS-NR_CONTROLS)/(BASE_INDEX+2)
  
  nr_col <- ncol(meta_dat)
  names(meta_dat) <- as.character(1:nr_col)
  
  meta_dat <- meta_dat %>%
    rename(Name = `1`)
  
  # Load NPX or QUANT data including the last rows of meta data
  dat <- readxl::read_excel(filename, skip = n_max_meta_data+2, col_names = F,
                            .name_repair="minimal", col_types = c('text'))
  
  nr_col <- ncol(dat)
  names(dat) <- as.character(1:nr_col)
  
  dat<-dat %>%
    rename(Name = `1`)
  
  # Calc nbr of plates
  plates <- dat[,nr_col-nr_panel] %>% distinct() %>% na.omit() %>% pull()
  nr_plates <- length(plates)
  
  # Extract the meta data from the last rows of data
  missfreq<-dat %>% filter(stringr::str_detect(Name, "Missing Data freq."))
  norm_method <- dat %>% filter(stringr::str_detect(Name, "Normalization"))
  if (!is_npx_data) {
    assay_warning <- dat %>% filter(stringr::str_detect(Name, "Assay warning"))
    Plate_LQL <- dat %>% filter(stringr::str_detect(Name, 
                                                    "Lowest quantifiable level"))
    LOD <- dat %>% filter(stringr::str_detect(Name, "Plate LOD"))
    LLOQ <- dat %>% filter(stringr::str_detect(Name, "LLOQ"))
    ULOQ <- dat %>% filter(stringr::str_detect(Name, "ULOQ"))
  } else {
    LOD <- dat %>% filter(stringr::str_detect(Name, "LOD"))
  }
  
  # Add the new meta data to Â´meta_datÂ´
  meta_dat <- rbind(meta_dat,missfreq)
  if (!is_npx_data) {
    meta_dat <- rbind(meta_dat,LLOQ,ULOQ,assay_warning,Plate_LQL)
  }
  meta_dat <- rbind(meta_dat,LOD,norm_method)
  
  # Remove the meta data from dat
  if (is_npx_data) {
    nbr_meta_data_rows_bottom <- 3
  } else {
    nbr_meta_data_rows_bottom <- 4+3*nr_plates
  }
  if (nrow(norm_method) == 0) {
    nbr_meta_data_rows_bottom <- nbr_meta_data_rows_bottom - 1
  } else {
    NORM_FLAG <- T
  }
  dat <- dat[c(-1*(nrow(dat) - nbr_meta_data_rows_bottom):nrow(dat)),]
  
  # Create index vector
  SampleID <- dat$Name
  Index_nr <- c(1:length(SampleID))
  
  # Initiate lists for later use
  panel_data <- list() ##NPX values to be stored
  QC_list <- list()    ##QC data
  meta_data_list <- list() ## meta data
  panel_list <- list()  ## combination of panel data and QC
  assay_name_list <- list()
  panel_list_long <- list()
  deviations_list <- list()
  
  if (NR_CONTROLS > 0) {
    BASE_INDEX <- BASE_INDEX + NR_CONTROLS/nr_panel
  }
  
  # Construct a list of tibbles that match the long format
  for (i in 1:nr_panel) {
    
    panel_data[[i]]<-dat[,(2+((i-1)*BASE_INDEX)):((BASE_INDEX+1)+((i-1)*BASE_INDEX))]
    
    if (NR_DEVIATIONS == 0) {
      QC_list[[i]]<-dat[,c((2+((nr_panel)*BASE_INDEX)+(i-1)),
                           (2+((nr_panel)*BASE_INDEX)+(i-1))+nr_panel)]
      
      meta_data_list[[i]]<-meta_dat[,c((2+((i-1)*BASE_INDEX)):((BASE_INDEX+1)+((i-1)*BASE_INDEX)),
                                       (2+((nr_panel)*BASE_INDEX)+(i-1)),
                                       (2+((nr_panel)*BASE_INDEX)+(i-1))+nr_panel)]
      
      
    } else {
      
      QC_list[[i]]<-dat[,c((2+((nr_panel)*BASE_INDEX)+(i-1)),
                           (2+((nr_panel)*BASE_INDEX)+(i-1))+nr_panel, 
                           (2+((nr_panel)*BASE_INDEX)+(i-1))+2*nr_panel+(i-1), 
                           (2+((nr_panel)*BASE_INDEX)+(i-1))+2*nr_panel+(i-1)+1)]
      
      meta_data_list[[i]]<-meta_dat[,c((2+((i-1)*BASE_INDEX)):((BASE_INDEX+1)+((i-1)*BASE_INDEX)),
                                       (2+((nr_panel)*BASE_INDEX)+(i-1)),
                                       (2+((nr_panel)*BASE_INDEX)+(i-1))+nr_panel,
                                       (2+((nr_panel)*BASE_INDEX)+(i-1))+2*nr_panel, 
                                       (2+((nr_panel)*BASE_INDEX)+(i-1))+3*nr_panel)]
      
      meta_data_list[[i]][4,(BASE_INDEX+3)] <- "QC Deviation Inc Ctrl"
      meta_data_list[[i]][4,(BASE_INDEX+4)] <- "QC Deviation Det Ctrl"
      
      
    }
    
    meta_data_list[[i]][4,(BASE_INDEX+1)] <- meta_data_list[[i]][2,(BASE_INDEX+1)]
    meta_data_list[[i]][4,(BASE_INDEX+2)] <- meta_data_list[[i]][2,(BASE_INDEX+2)]
    
    
    panel_list[[i]]<-cbind(panel_data[[i]],QC_list[[i]])
    
    colnames(panel_list[[i]]) <- unlist(meta_data_list[[i]][4,])
    
    panel_list[[i]][,c(-(BASE_INDEX+1),-(BASE_INDEX+2))] <- lapply(panel_list[[i]][,c(-(BASE_INDEX+1),-(BASE_INDEX+2))],
                                                                   function(x) as.numeric(stringr::str_replace_all(x, c('#' = '', ',' = '.', 'No Data' = NA, '> ULOQ' = NA, '< LLOQ' = NA))))
    
    # Remove the last two columns since they contain redundant meta data and
    # will only cause warnings
    meta_data_list[[i]] <- meta_data_list[[i]][,c(-(BASE_INDEX+1),-(BASE_INDEX+2))]
    
    assay_name_list[[i]]<-tibble(ID=c(t(meta_data_list[[i]][4,])),
                                 Name=c(t(meta_data_list[[i]][2,])),
                                 UniProt = c(t(meta_data_list[[i]][3,])),
                                 Panel=c(t(meta_data_list[[i]][1,])))
    
    if (is_npx_data) {
      assay_name_list[[i]] <- bind_cols(assay_name_list[[i]],
                                        MissingFreq=c(t(meta_data_list[[i]][5,])),
                                        LOD = as.numeric(c(t(meta_data_list[[i]][6,]))))
      
      if (NORM_FLAG == TRUE) {
        assay_name_list[[i]] <- bind_cols(assay_name_list[[i]], 
                                          Normalization = c(t(meta_data_list[[i]][7,])))
      }
      
      panel_list_long[[i]]<- panel_list[[i]] %>%
        mutate(SampleID = SampleID) %>%
        mutate(Index = Index_nr) %>%
        gather(Assay, NPX, -SampleID,-`QC Warning`,-`Plate ID`,-Index,-matches("QC Deviation Inc Ctrl"), -matches("QC Deviation Det Ctrl")) %>%
        left_join(assay_name_list[[i]], by = c('Assay' = 'ID')) %>%
        select(SampleID, Index, Assay, UniProt, Name, MissingFreq, Panel,`Plate ID`,`QC Warning`, LOD, NPX, matches("Assay_Warning"), matches("Normalization"), matches("QC Deviation Inc Ctrl"), matches("QC Deviation Det Ctrl")) %>%
        rename(PlateID = `Plate ID`) %>%
        rename(QC_Warning = `QC Warning`) %>%
        rename(OlinkID = Assay, Assay = Name) 
    } else {
      for (j in 1:nr_plates) {
        assay_name_by_plate <- bind_cols(assay_name_list[[i]],
                                         Unit=c(t(meta_data_list[[i]][5,])),
                                         MissingFreq=c(t(meta_data_list[[i]][6,])),
                                         LLOQ = as.numeric(c(t(meta_data_list[[i]][7,]))),
                                         ULOQ = as.numeric(c(t(meta_data_list[[i]][8,]))),
                                         Assay_Warning = c(t(meta_data_list[[i]][(9+(j-1)),])),
                                         Plate_LQL = as.numeric(c(t(meta_data_list[[i]][(9+nr_plates+(j-1)),]))),
                                         LOD = as.numeric(c(t(meta_data_list[[i]][(9+2*nr_plates+(j-1)),]))))
        if(NORM_FLAG == T){
          assay_name_by_plate <- bind_cols(assay_name_by_plate, 
                                           Normalization = c(t(meta_data_list[[i]][9+3*nr_plates,])))
        }
        panel_list_long[[(i-1)*j+j]]<- panel_list[[i]] %>%
          mutate(SampleID = SampleID) %>%
          mutate(Index = Index_nr) %>%
          gather(Assay, NPX, -SampleID,-`QC Warning`,-`Plate ID`,-Index,-matches("QC Deviation Inc Ctrl"), -matches("QC Deviation Det Ctrl")) %>%
          filter(`Plate ID` == plates[[j]]) %>% 
          left_join(assay_name_by_plate, by = c('Assay' = 'ID')) %>%
          select(SampleID, Index, Assay, UniProt, Name, MissingFreq, Panel,`Plate ID`,`QC Warning`, LOD, NPX, Unit, ULOQ, LLOQ, Plate_LQL, Assay_Warning, matches("Normalization"), matches("QC Deviation Inc Ctrl"), matches("QC Deviation Det Ctrl")) %>%
          rename(PlateID = `Plate ID`) %>%
          rename(QC_Warning = `QC Warning`) %>%
          rename(OlinkID = Assay, Assay = Name)
      }
    }
  }
  
  if (!is_npx_data) {
    for (i in 1:(nr_panel*nr_plates)) {
      panel_list_long[[i]] <- panel_list_long[[i]] %>%
        rename(Plate_LOD = LOD) %>%
        rename(Quantified_value = NPX)
    }
  }
  
  bind_rows(panel_list_long) %>%
    filter(!is.na(SampleID)) %>%
    as_tibble() %>% 
    mutate(Panel_Version = gsub(".*\\(","",Panel)) %>% 
    mutate(Panel_Version = gsub("\\)","",Panel_Version)) %>% 
    mutate(Panel =  gsub("\\(.*\\)","",Panel)) %>% 
    select(SampleID, Index, OlinkID, UniProt, Assay, MissingFreq, Panel,Panel_Version,PlateID, QC_Warning,matches("Plate_LQL"),matches("LOD"),matches("Plat_LOD"),matches("LLOQ"),matches("ULOQ"),matches("NPX"),matches("Quantified_value"),matches("Unit"),matches("Assay_Warning"),matches("Normalization"), matches("*Inc Ctrl*"), matches("*Det Ctrl*"))
  
}

# Create olink_normalization function that will perform this batch correction 
olink_normalization <- function(df1,
                                df2 = NULL,
                                overlapping_samples_df1,
                                overlapping_samples_df2 = NULL,
                                df1_project_nr = 'P1',
                                df2_project_nr = 'P2',
                                reference_project = 'P1',
                                reference_medians = NULL) {
  
  #Filtering on valid OlinkID
  df1 <- df1 %>%
    filter(stringr::str_detect(OlinkID,
                               "OID[0-9]{5}"))
  
  if(!is.null(df2)){
    
    df2 <- df2 %>%
      filter(stringr::str_detect(OlinkID,
                                 "OID[0-9]{5}"))
  }
  
  #median of difference flag
  MOD_FLAG <- T
  
  if(!missing(overlapping_samples_df1)){
    
    if(!all(overlapping_samples_df1 %in% df1$SampleID)){
      
      missing.sample.ids <- setdiff(overlapping_samples_df1,df1$SampleID)
      stop(paste0("SampleID's not found in df1: ",paste(missing.sample.ids,collapse=", ")))
      
    }
    
    if(is.null(overlapping_samples_df2)){
      
      #Testing for reference median normalization, working only on df1
      if(!is.null(reference_medians)){
        
        message('Reference median normalization will be performed.')
        
        shared_oids <- intersect(df1$OlinkID, reference_medians$OlinkID)
        
        missing_df1_oid <- setdiff(df1$OlinkID, shared_oids)
        
        if(length(missing_df1_oid) > 0){
          warning(paste0("There are no reference medians for these assays in df1: ",
                         paste(missing_df1_oid,collapse=", "),
                         '. They will not be normalized.'))
          
        }
        
        
        missing_ref_oid <- setdiff(reference_medians$OlinkID, shared_oids)
        
        if(length(missing_ref_oid) > 0){
          warning(paste0("These OlinkID:s from the reference_medians are not contained in df1: ",
                         paste(missing_df1_oid,collapse=", "),
                         '. They will not be used for normalization.'))
        }
        
        
        adj_factor_df <- df1 %>%
          filter(SampleID %in% overlapping_samples_df1) %>%
          left_join(reference_medians, by= c('OlinkID')) %>%
          group_by(OlinkID) %>%
          mutate(Assay_Median=median(NPX, na.rm = T)) %>%
          ungroup() %>%
          mutate(Adj_factor = if_else(is.na(Reference_NPX - Assay_Median),
                                      0,
                                      Reference_NPX - Assay_Median)) %>%
          select(OlinkID, Adj_factor) %>%
          distinct()
        
        
        df_adjusted_data <- df1 %>%
          mutate(Panel=str_replace(Panel,'\\(.+', '')) %>%
          left_join(adj_factor_df,by='OlinkID') %>%
          mutate(NPX = NPX + Adj_factor) %>%
          mutate(LOD = LOD + Adj_factor)
        
        return(df_adjusted_data)
        
      }
      
      #MOD bridge normalization
      if(!all(overlapping_samples_df1 %in% df2$SampleID)){
        
        missing.sample.ids <- setdiff(overlapping_samples_df1,df2$SampleID)
        stop(paste0("SampleID's not found in df2: ",paste(missing.sample.ids,collapse=", ")))
      }
      
    }else{
      
      MOD_FLAG <- F
      
      if(!all(overlapping_samples_df2 %in% df2$SampleID)){
        
        missing.sample.ids <- setdiff(overlapping_samples_df2, df2$SampleID)
        stop(paste0("SampleID's not found in df2: ",paste(missing.sample.ids,collapse=", ")))
      }
    }
    
  }else{
    
    stop('An array of SampleID:s (overlapping_samples_df1) needs to be provided for normalization. ')
    
  }
  
  
  if(!(reference_project == df1_project_nr) & !(reference_project == df2_project_nr)){
    
    stop('The reference project needs to be one of the included project numbers in variables df1_project_nr or df2_project_nr.')
    
  }
  
  
  if(MOD_FLAG){
    
    message('Bridging normalization with overlapping samples will be performed.')
    
  }else{
    
    message('Bridging normalization with subcohorts or intensity normalization will be performed.')
    
  }
  
  df1<-df1 %>%
    mutate(Project = df1_project_nr)
  
  df2<-df2 %>%
    mutate(Project = df2_project_nr)
  
  if(MOD_FLAG){
    
    #Median of differences
    
    #Calculate adjustment factors
    adj_factor_df<-df1 %>%
      rbind(df2) %>%
      filter(SampleID %in% overlapping_samples_df1) %>%
      select(SampleID,OlinkID,UniProt,NPX,Project) %>%
      spread(Project,NPX)
    
    if (reference_project == df1_project_nr) {
      adj_factor_df <- adj_factor_df %>%
        mutate(Diff = !!rlang::ensym(df1_project_nr) - !!rlang::ensym(df2_project_nr))
    } else {
      adj_factor_df <- adj_factor_df %>%
        mutate(Diff = !!rlang::ensym(df2_project_nr) - !!rlang::ensym(df1_project_nr))
    }
    
    adj_factor_df <- adj_factor_df %>%
      group_by(OlinkID) %>%
      summarise(Adj_factor = if_else(is.na(median(Diff, na.rm = T)),
                                     0,
                                     median(Diff, na.rm = T)))
    
  }else{
    
    #Difference of medians
    
    #Calculate adjustment factors
    adj_factor_df <- df1 %>%
      filter(SampleID %in% overlapping_samples_df1) %>%
      rbind(df2 %>% filter(SampleID %in% overlapping_samples_df2)) %>%
      select(SampleID,OlinkID,UniProt,NPX,Project) %>%
      group_by(Project, OlinkID) %>%
      summarise(Median=median(NPX,na.rm = T)) %>%
      ungroup() %>%
      spread(Project,Median)
    
    if (reference_project == df1_project_nr) {
      adj_factor_df <- adj_factor_df %>%
        mutate(Adj_factor = if_else(is.na(!!rlang::ensym(df1_project_nr) - !!rlang::ensym(df2_project_nr)),
                                    0,
                                    !!rlang::ensym(df1_project_nr) - !!rlang::ensym(df2_project_nr)))
    } else {
      adj_factor_df <- adj_factor_df %>%
        mutate(Adj_factor = if_else(is.na(!!rlang::ensym(df2_project_nr) - !!rlang::ensym(df1_project_nr)),
                                    0,
                                    !!rlang::ensym(df2_project_nr) - !!rlang::ensym(df1_project_nr)))
    }
    
    adj_factor_df <- adj_factor_df %>%
      select(OlinkID, Adj_factor)
  }
  
  
  #Apply adjustment factors
  df_adjusted_data<-df1 %>%
    rbind(df2) %>%
    mutate(Panel=str_replace(Panel,'\\(.+', '')) %>%
    left_join(adj_factor_df,by='OlinkID') %>%
    mutate(Adj_factor = if_else(Project == reference_project,0,Adj_factor)) %>%
    mutate(NPX = NPX + Adj_factor) %>%
    mutate(LOD = LOD + Adj_factor)
  
  
  return(df_adjusted_data)
  
}
#first a function to generalize control labels as control1, control2, ...
# assume that in the original data, the control values are labelled with string containing 'control' characters
#create function
avg_control <- function(df){
  cont.id <- df %>% filter(grepl("control", SampleID, ignore.case = T))
  df <- aggregate(NPX ~ OlinkID+UniProt+Assay+Panel+QC_Warning, cont.id, mean)
  # df <- aggregate(NPX ~ Assay+Index+OlinkID+UniProt+MissingFreq+Panel+Panel_Version+PlateID+QC_Warning+LOD, cont.id, mean)
  df$SampleID <- "AverageControl"
  df <- df %>% mutate(Index=NA, MissingFreq=NA, Panel_Version=NA, PlateID=NA, LOD=NA)
  df <- df[,c(colnames(cont.id))]
  
  return(df)
}



# bootstrap method for normal values
boot_int <- function(data, p1, p2){
  boot_sample <- sample(data, size = length(data), replace = TRUE)
  bsample_sorted <- sort(boot_sample)
  q1 <- p1 * (length(data)+1)
  q2 <- p2 * (length(data)+1)
  lower <- bsample_sorted[round(q1)]
  upper <- bsample_sorted[round(q2)]
  conf <- cbind(lower, upper)
  return(conf)
}

boot_NP <- function(B, data, p1=0.025, p2=0.975, bootfun){
  lower <- rep(0, B)
  upper <- rep(0, B)
  
  bs_rep <- data.frame(lower, upper)
  
  set.seed(79462891)
  for (i in 1:B) {
    bs_rep[i,] <- data.frame(bootfun(sort(data), p1, p2))
  }
  
  #get the interval
  ri <- colMeans(bs_rep)
  return(ri)
}


run_fun <- function(filename1, filename2){
  #### SECOND STEP: perform batch correction across different batches in your data cohort, this in order to ...
  # To define batches Batch 1 should be the biggest batch, or the first batch sampled if there is no difference in size 
  my_files<-list.files("./data", pattern = "\\.xlsx$", full.names=TRUE)
  df<-list()
  for (i in 1:length(my_files)) {
    df.read<-read_NPX(my_files[i])
    df[[i]]<-df.read[, 1:min(12, ncol(df.read))]
  }
  #now each data can be called by df[[1]], df[[2]], ..., so the lines underneeth can be removed
  
  #longdata_b1 <- read_NPX(filename1)
  #longdata_b2 <- read_NPX(filename2)
  #IAF_longformatb1 <- read_NPX(filename3)
  
  #df1 <- longdata_b1
  #df2 <- longdata_b2
  
  # Find overlaping samples between batches to perform the batch correction, but exclude Olink controls, ADD ignore case also for intersect function
  overlap_samples <- intersect((data.frame(df[1]) %>% filter(!grepl("control", SampleID, ignore.case = T)))$SampleID,
                               (data.frame(df[2]) %>% filter(!grepl("control", SampleID, ignore.case= T)))$SampleID)
  #df2$Normalization <- NULL
  
  # code draft adjustmnent for automating the batch correction process 
  #if it doesn't work the way it is now here subsetting DF[1] from the rest and renaming df [2] to end as Df[i] also allows iterating over all the non batch1 things. 
  
  if (length(overlap_samples) > 0) {
    # Perform batch correction only if there are overlapping samples
    batch_corrected_data <- olink_normalization(df1 = data.frame(df[[1]]),
                                                df2 = data.frame(df[[2]]),
                                                overlapping_samples_df1 = overlap_samples,
                                                df1_project_nr = 'P1',
                                                df2_project_nr = 'P2',
                                                reference_project = 'P1')
    
    # Update df with batch corrected data
    df[[2]] <- batch_corrected_data[1:12]
  } else {
    message("No overlapping samples found. Skipping batch correction.")
  }
  
  
  
  #### THIRD STEP: Bridging across different olink studies by using olink external control samples 
  ## still need to automate iteration for all created files in the next 15 lines
  #first to null columning not absolutely crucial 
  #batch_corrected_data$Project <- NULL
  #batch_corrected_data$Adj_factor <- NULL
  #df1 <- batch_corrected_data
  #df2 <- IAF_longformatb1
  
  ### function that applies uniform external control namings
  for (i in 1:length(df)) { 
    
    control_aligned_data <- rbind(avg_control(data.frame(df[i])), data.frame(df[i]))
    df[[i]] <- control_aligned_data	
  }
  

  
  
  #df3 <- rbind(avg_control(df1), df1)
  #df4 <- rbind(avg_control(df2), df2)
  #' # Find overlaping samples, olink controls --> overlap controls are the same controls?
  #' what random means?     
  overlap_samples <- intersect((data.frame(df[2]) %>% filter(grepl("AverageControl", SampleID, ignore.case = T)))$SampleID,
                               (data.frame(df[[length(df)]]) %>% filter(grepl("control", SampleID, ignore.case= T)))$SampleID)
  
  # Function loaded on line  316 will also perform this bridging , but with different type of ...
  # ... overlap samples, no need to load the olink_normalization function again. you can directly run
  
  # for (i in 1:length(df)) { 
  
  
 

    output <- olink_normalization(df1 = data.frame(df[1]),
                                df2 = data.frame(df[[length(df)]]),
                                overlapping_samples_df1 = overlap_samples,
                                df1_project_nr = 'P1',
                                df2_project_nr = 'P2',
                                reference_project = 'P1')
    
    pooleffect_scriptinput <- as.data.frame(read_excel(filename4))
    
   #This is were we recognize whether it is from the old or the new pool
     pooleffect <- function(df) {
      if(any(grepl('SAMPLE_AS', df[[1]]$SampleID))) {
        output <- olink_normalization(df1 = data.frame(df[[length(df)]]),
                                      df2 = data.frame(df[2]),
                                      overlapping_samples_df1 = overlap_samples,
                                      df1_project_nr = 'P1',
                                      df2_project_nr = 'P2',
                                      reference_project = 'P1')
        
        matching_indices <- match(output$Assay, pooleffect_scriptinput$Assay)
        output$NPX_pooleffect <- output$NPX + pooleffect_scriptinput$pooleffect[matching_indices]
        output$pooleffect <- pooleffect_scriptinput$pooleffect[matching_indices]
      } else {
        output <- output  # Or any other operation you wish to perform if the condition is not met
      }
      return(output)
    }

  
  output <- pooleffect(df)
  # df[[i]] <- batch_corrected_data[1:12]
  #}     
  
  # #output <- olink_normalization(df1 = df4,
  #                                 df2 = df3,
  #                                 overlapping_samples_df1 = overlap_samples,
  #                                 df1_project_nr = 'P1',
  #                                 df2_project_nr = 'P2',
  #                                 reference_project = 'P1')
  ########## FOURTH STEP: cosine similarity plot to assess the previous steps and a last control before ...
  # ... normal thresholds can be computed.
  #prepare data for similarity assesment
  # the IAF data is not changed, we brought the NPX values of the external scallop study closer to the IAF data, so change df1 (ext. scallop) values into the ones generated in ...
  # ... the output file
  #df1 <- output[output$Project == 'P2',]
  #d#f1 <- subset(df1, select = -c(Project, Adj_factor)) # can we ensure that the column structure would be the same with other external data?
  #df2 <- df4
  #and continue with the bridged data
  #scallop_data <- df1 %>% filter(grepl("control", SampleID, ignore.case = T)) %>%
  # select(., c(SampleID, NPX, Assay)) %>%
  #pivot_wider(names_from = Assay, values_from = NPX) %>%
  #t(.) 
  
  #colnames(scallop_data) <- scallop_data[1,]
  #scallop_data <- scallop_data[-1,]
  #scallop_data <- matrix(as.numeric(scallop_data),    # Convert to numeric matrix
  #                     ncol = ncol(scallop_data))
  #colnames(scallop_data) <- c(paste0('yourstudy', seq(1:ncol(scallop_data))))
  
  # doing the same for the I AM Frontier (IAF) data (maybe just include the output of this as data for the collaborators)
  # for this example to match with sepsis do only the inflammation pannel cfr line 580 we should do this up front
  #IAF <- df2 %>% filter(grepl("control", SampleID, ignore.case = T)) %>%
  # filter(grepl("inflammation", Panel, ignore.case = T)) %>%
  #select(., c(SampleID, NPX, Assay)) %>%
  #pivot_wider(names_from = Assay, values_from = NPX) %>%
  #t(.) 
  
  #colnames(IAF) <- IAF[1,]
  #IAF <- IAF[-1,]
  #IAF = as.matrix(IAF)
  #IAF <- matrix(as.numeric(IAF),    # Convert to numeric matrix
  #             ncol = ncol(IAF))
  #colnames(IAF) <- c(paste0('IAF', seq(1:ncol(IAF))))
  
  ##still need to merge both numeric matrices and then you can easily compare the both 
  #  cohorts_combined <- cbind(scallop_data, IAF)
  # plotdata <- cosine(cohorts_combined)
  
  
  #final data formatting before you can obtain the plot,
  #plotdata <- reshape2::melt(plotdata)
  
  #for final step transpose the data so cosine measure can be calculated between all control samples as columns
  
  
  #similarity_plot<- ggplot(plotdata, aes(x=Var1, y=Var2, fill=value)) +
  #  
  # geom_tile() +
  
  # scale_fill_gradient2(low = "blue", high = "red", mid = "white",
  
  #                         midpoint = 0.99, limit = c(0.98, 1), space = "Lab",
  
  #                        name = "cosine similarity") +
  
  #  theme_minimal() +
  
  # xlab("") + ylab("") +
  
  #  theme(axis.text.x = element_text(angle = 90, vjust = 1,size = 8, hjust = 1),
  
  #       axis.text.y = element_text(size = 8))
  #  ggsave("similarity_plot.png")
  
  # then we should add another function that gives us the exact similarity scores between the two data
  
  # FINAL STEP :calculate population normal value for each assay in corrected data P2
  df1 <- output[output$Project == 'P2',] %>% 
    filter(!grepl("control", SampleID, ignore.case = T))
  
  df.wide <- dcast(df1, SampleID ~ UniProt, 
                   value.var = "NPX", fun.aggregate = mean)


  int <- matrix(0, nrow = ncol(df.wide)-1, ncol=2)
  normal_intervals <- cbind.data.frame(colnames(df.wide[,-1]), int)
  colnames(normal_intervals) <- c('uniprot', 'lower', 'upper')
  for (i in 1:nrow(normal_intervals)) {
    np <- boot_NP(B=1000, data = df.wide[, i+1], p1=0.025, p2=0.975, bootfun = boot_int)
    normal_intervals[i, 2] <- np[1]
    normal_intervals[i, 3] <- np[2]
  }
  save(normal_intervals, file = "normal_intervals.RData")
  
  #plotdata <- as.data.frame(plotdata)
  outputlist <- list(normal_intervals=normal_intervals)
  write.csv(output,"./data/output/cohort_correction.csv")
  write.csv(normal_intervals,"./data/output/intervals.csv")
  outputlist
}
