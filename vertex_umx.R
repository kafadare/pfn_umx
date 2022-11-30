# set libPaths & load pkgs
.libPaths(c( "/gpfs/fs001/cbica/projects/bgd-pfn/.miniconda3/envs/umx-4.15/lib/R/library", "/gpfs/fs001/cbica/projects/bgd-pfn/R/x86_64-conda-linux-gnu-library/4.1", "/usr/share/R/library" ,  "/usr/lib64/R/library" , "/gpfs/fs001/cbica/software/external/R" ))
require(umx)
require(tidyverse)
require(magrittr)
require(data.table)
require(stringr)
require(dplyr)
#require(ciftiTools)

#setwd('/cbica/projects/bgd-pfn/vertex_results')
#
#
#Functions
vertex_umx_fct <- function(pfn_no, vertex_no) {
  source("/cbica/projects/bgd-pfn/pfn_umx/all_fcts.R")
  # assign custom fcts to the umx environment to access hidden umx fcts
  environment(umxSummaryACE_EK) <- asNamespace("umx")
  assignInNamespace("umxSummary", umxSummaryACE_EK, ns = "umx")
  assignInNamespace("umxSummaryACE", umxSummaryACE_EK, ns = "umx")
  
  environment(umxReduceACE_EK) <- asNamespace("umx")
  assignInNamespace("umxReduceACE", umxReduceACE_EK, ns = "umx")
  
  environment(umxCompare_EK) <- asNamespace("umx")
  assignInNamespace("umxCompare", umxCompare_EK, ns = "umx")
  
  # Getting pathname for pfn/vertex column
  path_template1 <- "/cbica/projects/bgd-pfn/vertex_columns/PFN"
  path_template2 <- "V"
  path <- paste0(path_template1, pfn_no, path_template2, vertex_no, ".csv")
  col_v <- read.csv(path, header = T, sep = ",")
  name_vertex <- names(col_v)[2]
  print(name_vertex)
  
  ##defining vectors/templates
  ABCD_twins_pfn <- read.csv('/cbica/projects/bgd-pfn/ABCD_twins_PFN.csv')
  #ABCD_twins_pfn <- read.csv("/Users/ekafadar/Documents/ABCD_pfn_umx_update/ABCD_twins_PFN.csv")
  #PFN stuff to help with analysis later
  names_pfn_17 <- rep('PFN', 17) %>% paste(seq(1:17), sep = "")
  names_pfn_7 <- c('VIS', 'MOT', 'DAN', 'VAN', 'FPN', 'AUD', 'DMN')
  names_pfn_all <- c(names_pfn_17,names_pfn_7)
  names_all <- c(names_pfn_all, "varexp")
  
  #define some vectors for later
  unimodal <- c("PFN6", "PFN10", "PFN2", "PFN4", "PFN11", "PFN13", "PFN16", "VIS", "MOT", "AUD")
  heteromodal <- c("PFN12", "PFN8","PFN1", "PFN17", "PFN15", "PFN3", "PFN9", "PFN7", "PFN14", "PFN5", "PFN16", "PFN17", "FPN", "DAN", "VAN", "DMN")
  pfn_name_color <- data.table(PFN = c(paste("PFN", 1:17, sep = "")), Name = c("DMN","MOT", "FPN", "MOT", "DAN", "VIS", "VAN", "DMN", "VAN", "VIS", "MOT", "DMN", "MOT", "DAN", "FPN", "AUD", "FPN"), Color = c("#CB94A9", "#AEBACA", "#DFCA8B", "#7793BB", "#83AD73", "#785397", "#B97FE7", "#AF576F","#A241EA","#5A2C7B","#5E78AF", "#833852", "#445489", "#517549", "#BD986A", "#5647A1", "#AA7C45"))
  pfn_name_color$Modality <- ifelse(pfn_name_color$PFN %in% unimodal, "UNI", ifelse(pfn_name_color$PFN %in% heteromodal, "HETERO", NA))
  pfn7 <- data.table(PFN = c("VIS", "MOT", "DAN", "VAN", "FPN", "AUD", "DMN"), Name = c("VIS", "MOT", "DAN", "VAN", "FPN", "AUD", "DMN"), Color = rep(NA, 7), Modality = c("UNI", "UNI", "HETERO", "HETERO", "HETERO", "UNI", "HETERO"))
  pfn_name_color = rbind(pfn_name_color, pfn7)
  
  SA_Axis <- c(6,10,2,13,11,4,5,16,14,15,9,8,12,3,1,7,17)
  SA_Axis <- paste("PFN", SA_Axis, sep = "")
  
  #look at variance in network size & plot
  
  cv_pfns <- data.frame(matrix(ncol = 2, nrow = length(names_all) ) ) %>% setNames(c('PFN', 'coef_var'))
  
  for(i in seq(1:length(names_all))) {
    var <- names_all[i] #variable names to select
    cv_pfns[i,1] <- var
    cv_pfns[i,2] <- sd(ABCD_twins_pfn[,var]) / mean(ABCD_twins_pfn[,var])* 100
  }
  
  #load data & shape accordingly
  #ABCD_twins <- read.table('ABCD_release_2.0.1_r1.twin_table.20201109.txt', header = TRUE, sep = "", dec = ".")
  ABCD_twins_pfn$numTRs <- as.numeric(ABCD_twins_pfn$numTRs)
  ABCD_twins_pfn$sex <- as.factor(ABCD_twins_pfn$sex)# code F as 1 M as 0
  
  
  #reshape into wide format
  ABCD_twins_vertex <- merge(ABCD_twins_pfn, col_v, by = c('subjectkey'))
  
  #reshape twins data back to paired
  
  var_names_long <- names(ABCD_twins_vertex)
  var_names_long <- var_names_long[-which(var_names_long %in% c("famID", "zygosity", "kinship", "X", "zyg"))]
  
  twins_vertex_wide <- umx_long2wide(ABCD_twins_vertex, famID = "famID", twinID = "twinno", vars2keep = var_names_long, zygosity = "zyg", passalong = c('zygosity', 'kinship', 'twinno'))
  
  twins_vertex_wide <- subset(twins_vertex_wide, !is.na(twins_vertex_wide$PFN1_T1) & !is.na(twins_vertex_wide$PFN1_T2)) #remove if one twin is missing PFN data
  
  rm(var_names_long)
  
  #scale trait variable as well as numerical covars :)
  scale_vars <- c('age', 'meanFD', 'numTRs')
  twins_vertex_wide <- umx_scale_wide_twin_data(varsToScale = c(name_vertex, scale_vars), sep = "_T", twins_vertex_wide, twins = 1:2)
  variance <- var(twins_vertex_wide[, paste0(name_vertex, "_T1")])
  selCovs <- c('age', 'meanFD', 'numTRs', 'sex')
  if (variance != 0) {
    dz <- twins_vertex_wide[twins_vertex_wide$zygosity == "DZ",]
    mz <- twins_vertex_wide[twins_vertex_wide$zygosity == "MZ",]
    ace <- umx_modelCompare(ROI_all = name_vertex, mz = mz, dz = dz, sep = "_T", tryHard = "no", selCovs = selCovs)
    # ace <- umxACE(selDVs = name_vertex, mzData = mz, dzData = dz,sep = "_T", tryHard = "no", selCovs = selCovs )
    ace_Vdf <- umx_df(ace)
    ace_Vdf$vertex_no <- vertex_no
    ace_Vdf$pfn_no <- pfn_no
    fName <- paste0("/cbica/projects/bgd-pfn/vertex_results/",name_vertex, "_ace_df.csv" )
    write.csv(ace_Vdf, file = fName, row.names = F, col.names = T)
  } else{
    fName_txt <- paste0(name_vertex, "_NoWork.txt" )
    error <- paste0("Variance of", name_vertex, "after scaling is ", variance, " . This is < 0.1, so we're NOT running umx, aborting script.")
    writeLines(error, fName_txt)
    stop(paste0("Variance of phenotype after scaling is ", variance, " . This is < 0.1, so we're running NOT umx, aborting script."), call. = FALSE)
  }
  return(ace_Vdf)
  unlink(path)
}

# Get the arguments as a character vector.
if (exists("pfn_no") & exists("vertex_no")) {
  ace_Vdf = vertex_umx_fct(pfn_no, vertex_no)
} else if (length(commandArgs(trailingOnly=TRUE)) == 2) {
  myargs = commandArgs(trailingOnly=TRUE)
  #index = as.numeric(myargs[1])
  #args <- as.numeric(read.csv("pfn1_vertex_pairs.csv", skip = (index-1), nrow = 1))
  pfn_no <- as.numeric(myargs[1])
  vertex_no <- as.numeric(myargs[2])
  ace_Vdf = vertex_umx_fct(pfn_no, vertex_no)
} else {
  myargs = commandArgs(trailingOnly=TRUE)
  stop(" Usage: vertex_umx.R <alignment.fasta> <tree.nex> <output.txt>", call.=FALSE)
}
cat("You ran the program with ", myargs,"\n")


