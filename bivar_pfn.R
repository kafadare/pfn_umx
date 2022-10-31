# set libPaths & load pkgs
.libPaths(c( "/gpfs/fs001/cbica/projects/bgd-pfn/.miniconda3/envs/umx-4.15/lib/R/library", "/gpfs/fs001/cbica/projects/bgd-pfn/R/x86_64-conda-linux-gnu-library/4.1", "/usr/share/R/library" ,  "/usr/lib64/R/library" , "/gpfs/fs001/cbica/software/external/R" ))
require(umx)
require(tidyverse)
require(magrittr)
require(data.table)
require(stringr)
require(dplyr)
#require(ciftiTools)

source("all_fcts.R")

# assign custom fcts to the umx environment to access hidden umx fcts
environment(umxSummaryACE_EK) <- asNamespace("umx")
assignInNamespace("umxSummary", umxSummaryACE_EK, ns = "umx")
assignInNamespace("umxSummaryACE", umxSummaryACE_EK, ns = "umx")

environment(umxReduceACE_EK) <- asNamespace("umx")
assignInNamespace("umxReduceACE", umxReduceACE_EK, ns = "umx")

environment(umxCompare_EK) <- asNamespace("umx")
assignInNamespace("umxCompare", umxCompare_EK, ns = "umx")

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



#reshape twins data back to paired
var_names_long <- names(ABCD_twins_pfn) 
var_names_long <- var_names_long[-which(var_names_long %in% c("famID", "zygosity", "kinship", "X", "zyg"))]

twins_pfn_wide <- umx_long2wide(ABCD_twins_pfn, famID = "famID", twinID = "twinno", vars2keep = var_names_long, zygosity = "zyg", passalong = c('zygosity', 'kinship', 'twinno'))

twins_pfn_wide <- subset(twins_pfn_wide, !is.na(twins_pfn_wide$PFN1_T1) & !is.na(twins_pfn_wide$PFN1_T2)) #remove if one twin is missing PFN data

rm(var_names_long)
scale_vars <- c('age', 'meanFD', 'numTRs')
#scale varexp to increase variance for umx to run  better. currently variance is < 0.1. This fct z-scores the data.
twins_pfn_wide <- umx_scale_wide_twin_data(varsToScale = c(names_all, scale_vars), sep = "_T", twins_pfn_wide, twins = 1:2)
dz = twins_pfn_wide[twins_pfn_wide$zygosity == "DZ",]
mz = twins_pfn_wide[twins_pfn_wide$zygosity == "MZ",]

#get pfn pair from args
# Get the arguments as a character vector.
myargs = commandArgs(trailingOnly=TRUE)
cat(myargs, sep = "\n")
index = as.numeric(myargs[1])
# The first argument is the line of csv file to read.
cat(index)
pfn_nos <- as.numeric(read.csv("bivar_input.csv", skip = (index-1), nrow = 1))
if (pfn_nos[1] == pfn_nos[2]){
stop ("Same pfn with itself, bivar won't work. Skipping this.")
}
#pfn_nos = c(myargs[1], myargs[2])
pairs_pfn <- paste0("PFN", pfn_nos)
print(pairs_pfn)
var_name <- paste(pairs_pfn, collapse = '')
var_name_df <- paste(var_name, "df", sep = "_")
#run bivar
umx_set_optimizer("SLSQP")
selCovs <- c('age', 'meanFD', 'numTRs', 'sex')
shared_env <- c("c_r1c1", "c_r2c1", "c_r2c2")
add_gen <- c("a_r1c1", "a_r2c1", "a_r2c2")
bivar_model <- umx_modelCompare( ROI_all = pairs_pfn, mz = mz, dz = dz, selCovs = selCovs, tryHard = "no", sep = "_T")
assign(var_name, bivar_model)
bivar_df <- umx_df(bivar_model)
assign(var_name_df, bivar_df)

fName <- paste0("/cbica/projects/bgd-pfn/bivar_PFN_results/",var_name, "_ace")
#fName <- paste0("/Users/ekafadar/Documents/pfn_scripts/",var_name, "_ace")

#save files
write.csv(get(var_name_df), file = paste0(fName, "_df.csv"), row.names = F, col.names = T)
save(var_name, file = paste(fName,format(Sys.time(), "%m_%d_%y"),"RData", sep = "."))
