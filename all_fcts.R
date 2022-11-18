#Functions
## Functions - analysis
#quick little function to calculate correlations, Falconer estimates and means.
#input: mz and dz data frames, ROI name. 
require(umx)
falconer<-function(mz,dz,selDVs,sep,plot=F){
  
  corMZ<-cor(mz[,paste(selDVs,1,sep=sep)],mz[,paste(selDVs,2,sep=sep)],use="pair")
  corDZ<-cor(dz[,paste(selDVs,1,sep=sep)],dz[,paste(selDVs,2,sep=sep)],use="pair")
  
  MZm<-mean(c(mz[,paste(selDVs,1,sep=sep)],mz[,paste(selDVs,2,sep=sep)]),na.rm=T)
  DZm<-mean(c(dz[,paste(selDVs,1,sep=sep)],dz[,paste(selDVs,2,sep=sep)]),na.rm=T)
  
  a2<-2*(corMZ-corDZ)
  c2<-(2*corDZ)-corMZ
  e2<-(1-a2-c2)
  
  if(plot==T){
    r<-range(c(mz[,paste(selDVs,1,sep=sep)],mz[,paste(selDVs,2,sep=sep)],dz[,paste(selDVs,1,sep=sep)],dz[,paste(selDVs,2,sep=sep)]))
    plot(NA,NA,xlim=r,ylim=r,xlab=paste(selDVs,1),ylab=paste(selDVs,2))
    points(mz[,paste(selDVs,1,sep=sep)],mz[,paste(selDVs,2,sep=sep)],col="red",cex=.5,pch=15)
    points(dz[,paste(selDVs,1,sep=sep)],dz[,paste(selDVs,2,sep=sep)],col="blue",cex=.5,pch=16)
    legend("topright",legend=c("MZ","DZ"),col=c("red","blue"),bg="white",pch=c(15,16))
    legend("topleft",legend=c(paste("cor MZ =",round(corMZ,digits=2)),
                              paste("cor DZ =",round(corDZ,digits=2)),
                              paste("a2 =",round(a2,digits=2)),
                              paste("c2 =",round(c2,digits=2)),
                              paste("e2 =",round(e2,digits=2))),
           bg="white")
    
  }
  
  out<-c(corMZ,corDZ,a2,c2,e2,MZm,DZm)
  out<-data.frame(t(unlist(out)))
  names(out)<-c("corMZ","corDZ","a2","c2","e2","MZm","DZm")
  return(out)
}

#to get coefficients from umx model output
acevec<-function(mo){
  
  out<-c(
    invisible(umxSummary(mo)$a1),
    invisible(umxSummary(mo)$c1),
    invisible(umxSummary(mo)$e1),
    invisible(umxSummary(mo)$a1**2),
    invisible(umxSummary(mo)$c1**2),
    invisible(umxSummary(mo)$e1**2)
  )
  out<-data.frame(t(unlist(out)))
  names(out)<-c("a1","c1","e1","a-2","c-2","e-2")
  return(out)
}

aevec<-function(mo){
  
  out<-c(
    invisible(umxSummary(mo)$a1),
    invisible(umxSummary(mo)$e1),
    invisible(umxSummary(mo)$a1**2),
    invisible(umxSummary(mo)$e1**2)
  )
  out<-data.frame(t(unlist(out)))
  names(out)<-c("a1","e1","a-2","e-2")
  return(out)
}

#function for bivar analysis and formatting output
umx_modelCompare <- function(ROI_all, selCovs = NULL , tryHard = c("no","yes"), dz, mz, sep = "_T", sex = F, mzm = NULL, mzf = NULL, dzm = NULL, dzf = NULL, dzo = NULL, optimizer = NULL){
  selDVs<- c(ROI_all) # two ROIs to run together
  if (length(ROI_all) == 1) {
    falconer <- falconer(mz, dz, selDVs = ROI_all, sep = sep)
    showRg = FALSE
    out <- list()
    out.names <- c("ROIs", "falconer", "summary_ACE", "reduced_ACE","subModels_ACE", "summary_AE",  "submodels_AE")
    out$falconer <- falconer
  }else{
    showRg = TRUE
    out <- list()
    out.names <- c("ROIs", "summary_ACE", "reduced_ACE","subModels_ACE", "summary_AE", "submodels_AE")
  }
  shared_env <- c("c_r1c1", "c_r2c1", "c_r2c2")
  add_gen <- c("a_r1c1", "a_r2c1", "a_r2c2")
  #ACE model
  if (sex){
    #haven't tested this yet, EK 08/09/22
    ACE<- umx::umxSexLim(selDVs = selDVs, selCovs = selCovs, mzmData = mzm, dzmData = dzm, mzfData = mzf, dzfData = dzf, sep = sep,   tryHard = tryHard)
    #might need special summary fct for SexLim, this is not tested yet 08/09/22
    summary_sexACE <- umx::umxSummaryACE(ACE, showRg = showRg, sex = T, list_out = T, CIs = T)
  } else{
    ACE <- umx::umxACE(selDVs = selDVs, selCovs = selCovs, dzData = dz, mzData = mz, sep = sep, tryHard = tryHard, intervals = T, optimizer = optimizer)
    status <- ACE$output$status$code
    status_TH <- NA
    if (grepl("[2-9]", status)){
      tryHard = "yes"
      ACE <- umx::umxACE(selDVs = selDVs, selCovs = selCovs, dzData = dz, mzData = mz, sep = sep, tryHard = tryHard, intervals = T, optimizer = optimizer)
      status_TH = ACE$output$status$code
    }
    summary_ACE <- umx::umxSummaryACE(ACE, showRg = showRg, CIs = T, list_out = T, extended = T)
    summary_ACE$status_code <- status
    summary_ACE$status_code_TH <- status_TH
    
    #ACE full submodel analysis
    reduce_ACE <- umx::umxReduceACE(ACE, list_out = T, tryHard = tryHard, intervals = T, optimizer = optimizer)
  }
  #compare AE to E
  AE <- umxModify(ACE, regex = "c_r[0-9]+c[0-9]+", name = "AE", tryHard = tryHard, intervals = T)
  AE <- umxConfint(AE, parm = c("existing"), run = T) #drop shared env
  E <- umxModify(AE,regex = "a_r[0-9]+c[0-9]+", name = "E", tryHard = tryHard, intervals = T) #drop additive genetics
  
  sub_ACE <- umx::umxCompare(ACE,AE, report = "markdown", list_out = T, compareWeightedAIC = T) #compare ACE and AE models
  sub_AE <- umx::umxCompare(AE,E, report = "markdown", list_out = T, compareWeightedAIC = T) #compare AE and E models
  summary_AE <- umx::umxSummaryACE(AE, showRg = showRg, CIs = T, list_out = T, extended = T) #get genetic correlations & path loadings from AE model
  
  ROIs  <- paste(ROI_all, sep = "_")	
  out$ROIS <- ROIs
  out$summary_ACE <- summary_ACE
  out$reducedACE <- reduce_ACE
  out$summary_AE <- summary_AE
  out$subModelsACE <- sub_ACE
  out$subModelsAE <- sub_AE        
  return(out)
}  


#should customize this fct to work w/ multivar != 2, a1,a2,c1,c2 stuff should vary based of # of vars in ACE model
#use regexp to do this , it isn't hard LOL
umx_df <- function(umx_output){
  
  rois <- as.data.frame(strsplit(umx_output$ROIS, "_"))
  n <- c(1:dim(rois)[2])
  roi_no <- paste("ROI",n, sep = "_")
  names(rois) <- roi_no
  path_names_ACE <- names(umx_output$summary_ACE$Estimates)
  path_names_AE <- names(umx_output$summary_AE$Estimates)
  CI_names <- colnames(umx_output$summary_ACE$model$output$confidenceIntervals)
  if (dim(rois)[2] == 1) {
    a1ACE<-as.data.frame(unlist(umx_output$summary_ACE$Estimates[,which(names(umx_output$summary_ACE$Estimates)=="a1")]))
    names(a1ACE) <-paste ("a1_ACE",roi_no, sep = "_")
    c1ACE<-as.data.frame(unlist(umx_output$summary_ACE$Estimates[,which(names(umx_output$summary_ACE$Estimates)=="c1")]))
    names(c1ACE)<-paste("c1_ACE",roi_no,sep = "_")
    e1ACE<-as.data.frame(unlist(umx_output$summary_ACE$Estimates[,which(names(umx_output$summary_ACE$Estimates)=="e1")]))
    names(e1ACE)<-paste("e1_ACE",roi_no,sep = "_")
    ntest <- paste (sort(rep(path_names_ACE,length(roi_no))),roi_no,sep = "_")
    #a1ACE<- as.data.frame((unlist(umx_output$summary_ACE$model$output$confidenceIntervals[grepl("top.a_std", row.names(umx_output$summary_ACE$model$output$confidenceIntervals)),])))
    
    
    # names(a1ACE)<-paste("a1_ACE",CI_names,sep = "_")
    
    # paths <- names(test$summary_ACE$Estimates)
    # CI_names <- paste(sort(rep(paths,2)), c("upper", "lower"), sep = "_")
    # CIs <- as.data.frame(unlist(umx_output$summary_ACE$CIs)))
    # names(CIs) <- paste(paths, )
    
    a1AE<-as.data.frame(unlist(umx_output$summary_AE$Estimates[,which(names(umx_output$summary_AE$Estimates)=="a1")]))
    names(a1AE)<-paste("a1_AE",roi_no,sep = "_")
    e1AE<-as.data.frame(unlist(umx_output$summary_AE$Estimates[,which(names(umx_output$summary_AE$Estimates)=="e1")]))
    names(e1AE)<-paste("e1_AE",roi_no,sep = "_")
    
    
    falconer <- umx_output$falconer
  }
  else{
    a1ACE<-as.data.frame(t(unlist(umx_output$summary_ACE$Estimates[,which(names(umx_output$summary_ACE$Estimates)=="a1")])))
    names(a1ACE) <-paste ("a1_ACE",roi_no, sep = "_")
    a2ACE<-as.data.frame(t(unlist(umx_output$summary_ACE$Estimates[,which(names(umx_output$summary_ACE$Estimates)=="a2")])))
    names(a2ACE)<-paste("a2_ACE",roi_no,sep = "_")
    c1ACE<-as.data.frame(t(unlist(umx_output$summary_ACE$Estimates[,which(names(umx_output$summary_ACE$Estimates)=="c1")])))
    names(c1ACE)<-paste("c1_ACE",roi_no,sep = "_")
    c2ACE<-as.data.frame(t(unlist(umx_output$summary_ACE$Estimates[,which(names(umx_output$summary_ACE$Estimates)=="c2")])))
    names(c2ACE)<-paste("c2_ACE",roi_no,sep = "_")
    e1ACE<-as.data.frame(t(unlist(umx_output$summary_ACE$Estimates[,which(names(umx_output$summary_ACE$Estimates)=="e1")])))
    names(e1ACE)<-paste("e1_ACE",roi_no,sep = "_")
    e2ACE<-as.data.frame(t(unlist(umx_output$summary_ACE$Estimates[,which(names(umx_output$summary_ACE$Estimates)=="e2")])))
    names(e2ACE)<-paste("e2_ACE",roi_no,sep = "_")
    
    # paths <- names(test$summary_ACE$Estimates)
    # CI_names <- paste(sort(rep(paths,2)), c("upper", "lower"), sep = "_")
    # CIs <- as.data.frame(t(unlist(umx_output$summary_ACE$CIs)))
    # names(CIs) <- paste(paths, )
    
    rA1ACE <- as.data.frame(t(unlist(umx_output$summary_ACE$Correlations[,which(names(umx_output$summary_ACE$Correlations)=="rA1")])))
    names(rA1ACE)<-paste("rA1_ACE",roi_no, sep="_")
    rA2ACE <- as.data.frame(t(unlist(umx_output$summary_ACE$Correlations[,which(names(umx_output$summary_ACE$Correlations)=="rA2")])))
    names(rA2ACE)<-paste("rA2_ACE",roi_no,sep = "_")
    rC1ACE <- as.data.frame(t(unlist(umx_output$summary_ACE$Correlations[,which(names(umx_output$summary_ACE$Correlations)=="rC1")])))
    names(rC1ACE)<-paste("rC1_ACE",roi_no,sep = "_")
    rC2ACE <- as.data.frame(t(unlist(umx_output$summary_ACE$Correlations[,which(names(umx_output$summary_ACE$Correlations)=="rC2")])))
    names(rC2ACE)<-paste("rC2_ACE",roi_no,sep = "_")
    rE1ACE <- as.data.frame(t(unlist(umx_output$summary_ACE$Correlations[,which(names(umx_output$summary_ACE$Correlations)=="rE1")])))
    names(rE1ACE)<-paste("rE1_ACE",roi_no,sep = "_")
    rE2ACE <- as.data.frame(t(unlist(umx_output$summary_ACE$Correlations[,which(names(umx_output$summary_ACE$Correlations)=="rE2")])))
    names(rE2ACE)<-paste("rE2_ACE",roi_no,sep = "_")
    
    a1AE<-as.data.frame(t(unlist(umx_output$summary_AE$Estimates[,which(names(umx_output$summary_AE$Estimates)=="a1")])))
    names(a1AE)<-paste("a1_AE",roi_no,sep = "_")
    a2AE<-as.data.frame(t(unlist(umx_output$summary_AE$Estimates[,which(names(umx_output$summary_AE$Estimates)=="a2")])))
    names(a2AE)<-paste("a2_AE",roi_no,sep = "_")
    e1AE<-as.data.frame(t(unlist(umx_output$summary_AE$Estimates[,which(names(umx_output$summary_AE$Estimates)=="e1")])))
    names(e1AE)<-paste("e1_AE",roi_no,sep = "_")
    e2AE<-as.data.frame(t(unlist(umx_output$summary_AE$Estimates[,which(names(umx_output$summary_AE$Estimates)=="e2")])))
    names(e2AE)<-paste("e2_AE",roi_no,sep = "_")
    
    rA1AE <- as.data.frame(t(unlist(umx_output$summary_AE$Correlations[,which(names(umx_output$summary_AE$Correlations)=="rA1")])))
    names(rA1AE)<-paste("rA1_AE",roi_no,sep = "_")
    rA2AE <- as.data.frame(t(unlist(umx_output$summary_AE$Correlations[,which(names(umx_output$summary_AE$Correlations)=="rA2")])))
    names(rA2AE)<-paste("rA2_AE",roi_no,sep = "_")
    rE1AE <- as.data.frame(t(unlist(umx_output$summary_AE$Correlations[,which(names(umx_output$summary_AE$Correlations)=="rE1")])))
    names(rE1AE)<-paste("rE1_AE",roi_no,sep = "_")
    rE2AE <- as.data.frame(t(unlist(umx_output$summary_AE$Correlations[,which(names(umx_output$summary_AE$Correlations)=="rE2")])))
    names(rE2AE)<-paste("rE2_AE",roi_no,sep = "_")
  }   
  
  modname_ACE <- umx_output$reducedACE$modelCompare$Model	 
  
  AIC_ACE<-as.data.frame(t(unlist(umx_output$reducedACE$modelCompare[,which(names(umx_output$reducedACE$modelCompare)=="AIC")])))
  names(AIC_ACE)<-paste("AIC_ACE",modname_ACE,sep="_")
  dAIC_ACE<-as.data.frame(t(unlist(umx_output$reducedACE$modelCompare[,which(names(umx_output$reducedACE$modelCompare)=="Δ AIC")])))
  names(dAIC_ACE)<-paste("dAIC_ACE",modname_ACE,sep="_")
  X2_ACE<-as.data.frame(t(unlist(umx_output$reducedACE$modelCompare[,which(names(umx_output$reducedACE$modelCompare)=="Δ Fit")])))
  names(X2_ACE)<-paste("X2_ACE",modname_ACE,sep="_")
  dDF_ACE<- as.data.frame(t(unlist(umx_output$reducedACE$modelCompare[,which(names(umx_output$reducedACE$modelCompare)=="Δ df")])))
  names(dDF_ACE)<-paste("dDF_ACE", modname_ACE,sep="_")  	   		
  p_ACE<-as.data.frame(t(unlist(umx_output$reducedACE$modelCompare[,which(names(umx_output$reducedACE$modelCompare)=="p")])))
  names(p_ACE)<-paste("ptim_ACE",modname_ACE,sep="_")
  
  
  modname_AE <- umx_output$subModelsAE$tablePub$Model
  
  AIC_AE<-as.data.frame(t(unlist(umx_output$subModelsAE$tablePub[,which(names(umx_output$subModelsAE$tablePub)=="AIC")])))
  names(AIC_AE)<-paste("AIC_AE",modname_AE,sep="_")
  dAIC_AE<-as.data.frame(t(unlist(umx_output$subModelsAE$tablePub[,which(names(umx_output$subModelsAE$tablePub)=="Δ AIC")])))
  names(dAIC_AE)<-paste("dAIC_AE",modname_AE,sep="_")
  X2_AE<-as.data.frame(t(unlist(umx_output$subModelsAE$tablePub[,which(names(umx_output$subModelsAE$tablePub)=="Δ Fit")])))
  names(X2_AE)<-paste("X2_AE", modname_AE,sep="_")
  dDF_AE<- as.data.frame(t(unlist(umx_output$subModelsAE$tablePub[,which(names(umx_output$subModelsAE$tablePub)=="Δ df")])))
  names(dDF_AE)<-paste("dDF_AE", modname_AE,sep="_")
  p_AE<-as.data.frame(t(unlist(umx_output$subModelsAE$tablePub[,which(names(umx_output$subModelsAE$tablePub)=="p")])))
  names(p_AE)<-paste("ptim_AE",modname_AE,sep="_")
  
  modname_AEvACE <- umx_output$subModelsACE$tablePub$Model
  
  AIC_AEvACE<-as.data.frame(t(unlist(umx_output$subModelsACE$tablePub[,which(names(umx_output$subModelsACE$tablePub)=="AIC")])))
  names(AIC_AEvACE)<-paste("AIC_AEvACE",modname_AEvACE,sep="_")
  dAIC_AEvACE<-as.data.frame(t(unlist(umx_output$subModelsACE$tablePub[,which(names(umx_output$subModelsACE$tablePub)=="Δ AIC")])))
  names(dAIC_AEvACE)<-paste("dAIC_AEvACE",modname_AEvACE,sep="_")
  X2_AEvACE<-as.data.frame(t(unlist(umx_output$subModelsACE$tablePub[,which(names(umx_output$subModelsACE$tablePub)=="Δ Fit")])))
  names(X2_AEvACE)<-paste("X2_AEvACE", modname_AEvACE,sep="_")
  dDF_AEvACE<- as.data.frame(t(unlist(umx_output$subModelsACE$tablePub[,which(names(umx_output$subModelsACE$tablePub)=="Δ df")])))
  names(dDF_AEvACE)<-paste("dDF_AEvACE", modname_AEvACE,sep="_")
  p_AEvACE<-as.data.frame(t(unlist(umx_output$subModelsACE$tablePub[,which(names(umx_output$subModelsACE$tablePub)=="p")])))
  names(p_AEvACE)<-paste("ptim_AEvACE",modname_AEvACE,sep="_")
  
  bestModels <- c(umx_output$reducedACE$bestModel$name, umx_output$subModelsACE$bestModel$name,   umx_output$subModelsAE$bestModel$name)
  names(bestModels) <- c("bestModel_ACE_rd", "bestModel_ACE", "bestModel_AE")
  status_code_ACE <- umx_output$summary_ACE$status_code
  status_code_ACE_TH <- umx_output$summary_ACE$status_code_TH
  names(status_code_ACE) <- "status_code_ACE"
  names(status_code_ACE_TH) <- "status_code_TH_ACE"
  
  if(dim(rois)[2] == 1){
    df <- as.data.frame(t(unlist(c(rois,AIC_ACE,dAIC_ACE,X2_ACE,dDF_ACE, p_ACE, bestModels, status_code_ACE, a1ACE, c1ACE,  e1ACE,  AIC_AE, dAIC_AE, X2_AE, dDF_AE, p_AE, AIC_AEvACE, dAIC_AEvACE, X2_AEvACE, dDF_AEvACE, p_AEvACE, a1AE, e1AE, falconer))))	
  }else{
    df <- as.data.frame(t(unlist(c(rois,AIC_ACE,dAIC_ACE,X2_ACE,dDF_ACE, p_ACE, bestModels, status_code_ACE, a1ACE, a2ACE, c1ACE, c2ACE, e1ACE, e2ACE, rA1ACE, rA2ACE, rC1ACE, rC2ACE, rE1ACE, rE2ACE, AIC_AE, dAIC_AE, X2_AE, dDF_AE, p_AE, AIC_AEvACE, dAIC_AEvACE, X2_AEvACE, dDF_AEvACE, p_AEvACE, a1AE, a2AE, e1AE, e2AE, rA1AE, rA2AE, rE1AE, rE2AE))))	
  } 
  
  #remove < signs from p-value columns for easier quant analysis
  cols_ptim <- names(df)[str_which(names(df), "ptim")]
  df<- df %>% mutate_at(vars(contains("ptim")),.funs = list(clean = ~gsub("< ", "", .))) %>% select(!cols_ptim)
  names(df) <- sub("_clean", "",names(df))
  
  # make all of them numeric
  is_num <- names(df)[names(df) %!in% c("ROI_1","ROI_2", "bestModel_ACE_rd", "bestModel_ACE", "bestModel_AE")]
  df<- df %>% mutate_at(is_num,.funs = list(num = ~as.numeric(.))) %>% select(!is_num)
  names(df) <- sub("_num", "",names(df))
  
  #convert path coefs to variances (square them)
  df<- df %>% mutate_at(vars(contains("_ROI")),.funs = list(sq = ~.^2)) %>%
    rename_at(vars(contains( "_sq") ), list( ~paste("SQ", gsub(c("_sq"), "", .), sep = "_") ) )
  
  if(dim(rois)[2] == 1){
    #calc heritability
    df <- df %>% mutate(h2_ACE = SQ_a1_ACE_ROI_1/(SQ_a1_ACE_ROI_1+ SQ_c1_ACE_ROI_1 + SQ_e1_ACE_ROI_1), h2_AE = SQ_a1_AE_ROI_1/(SQ_a1_AE_ROI_1 + SQ_e1_AE_ROI_1))
  }
  else{
    #calc heritability
    df <- df %>% mutate(
      h2_ACE_ROI_1 = SQ_a1_ACE_ROI_1 / (SQ_a1_ACE_ROI_1+ SQ_c1_ACE_ROI_1 + SQ_e1_ACE_ROI_1), 
      h2_ACE_ROI_2 = (SQ_a1_ACE_ROI_2 + SQ_a2_ACE_ROI_2) / (SQ_a1_ACE_ROI_2 + SQ_a2_ACE_ROI_2 + SQ_c1_ACE_ROI_2 +   SQ_c2_ACE_ROI_2 + SQ_e1_ACE_ROI_2 + SQ_e2_ACE_ROI_2),
      
      h2_AE_ROI_1 = SQ_a1_AE_ROI_1 / (SQ_a1_AE_ROI_1 + SQ_e1_AE_ROI_1), 
      h2_AE_ROI_2 = (SQ_a1_AE_ROI_2 + SQ_a2_AE_ROI_2) / (SQ_a1_AE_ROI_2 + SQ_a2_AE_ROI_2 + SQ_e1_AE_ROI_2 + SQ_e2_AE_ROI_2))
    
    #calculate additive genetic correlation
    df <- df %>% mutate( 
      a_corr_AE = rA1_AE_ROI_2 * sqrt(h2_AE_ROI_1)*sqrt(h2_AE_ROI_2),
      a_corr_ACE = rA1_ACE_ROI_2 * sqrt(h2_ACE_ROI_1)*sqrt(h2_ACE_ROI_2))
  }
  
  return(df)
}

pfn_df <- function(df, save = FALSE) {
  #define some vectors for later
  unimodal <- c("PFN6", "PFN10", "PFN2", "PFN4", "PFN11", "PFN13", "PFN16", "VIS", "MOT", "AUD")
  heteromodal <- c("PFN12", "PFN8","PFN1", "PFN17", "PFN15", "PFN3", "PFN9", "PFN7", "PFN14", "PFN5", "PFN16", "PFN17", "FPN", "DAN", "VAN", "DMN")
  pfn_name_color <- data.table(PFN = c(paste("PFN", 1:17, sep = "")), Name = c("DMN","MOT", "FPN", "MOT", "DAN", "VIS", "VAN", "DMN", "VAN", "VIS", "MOT", "DMN", "MOT", "DAN", "FPN", "AUD", "FPN"), Color = c("#CB94A9", "#AEBACA", "#DFCA8B", "#7793BB", "#83AD73", "#785397", "#B97FE7", "#AF576F","#A241EA","#5A2C7B","#5E78AF", "#833852", "#445489", "#517549", "#BD986A", "#5647A1", "#AA7C45"))
  pfn_name_color$Modality <- ifelse(pfn_name_color$PFN %in% unimodal, "UNI", ifelse(pfn_name_color$PFN %in% heteromodal, "HETERO", NA))
  
  pfn7 <- data.table(PFN = c("VIS", "MOT", "DAN", "VAN", "FPN", "AUD", "DMN"), Name = c("VIS", "MOT", "DAN", "VAN", "FPN", "AUD", "DMN"), Color = rep(NA, 7), Modality = c("UNI", "UNI", "HETERO", "HETERO", "HETERO", "UNI", "HETERO"))
  pfn_name_color = rbind(pfn_name_color, pfn7)
  
  SA_Axis <- c(6,10,2,13,11,4,5,16,14,15,9,8,12,3,1,7,17)
  SA_Axis <- paste("PFN", SA_Axis, sep = "")
  
  #get name of input dataframe
  df_name <- deparse(substitute(df))
  if (class(df) %in% c("data.table", "data.frame")) {
    out_df = as.data.table(df)
  }
  
  ## for univariate
  if (sum(startsWith(names(out_df), "ROI")) == 1){
    #get SA order
    out_df$SA_order <- match(out_df$ROI_1, SA_Axis)
    
    #assign some important variables, if adding more non-numeric vars add here. Numeric ones add above.
    out_df <- cbind(out_df, pfn_name_color[match(out_df$ROI_1, pfn_name_color$PFN, NA),4])
    out_df <- cbind(out_df, pfn_name_color[match(out_df$ROI_1, pfn_name_color$PFN, NA),2])
    out_df <- cbind(out_df, pfn_name_color[match(out_df$ROI_1, pfn_name_color$PFN, NA),3])
    
    #get PFN no only as column to match with xifti file later
    out_df$pfn_no <- gsub("PFN","", out_df$ROI_1)
  } 
  #for bivariate
  else if(sum(startsWith(names(out_df), "ROI")) == 2){
    #calc SA distance variable btw bivar pairs
    out_df$SA_ROI1 <- as.numeric(match(out_df$ROI_1, SA_Axis, NA))
    out_df$SA_ROI2 <- as.numeric(match(out_df$ROI_2, SA_Axis, NA))
    out_df$SA_dist <- abs(out_df$SA_ROI1 - out_df$SA_ROI2)
    
    #naming column for ROI1 and ROI2
    out_df$ROI1_Name <- pfn_name_color[match(out_df$ROI_1, pfn_name_color$PFN, NA),"Name"]
    out_df$ROI2_Name <- pfn_name_color[match(out_df$ROI_2, pfn_name_color$PFN, NA),"Name"]
    out_df$ROI1_Color <- pfn_name_color[match(out_df$ROI_1, pfn_name_color$PFN, NA),"Color"]
    out_df$ROI2_Color <- pfn_name_color[match(out_df$ROI_2, pfn_name_color$PFN, NA),"Color"]
    out_df$ROI1_Mod <- pfn_name_color[match(out_df$ROI_1, pfn_name_color$PFN, NA),"Modality"]
    out_df$ROI2_Mod <- pfn_name_color[match(out_df$ROI_2, pfn_name_color$PFN, NA),"Modality"]
    
    #grouping variables for modality & network assignment
    out_df$class_same <- ifelse(out_df$ROI1_Name == out_df$ROI2_Name, TRUE, FALSE)
    out_df$mod_same <- ifelse(out_df$ROI1_Mod == out_df$ROI2_Mod, TRUE, FALSE)
    out_df$roi_both <- paste(out_df$ROI_1, out_df$ROI_2, sep = "_")
    out_df$roi_class_both <- paste(out_df$ROI1_Name, out_df$ROI2_Name, sep = "_")
  }
  else{stop(paste0("Currently fct only works for univar and bivar, might work on expanding it later. Or there are no ROIs in your data, which is another issue altogether."), call. = FALSE)}
  if(save == TRUE){
    filename = paste(df_name, format(Sys.time(), "%m_%d_%y"),"RData", sep = ".")
    assign(df_name, out_df)
    save(df_name, file = filename)
  }
  return(out_df)
}




##Functions - umx related

umxSummaryACE_EK <- function(model, digits = 2, comparison = NULL, std = TRUE, showRg = FALSE, CIs = FALSE, report = c("markdown", "html"), file = getOption("umx_auto_plot"), returnStd = FALSE, extended = FALSE, zero.print = ".", list_out = F, ...) {
  report = match.arg(report)
  commaSep = paste0(umx_set_separator(silent=TRUE), " ")
  
  out <- vector("list", 5)
  names(out) <- c("model", "Estimates", "NonStdEstimates","Correlations", "CIs")
  out$model <- model
  report   = match.arg(report)
  commaSep = paste0(umx_set_separator(silent=TRUE), " ")
  
  if(typeof(model) == "list"){ # call self recursively
    for(thisFit in model) {
      message("Output for Model: ", thisFit$name)
      umxSummaryACE(thisFit, digits = digits, file = file, showRg = showRg, std = std, comparison = comparison, CIs = CIs, returnStd = returnStd, extended = extended, zero.print = zero.print, report = report)
    }
  } else {
    umx_has_been_run(model, stop = TRUE)
    xmu_show_fit_or_comparison(model, comparison = comparison, digits = digits)
    selDVs = xmu_twin_get_var_names(model, trim= TRUE, twinOneOnly= TRUE)
    nVar   = length(selDVs)
    # TODO umxSummaryACE these already exist if a_std exists..
    # TODO replace all this with xmu_standardizeACE
    # Calculate standardized variance components
    a = mxEval(top.a, model) # Path coefficients
    c = mxEval(top.c, model)
    e = mxEval(top.e, model)
    A = mxEval(top.A, model) # Variances
    C = mxEval(top.C, model)
    E = mxEval(top.E, model)
    
    if(std){
      caption = paste0("Standardized parameter estimates from a ", dim(a)[2], "-factor Cholesky ACE model. ")
      Vtot = A + C + E;            # Total variance
      I    = diag(nVar);           # nVar Identity matrix
      SD   = solve(sqrt(I * Vtot)) # Inverse of diagonal matrix of standard deviations
      # (same as "(\sqrt(I.Vtot))~"
      
      # Standardized _path_ coefficients ready to be stacked together
      a_std  = SD %*% a; # Standardized path coefficients
      c_std  = SD %*% c;
      e_std  = SD %*% e;
      aClean = a_std
      cClean = c_std
      eClean = e_std
    } else {
      caption = paste0("Raw parameter estimates from a ", dim(a)[2], "-factor Cholesky ACE model. ")
      aClean = a
      cClean = c
      eClean = e
    }
    
    aClean[upper.tri(aClean)] = NA
    cClean[upper.tri(cClean)] = NA
    eClean[upper.tri(eClean)] = NA
    Estimates = data.frame(cbind(aClean, cClean, eClean), row.names = selDVs, stringsAsFactors = FALSE);
    
    if(model$top$dzCr$values == .25){
      colNames = c("a", "d", "e")
      caption = paste0(caption, "A: additive genetic; D: dominance effects; E: unique environment.")
    } else {
      colNames = c("a", "c", "e")
      caption = paste0(caption, "A: additive genetic; C: common environment; E: unique environment.")
    }
    names(Estimates) = paste0(rep(colNames, each = nVar), rep(1:nVar));
    out$Estimates <- Estimates
    umx_print(Estimates, digits = digits, caption = caption, report = report, zero.print = zero.print)
    xmu_twin_print_means(model = model, report = report)
    
    if(extended == TRUE) {
      aClean = a
      cClean = c
      eClean = e
      aClean[upper.tri(aClean)] = NA
      cClean[upper.tri(cClean)] = NA
      eClean[upper.tri(eClean)] = NA
      unStandardizedEstimates = data.frame(cbind(aClean, cClean, eClean), row.names = selDVs);
      names(unStandardizedEstimates) = paste0(rep(colNames, each = nVar), rep(1:nVar));
      out$NonStdEstimates <- unStandardizedEstimates
      umx_print(unStandardizedEstimates, caption = "Unstandardized Cholesky ACE model path coefficients", digits = digits, zero.print = zero.print)
    }
    
    if(showRg) {
      # Pre & post multiply covariance matrix by inverse of standard deviations
      NAmatrix = matrix(NA, nVar, nVar);
      rA = tryCatch(solve(sqrt(I*A)) %*% A %*% solve(sqrt(I*A)), error = function(err) return(NAmatrix)); # genetic correlations
      rC = tryCatch(solve(sqrt(I*C)) %*% C %*% solve(sqrt(I*C)), error = function(err) return(NAmatrix)); # C correlations
      rE = tryCatch(solve(sqrt(I*E)) %*% E %*% solve(sqrt(I*E)), error = function(err) return(NAmatrix)); # E correlations
      rAClean = rA
      rCClean = rC
      rEClean = rE
      rAClean[upper.tri(rAClean)] = NA
      rCClean[upper.tri(rCClean)] = NA
      rEClean[upper.tri(rEClean)] = NA
      genetic_correlations = data.frame(cbind(rAClean, rCClean, rEClean), row.names = selDVs);
      names(genetic_correlations) = selDVs
      # Make a nice table.
      names(genetic_correlations) = paste0(rep(c("rA", "rC", "rE"), each = nVar), rep(1:nVar));
      out$Correlations <- genetic_correlations
      umx_print(genetic_correlations, caption = "Genetic correlations", digits = digits, zero.print = zero.print)
    }
    hasCIs = umx_has_CIs(model)
    if(hasCIs && CIs) {
      # TODO umxACE CI code: Need to refactor into some function calls...
      # TODO and then add to umxSummaryIP and CP
      message("Creating CI-based report!")
      # CIs exist, get lower and upper CIs as a dataframe
      CIlist <- data.frame(model$output$confidenceIntervals)
      # Drop rows fixed to zero
      #CIlist = CIlist[(CIlist$lbound != 0 & CIlist$ubound != 0),]
      # Discard rows named NA
      #CIlist = CIlist[!grepl("^NA", row.names(CIlist)), ]
      # TODO fix for singleton CIs
      # THIS IS NOT NEEDED: confidenceIntervals come with estimate in the middle now...
      # These can be names ("top.a_std[1,1]") or labels ("a_r1c1")
      # imxEvalByName finds them both
      # outList = c();
      # for(aName in row.names(CIlist)) {
      # 	outList = append(outList, imxEvalByName(aName, model))
      # }
      # # Add estimates into the CIlist
      # CIlist$estimate = outList
      # reorder to match summary
      # CIlist = CIlist[, c("lbound", "estimate", "ubound")]
      CIlist$fullName = row.names(CIlist)
      # Initialise empty matrices for the CI results
      rows = dim(model$top$matrices$a$labels)[1]
      cols = dim(model$top$matrices$a$labels)[2]
      a_CI = c_CI = e_CI = matrix(NA, rows, cols)
      
      # Iterate over each CI
      labelList = imxGenerateLabels(model)	
      rowCount  = dim(CIlist)[1]
      # return(CIlist)
      for(n in 1:rowCount) { # n = 1
        thisName = row.names(CIlist)[n] # thisName = "a11"
        # Convert labels to [bracket] style
        if(!umx_has_square_brackets(thisName)) {
          nameParts = labelList[which(row.names(labelList) == thisName),]
          CIlist$fullName[n] = paste(nameParts$model, ".", nameParts$matrix, "[", nameParts$row, ",", nameParts$col, "]", sep = "")
        }
        fullName = CIlist$fullName[n]
        
        thisMatrixName = sub(".*\\.([^\\.]*)\\[.*", replacement = "\\1", x = fullName) # .matrix[
        thisMatrixRow  = as.numeric(sub(".*\\[(.*),(.*)\\]", replacement = "\\1", x = fullName))
        thisMatrixCol  = as.numeric(sub(".*\\[(.*),(.*)\\]", replacement = "\\2", x = fullName))
        #CIparts    = round(CIlist[n, c("estimate", "lbound", "ubound")], digits)
        CIparts    = CIlist
        thisString = paste0(CIparts[1], " [",CIparts[2], commaSep, CIparts[3], "]")
        
        if(grepl("^a", thisMatrixName)) {
          a_CI[thisMatrixRow, thisMatrixCol] = thisString
        } else if(grepl("^c", thisMatrixName)){
          c_CI[thisMatrixRow, thisMatrixCol] = thisString
        } else if(grepl("^e", thisMatrixName)){
          e_CI[thisMatrixRow, thisMatrixCol] = thisString
        } else{
          stop(paste("Illegal matrix name: must begin with a, c, or e. You sent: ", thisMatrixName))
        }
      }
      # TODO Check the merge of a_, c_ and e_CI INTO the output table works with more than one variable
      # TODO umxSummaryACE: Add option to use mxSE
      CI_Estimates = data.frame(cbind(a_CI, c_CI, e_CI), row.names = selDVs, stringsAsFactors = FALSE)
      names(CI_Estimates) = paste0(rep(colNames, each = nVar), rep(1:nVar));
      umx_print(CI_Estimates, digits = digits, zero.print = zero.print,  report=report, file = "tmpCI.html")
      xmu_twin_print_means(model, digits = digits, report = report)
      CI_Fit = model
      CI_Fit$top$a$values = a_CI
      CI_Fit$top$c$values = c_CI
      CI_Fit$top$e$values = e_CI
      out$CIs <- CIlist
    } # end Use CIs
  } # end list catcher?
  
  if(!is.na(file)) {
    # message("making dot file")
    if(hasCIs & CIs){
      umxPlotACE(CI_Fit, file = file, std = FALSE)
    } else {
      umxPlotACE(model, file = file, std = std)
    }
  }
  if(returnStd) {
    if(CIs){
      message("If you asked for CIs, returned model is not runnable (contains CIs not parameter values)")
    }
    xmu_standardize_ACE(model)
  }else if (list_out == F) {
    invisible(Estimates)
  }
  else if (list_out == T){
    return(out)
  }
}


umxReduceACE_EK <- function(model, report = c("markdown", "inline", "html", "report"), intervals = TRUE, baseFileName = "tmp", tryHard = c("yes", "no", "ordinal", "search"), silent=FALSE, digits = 2, list_out = F, ...) {
  out <- vector("list", 2)
  names(out) <- c("bestModel", "modelCompare")
  report  = match.arg(report)
  tryHard = match.arg(tryHard)
  if(silent){
    oldSilent = umx_set_silent(TRUE)
  }else{
    oldSilent = FALSE
  }
  oldAutoPlot = umx_set_auto_plot(FALSE, silent = TRUE)
  if(model$top$dzCr$values == 1){
    message("You gave me an ACE model")		
    ACE = model
    ADE = umxModify(model, 'dzCr_r1c1', value = .25, name = "ADE", tryHard = tryHard)
    if(-2*logLik(ACE) > -2*logLik(ADE)){
      CE = umxModify(ADE, regex = "a_r[0-9]+c[0-9]+" , name = "DE", tryHard = tryHard)
      AE = umxModify(ADE, regex = "c_r[0-9]+c[0-9]+" , name = "AE", tryHard = tryHard)
      E = umxModify( AE, regex = "a_r[0-9]+c[0-9]+" , name =  "E", tryHard = tryHard)
      message("A dominance model is preferred, set dzCr = 0.25")
    }else{
      CE = umxModify(ACE, regex = "a_r[0-9]+c[0-9]+" , name = "CE", tryHard = tryHard)
      AE = umxModify(ACE, regex = "c_r[0-9]+c[0-9]+" , name = "AE", tryHard = tryHard)
      E = umxModify( AE, regex = "a_r[0-9]+c[0-9]+" , name =  "E", tryHard = tryHard)
    }
  }else if(model$top$dzCr$values == .25){
    if(model$name=="ACE"){
      message("You gave me an ADE model, but it was called 'ACE'. I have renamed it ADE for the purposes of clarity in model reduction.")
      model = mxRename(model, newname = "ADE", oldname = "ACE")
    } else {
      message("You gave me an ADE model.")
    }
    ADE = model
    ACE = umxModify(ADE, 'dzCr_r1c1', value = 1, name = "ACE", tryHard = tryHard)
    AE  = umxModify(ADE, regex = "c_r[0-9]+c[0-9]+" , name = "AE", tryHard = tryHard)
    E  = umxModify( AE, regex = "a_r[0-9]+c[0-9]+" , name =  "E", tryHard = tryHard)
    if(-2*logLik(ADE) > -2*logLik(ACE)){
      CE = umxModify(ACE, regex = "a_r[0-9]+c[0-9]+" , name = "CE", tryHard=tryHard)
      message("An ACE model is preferred, set dzCr = 1.0")
    }else{
      CE = umxModify(ADE, regex = "a_r[0-9]+c[0-9]+" , name = "DE", tryHard=tryHard)
    }
  }else{
    stop(model$top$dzCr$values, " is an odd number for dzCr, isn't it? I was expecting 1 (C) or .25 (D)",
         "\nPerhaps you're a John Loehlin (requiescat in pace :-( ) fan, and are doing an assortative mating test? e-mail me to get this added here.")
    # TODO umxReduceACE handle odd values of dzCr as assortative mating etc.?
    bestModel = model
    out$bestModel <- bestModel
  }
  # = Show fit table =
  tmp = data.frame(matrix(nrow=5,ncol=4))
  names(tmp) = c("a"              , "c"                 , "e"                 , "d")
  tmp[1,] = c(ACE$top$a_std$result, ACE$top$c_std$result, ACE$top$e_std$result, NA)
  tmp[2,] = c(ADE$top$a_std$result, NA                  , ADE$top$e_std$result, ADE$top$c_std$result)
  tmp[3,] = c( NA                 ,  CE$top$c_std$result,  CE$top$e_std$result, NA)
  tmp[4,] = c( AE$top$a_std$result, NA                  ,  AE$top$e_std$result, NA)
  tmp[5,] = c( NA                 , NA                  ,   E$top$e_std$result, NA)
  
  biggles = umxCompare(ACE, c(ADE, CE, AE, E), all = TRUE, report = report, silent=TRUE)
  tmp2 = cbind(biggles[, 1, drop = FALSE], tmp, biggles[, 2:dim(biggles)[2] ] )
  out$modelCompare <- tmp2
  umx_print(tmp2, digits = digits, report = report)
  whichBest = which.min(AIC(ACE, ADE, CE, AE)[,"AIC"])[1]
  bestModel = list(ACE, ADE, CE, AE)[[whichBest]]
  out$bestModel <- bestModel
  message("Among ACE, ADE, CE, and AE models ", omxQuotes(bestModel$name), " fit best according to AIC.")
  # Probabilities according to AIC MuMIn::Weights (Wagenmakers et al https://pubmed.ncbi.nlm.nih.gov/15117008/ )
  aic.weights = round(Weights(AIC(ACE, ADE, CE, AE)[,"AIC"]), 2)
  aic.names   = namez(c(ACE, ADE, CE, AE))
  message("Conditional AIC probability {Wagenmakers, 2004, 192-196}  indicates relative model support as", 
          omxQuotes(aic.names), " respectively are: ", 
          omxQuotes(aic.weights), " Using MuMIn::Weights(AIC()).")
  message(paste0(aic.names," (", aic.weights, "%)"))
  if(intervals){
    bestModel = mxRun(bestModel, intervals = intervals)
  }
  out$bestModel <- bestModel
  umx_set_auto_plot(oldAutoPlot, silent = TRUE)
  umx_set_silent(oldSilent)
  if (list_out == F) {
    invisible(bestModel)
  } else if(list_out == T) {
    return(out)
  }
}

umxCompare_EK <- function(base = NULL, comparison = NULL, all = TRUE, digits = 3, report = c("markdown", "html", "inline"), compareWeightedAIC = FALSE, silent = FALSE, file = "tmp.html", list_out = F) {
  out <- vector("list", 3)
  names(out) <- c("bestModel", "modelCompare", "tablePub")
  report = match.arg(report)
  if(umx_is_MxModel(all)){
    stop("Provide all comparison models as a c() (You provided a model as input to 'all', and I'm guessing that's a mistake)")
  }
  if(is.null(comparison)){
    comparison = base
  } else if (is.null(base)) {
    stop("You must provide at least a base model for umxCompare")
  }
  if(length(base) == 1) {
    if(typeof(base) == "list"){
      base = base[[1]]
    }
    if(!umx_has_been_run(base)){
      warning("Base model not run yet!")		
    }
  }
  if(length(comparison) == 1) {
    if(typeof(comparison) == "list"){
      comparison = comparison[[1]]
    }
    if(!umx_has_been_run(comparison)){
      warning("Comparison model has not been run!")		
    }
  }
  tableOut = mxCompare(base = base, comparison = comparison, all = all)
  tableOut = as.data.frame(tableOut)
  out$modelCompare <- tableOut
  # | base    | comparison    | ep | minus2LL | df  | AIC      | diffLL   | diffdf |p     |fit       |fitUnits |diffFit |chisq     |SBchisq |
  # |:--------|:--------------|---:|:---------|:----|---------:|:---------|:-------|:-----|:---------|:--------|:-------|:---------|:-------|
  # |DWLS     |               |  6 |          |0    | 12.00000 |          |        |      |0         |r'Wr     |        |0         |        |
  # |DWLS     |drop_l2mpg     |  5 |          |1    | 14.49542 |          |1       |      |4.4954186 |r'Wr     |        |4.4954186 |        |
  
  # | base    | comparison    | ep | minus2LL | df  | AIC      | diffLL   | diffdf | p    |
  #    1            2           3     4          5     6          7          8        9     
  # | twinSat | <NA>          | 13 | 333.0781 | 149 | 35.07809 | NA       | NA     | NA   |
  # | twinSat | betaSetToZero | 10 | 351.6486 | 152 | 47.64858 | 18.57049 | 3      | 0.01 |
  
  # Pre Feb 2021 version 2.18.1.233
  if(packageVersion("OpenMx") < "2.18.1.233"){
    # old format mxCompare
    tablePub = tableOut[, c("comparison", "ep", "diffLL" , "diffdf", "p", "AIC", "base")]
    names(tablePub)     = c("comparison", "ep", "diffFit", "diffdf", "p", "AIC", "base")
    tablePub$fitUnits = ""
  } else {
    # new format mxCompare
    tablePub = tableOut[, c("comparison", "ep", "diffFit", "diffdf", "p", "AIC", "base", "fitUnits")]
    # str(tmp@results)
    # 'data.frame':	2 obs. of  13 variables:
    #  $ base      : chr  "tim" "tim"
    #  $ comparison: chr  NA "tim"
    #  $ ep        : num  9 9
    #  $ df        : num  87 87
    #  $ diffdf    : num  NA 0
    #  $ fit       : num  330 330
    #  $ fitUnits  : chr  "-2lnL" "-2lnL"
    #  $ diffFit   : num  NA 0
    #  $ AIC       : num  348 348
    #  $ p         : num  NA NA
    #  $ minus2LL  : num  330 330
    #  $ diffLL    : num  NA 0
    #  $ SBchisq   : num  NA NA
  }
  
  # Subtract row-1 AIC from all values and place the resulting deltaAIC column after AIC 
  tablePub$deltaAIC = tablePub[, "AIC"] - tablePub[1, "AIC"]
  tablePub = tablePub[,c("comparison", "ep", "diffFit", "diffdf", "p", "AIC", "deltaAIC", "base", "fitUnits")]
  
  # c("1: Comparison", "2: Base", "3: EP", "4: AIC", "5: &Delta; -2LL", "6: &Delta; df", "7: p")
  # U+2206 = math delta
  # Fix problem where base model has compare set to its own name, and name set to NA
  nRows = dim(tablePub)[1]
  for (i in 1:nRows) {
    if(is.na(tablePub[i, "comparison"])){
      tablePub[i, "comparison"] = tablePub[i, "base"]
      tablePub[i, "base"] = NA
    }
  }
  tablePub[, "p"] = umx_APA_pval(tablePub[, "p"], min = (1/ 10^3), digits = digits, addComparison = NA)
  if(report == "inline"){
    n_rows = dim(tablePub)[1]
    for (i in 1:n_rows) {
      thisPValue = tableOut[i, "p"]
      if(!is.na(thisPValue) && !is.nan(thisPValue)){
        if(tableOut[i, "p"] < .05){
          this = ". This caused a significant loss of fit "
        } else {
          this = ". This did not lower fit significantly "
        }
        inlineMsg = paste0("The hypothesis that ", omxQuotes(tablePub[i, "comparison"]), 
                           " was tested by dropping ", tablePub[i, "comparison"],
                           " from ", omxQuotes(tablePub[i, "base"]), 
                           this, "(\u03C7\u00B2(", tablePub[i, "diffdf"], ") = ", round(tablePub[i, "diffFit"], 2), # \u03A7 = Chi \u00B2 = superscript 2
                           ", p = ", tablePub[i, "p"], ": AIC = ", round(tablePub[i, "AIC"], digits), " change in AIC = ", round(tablePub[i, "deltaAIC"], digits), ")."
        )
        if(!silent){
          cat(inlineMsg)
        }
      }
    }
  }
  
  # Rename for printing to console
  names(tablePub) = c("Model", "EP", "\u0394 Fit" , "\u0394 df" , "p", "AIC", "\u0394 AIC", "Compare with Model", "Fit units")
  
  if(report == "inline"){ report= "markdown"}
  out$tablePub <- tablePub
  if(!silent){
    umx_print(tablePub, digits = digits, zero.print = "0", caption = "Table of Model Comparisons", report = report)
  }
  # htmlNames       = c("Model", "EP", "&Delta; -2LL", "&Delta; df", "p", "AIC", "&Delta AIC", "Compare with Model")
  # if(report == "html"){
  # 	tableHTML = tablePub
  # 	names(tableHTML) = htmlNames
  # 	print(xtable::xtable(tableHTML), type = "HTML", file = file, sanitize.text.function = function(x){x})
  # 	umx_open(file)
  # }
  if(compareWeightedAIC){
    modelList = c(base, comparison)
    # get list of AICs
    AIClist = c()
    for (i in modelList) {
      AIClist = c(AIClist, AIC(i))
    }
    whichBest = which.min(AIClist)
    bestModel = modelList[[whichBest]]
    out$bestModel <- bestModel
    # Probabilities according to AIC MuMIn::Weights (Wagenmakers et al https://pubmed.ncbi.nlm.nih.gov/15117008/ )
    aic.weights = round(Weights(AIClist), 2)
    if(!silent){
      cat("The ", omxQuotes(bestModel$name), " model is the best fitting model according to AIC.")
      cat("AIC weight-based  {Wagenmakers, 2004, 192-196} conditional probabilities of being the best model for ", 
          omxQuotes(namez(modelList)), " respectively are: ", 
          omxQuotes(aic.weights), " Using MuMIn::Weights(AIC()).")	
    }
    
  }
  if(list_out == F){
    invisible(tablePub)
  } else if(list_out == T){
    return(out)
  }
}

umxSexLim_EK <- function (name = "sexlim", selDVs, mzmData, dzmData, mzfData, dzfData, dzoData = NULL , sep = NA, A_or_C = c("A", "C"), sexlim = c("Nonscalar", "Scalar", "Homogeneity"), dzAr = 0.5, dzCr = 1, autoRun = getOption("umx_auto_run"),  tryHard = c("no", "yes", "ordinal", "search"), optimizer = NULL) {
  message("umxSexLim is a beta feature. Some things are broken. If any desired stats are not presented, let me know what's missing")
  A_or_C = match.arg(A_or_C)
  sexlim = match.arg(sexlim)
  tryHard = match.arg(tryHard)
  nSib = 2
  xmu_twin_check(selDVs = selDVs, sep = sep, dzData = dzmData, 
                 mzData = mzmData, enforceSep = TRUE, nSib = nSib, optimizer = optimizer)
  if (name == "sexlim") {
    if (dzCr == 0.25) {
      name = paste0(sexlim, "ADE")
    }
    else {
      name = paste0(sexlim)
    }
  }
  nVar = length(selDVs)
  selVars = umx_paste_names(selDVs, sep = sep, suffixes = 1:2)
  umx_check_names(selVars, data = mzmData, die = TRUE)
  mzmData = mzmData[, selVars]
  umx_check_names(selVars, data = dzmData, die = TRUE)
  dzmData = dzmData[, selVars]
  umx_check_names(selVars, data = mzfData, die = TRUE)
  mzfData = mzfData[, selVars]
  umx_check_names(selVars, data = dzfData, die = TRUE)
  dzfData = dzfData[, selVars]
  if (!is.null(dzoData)){
    umx_check_names(selVars, data = dzoData, die = TRUE)
    dzoData = dzoData[, selVars]
  }
  obsMean = umx_means(mzmData[, selVars[1:nVar], drop = FALSE])
  varStarts = umx_var(mzmData[, selVars[1:nVar], drop = FALSE], 
                      format = "diag", ordVar = 1, use = "pairwise.complete.obs")
  if (nVar == 1) {
    varStarts = sqrt(varStarts)/3
  }
  else {
    varStarts = t(chol(diag(varStarts/3)))
  }
  varStarts = matrix(varStarts, nVar, nVar)
  if (A_or_C == "A") {
    Rao = umxMatrix("Rao", "Full", nrow = nVar, ncol = nVar, 
                    free = TRUE, values = 1, lbound = -1, ubound = 1)
    Rco = umxMatrix("Rco", "Stand", nrow = nVar, ncol = nVar, 
                    free = TRUE, values = 0.4, lbound = -1, ubound = 1)
  }
  else if (A_or_C == "C") {
    Rao = umxMatrix("Rao", "Stand", nrow = nVar, ncol = nVar, 
                    free = TRUE, values = 0.4, lbound = -1, ubound = 1)
    Rco = umxMatrix("Rco", "Full", nrow = nVar, ncol = nVar, 
                    free = TRUE, values = 1, lbound = -1, ubound = 1)
  }
  if(!is.null(dzoData)) {
    model = mxModel(name, mxModel("top", umxMatrix("dzAr", "Full", 
                                                   1, 1, free = FALSE, values = dzAr), umxMatrix("dzCr", 
                                                                                                 "Full", 1, 1, free = FALSE, values = dzCr), umxMatrix("am", 
                                                                                                                                                       "Diag", nrow = nVar, free = TRUE, values = varStarts, 
                                                                                                                                                       lbound = 1e-04), umxMatrix("cm", "Diag", nrow = nVar, 
                                                                                                                                                                                  free = TRUE, values = varStarts, lbound = 1e-04), umxMatrix("em", 
                                                                                                                                                                                                                                              "Diag", nrow = nVar, free = TRUE, values = varStarts, 
                                                                                                                                                                                                                                              lbound = 1e-04), umxMatrix("af", "Diag", nrow = nVar, 
                                                                                                                                                                                                                                                                         free = TRUE, values = varStarts, lbound = 1e-04), umxMatrix("cf", 
                                                                                                                                                                                                                                                                                                                                     "Diag", nrow = nVar, free = TRUE, values = varStarts, 
                                                                                                                                                                                                                                                                                                                                     lbound = 1e-04), umxMatrix("ef", "Diag", nrow = nVar, 
                                                                                                                                                                                                                                                                                                                                                                free = TRUE, values = varStarts, lbound = 1e-04), umxMatrix("Ram", 
                                                                                                                                                                                                                                                                                                                                                                                                                            "Stand", nrow = nVar, free = TRUE, values = 0.4, lbound = -1, 
                                                                                                                                                                                                                                                                                                                                                                                                                            ubound = 1), umxMatrix("Rcm", "Stand", nrow = nVar, free = TRUE, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                   values = 0.4, lbound = -1, ubound = 1), umxMatrix("Rem", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     "Stand", nrow = nVar, free = TRUE, values = 0.4, lbound = -1, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     ubound = 1), umxMatrix("Raf", "Stand", nrow = nVar, free = TRUE, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            values = 0.4, lbound = -1, ubound = 1), umxMatrix("Rcf", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              "Stand", nrow = nVar, free = TRUE, values = 0.4, lbound = -1, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              ubound = 1), umxMatrix("Ref", "Stand", nrow = nVar, free = TRUE, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     values = 0.4, lbound = -1, ubound = 1), Rao, Rco, mxAlgebra(name = "Am", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 Ram %&% am), mxAlgebra(name = "Cm", Rcm %&% cm), mxAlgebra(name = "Em", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            Rem %&% em), mxAlgebra(name = "Af", Raf %&% af), mxAlgebra(name = "Cf", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       Rcf %&% cf), mxAlgebra(name = "Ef", Ref %&% ef), mxAlgebra(name = "Amf", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  am %*% (Rao) %*% t(af)), mxAlgebra(name = "Cmf", cm %*% 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       (Rco) %*% t(cf)), umxMatrix("pos1by6", "Full", nrow = 1, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   ncol = 6, free = FALSE, values = 1e-04), mxAlgebra(name = "minCor", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      cbind(min(eigenval(Ram)), min(eigenval(Rcm)), min(eigenval(Rem)), 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            min(eigenval(Raf)), min(eigenval(Rcf)), min(eigenval(Ref)))), 
                                  mxConstraint(name = "Keep_it_Positive", minCor > pos1by6), 
                                  umxMatrix("I", "Iden", nrow = nVar), mxAlgebra(name = "Vm", 
                                                                                 Am + Cm + Em, list(selDVs, selDVs)), mxAlgebra(name = "Vf", 
                                                                                                                                Af + Cf + Ef, list(selDVs, selDVs)), mxAlgebra(name = "iSDm", 
                                                                                                                                                                               solve(sqrt(I * Vm))), mxAlgebra(name = "iSDf", solve(sqrt(I * 
                                                                                                                                                                                                                                           Vf))), mxAlgebra(name = "AmStd", Am/Vm, dimnames = list(selDVs, 
                                                                                                                                                                                                                                                                                                   paste0(selDVs, "AmStd"))), mxAlgebra(name = "CmStd", 
                                                                                                                                                                                                                                                                                                                                        Cm/Vm, dimnames = list(selDVs, paste0(selDVs, "CmStd"))), 
                                  mxAlgebra(name = "EmStd", Em/Vm, dimnames = list(selDVs, 
                                                                                   paste0(selDVs, "EmStd"))), mxAlgebra(name = "AfStd", 
                                                                                                                        Af/Vf, dimnames = list(selDVs, paste0(selDVs, "AfStd"))), 
                                  mxAlgebra(name = "CfStd", Cf/Vf, dimnames = list(selDVs, 
                                                                                   paste0(selDVs, "CfStd"))), mxAlgebra(name = "EfStd", 
                                                                                                                        Ef/Vf, dimnames = list(selDVs, paste0(selDVs, "EfStd"))), 
                                  umxMatrix("expMeanGm", "Full", nrow = 1, ncol = nVar * 
                                              2, free = TRUE, values = obsMean, labels = paste0(selDVs, 
                                                                                                "_mean_m")), umxMatrix("expMeanGf", "Full", nrow = 1, 
                                                                                                                       ncol = nVar * 2, free = TRUE, values = obsMean, labels = paste0(selDVs, 
                                                                                                                                                                                       "_mean_f")), umxMatrix("expMeanGo", "Full", nrow = 1, 
                                                                                                                                                                                                              ncol = nVar * 2, free = TRUE, values = obsMean, labels = paste0(selDVs, 
                                                                                                                                                                                                                                                                              rep(c("_mean_m", "_mean_f"), each = nVar))), 
                                  mxAlgebra(name = "expCovMZm", rbind(cbind(Vm, Am + Cm), 
                                                                      cbind(Am + Cm, Vm))), mxAlgebra(name = "expCovDZm", 
                                                                                                      rbind(cbind(Vm, dzAr %x% Am + dzCr %x% Cm), cbind(dzAr %x% 
                                                                                                                                                          Am + dzCr %x% Cm, Vm))), mxAlgebra(name = "expCovMZf", 
                                                                                                                                                                                             rbind(cbind(Vf, Af + Cf), cbind(Af + Cf, Vf))), mxAlgebra(name = "expCovDZf", 
                                                                                                                                                                                                                                                       rbind(cbind(Vf, dzAr %x% Af + dzCr %x% Cf), cbind(dzAr %x% 
                                                                                                                                                                                                                                                                                                           Af + dzCr %x% Cf, Vf))), mxAlgebra(name = "expCovDZo", 
                                                                                                                                                                                                                                                                                                                                              rbind(cbind(Vm, dzAr %x% Amf + dzCr %x% Cmf), cbind(dzAr %x% 
                                                                                                                                                                                                                                                                                                                                                                                                    t(Amf) + t(Cmf), Vf)))), mxModel("MZm", mxExpectationNormal("top.expCovMZm", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                means = "top.expMeanGm", dimnames = selVars), mxFitFunctionML(), 
                                                                                                                                                                                                                                                                                                                                                                                                                                     mxData(mzmData, type = "raw")), mxModel("DZm", mxExpectationNormal("top.expCovDZm", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        means = "top.expMeanGm", dimnames = selVars), mxFitFunctionML(), 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                             mxData(dzmData, type = "raw")), mxModel("MZf", mxExpectationNormal("top.expCovMZf", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                means = "top.expMeanGf", dimnames = selVars), mxFitFunctionML(), 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     mxData(mzfData, type = "raw")), mxModel("DZf", mxExpectationNormal("top.expCovDZf", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        means = "top.expMeanGf", dimnames = selVars), mxFitFunctionML(), 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             mxData(dzfData, type = "raw")), mxModel("DZo", mxExpectationNormal("top.expCovDZo", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                means = "top.expMeanGo", dimnames = selVars), mxFitFunctionML(), 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     mxData(dzoData, type = "raw")), mxFitFunctionMultigroup(c("MZf", "DZf", "MZm", "DZm", "DZo")))
  }
  else {
    model = mxModel(name, mxModel("top", umxMatrix("dzAr", "Full", 
                                                   1, 1, free = FALSE, values = dzAr), umxMatrix("dzCr", 
                                                                                                 "Full", 1, 1, free = FALSE, values = dzCr), umxMatrix("am", 
                                                                                                                                                       "Diag", nrow = nVar, free = TRUE, values = varStarts, 
                                                                                                                                                       lbound = 1e-04), umxMatrix("cm", "Diag", nrow = nVar, 
                                                                                                                                                                                  free = TRUE, values = varStarts, lbound = 1e-04), umxMatrix("em", 
                                                                                                                                                                                                                                              "Diag", nrow = nVar, free = TRUE, values = varStarts, 
                                                                                                                                                                                                                                              lbound = 1e-04), umxMatrix("af", "Diag", nrow = nVar, 
                                                                                                                                                                                                                                                                         free = TRUE, values = varStarts, lbound = 1e-04), umxMatrix("cf", 
                                                                                                                                                                                                                                                                                                                                     "Diag", nrow = nVar, free = TRUE, values = varStarts, 
                                                                                                                                                                                                                                                                                                                                     lbound = 1e-04), umxMatrix("ef", "Diag", nrow = nVar, 
                                                                                                                                                                                                                                                                                                                                                                free = TRUE, values = varStarts, lbound = 1e-04), umxMatrix("Ram", 
                                                                                                                                                                                                                                                                                                                                                                                                                            "Stand", nrow = nVar, free = TRUE, values = 0.4, lbound = -1, 
                                                                                                                                                                                                                                                                                                                                                                                                                            ubound = 1), umxMatrix("Rcm", "Stand", nrow = nVar, free = TRUE, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                   values = 0.4, lbound = -1, ubound = 1), umxMatrix("Rem", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     "Stand", nrow = nVar, free = TRUE, values = 0.4, lbound = -1, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     ubound = 1), umxMatrix("Raf", "Stand", nrow = nVar, free = TRUE, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            values = 0.4, lbound = -1, ubound = 1), umxMatrix("Rcf", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              "Stand", nrow = nVar, free = TRUE, values = 0.4, lbound = -1, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              ubound = 1), umxMatrix("Ref", "Stand", nrow = nVar, free = TRUE, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     values = 0.4, lbound = -1, ubound = 1), Rao, Rco, mxAlgebra(name = "Am", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 Ram %&% am), mxAlgebra(name = "Cm", Rcm %&% cm), mxAlgebra(name = "Em", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            Rem %&% em), mxAlgebra(name = "Af", Raf %&% af), mxAlgebra(name = "Cf", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       Rcf %&% cf), mxAlgebra(name = "Ef", Ref %&% ef), mxAlgebra(name = "Amf", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  am %*% (Rao) %*% t(af)), mxAlgebra(name = "Cmf", cm %*% 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       (Rco) %*% t(cf)), umxMatrix("pos1by6", "Full", nrow = 1, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   ncol = 6, free = FALSE, values = 1e-04), mxAlgebra(name = "minCor", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      cbind(min(eigenval(Ram)), min(eigenval(Rcm)), min(eigenval(Rem)), 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            min(eigenval(Raf)), min(eigenval(Rcf)), min(eigenval(Ref)))), 
                                  mxConstraint(name = "Keep_it_Positive", minCor > pos1by6), 
                                  umxMatrix("I", "Iden", nrow = nVar), mxAlgebra(name = "Vm", 
                                                                                 Am + Cm + Em, list(selDVs, selDVs)), mxAlgebra(name = "Vf", 
                                                                                                                                Af + Cf + Ef, list(selDVs, selDVs)), mxAlgebra(name = "iSDm", 
                                                                                                                                                                               solve(sqrt(I * Vm))), mxAlgebra(name = "iSDf", solve(sqrt(I * 
                                                                                                                                                                                                                                           Vf))), mxAlgebra(name = "AmStd", Am/Vm, dimnames = list(selDVs, 
                                                                                                                                                                                                                                                                                                   paste0(selDVs, "AmStd"))), mxAlgebra(name = "CmStd", 
                                                                                                                                                                                                                                                                                                                                        Cm/Vm, dimnames = list(selDVs, paste0(selDVs, "CmStd"))), 
                                  mxAlgebra(name = "EmStd", Em/Vm, dimnames = list(selDVs, 
                                                                                   paste0(selDVs, "EmStd"))), mxAlgebra(name = "AfStd", 
                                                                                                                        Af/Vf, dimnames = list(selDVs, paste0(selDVs, "AfStd"))), 
                                  mxAlgebra(name = "CfStd", Cf/Vf, dimnames = list(selDVs, 
                                                                                   paste0(selDVs, "CfStd"))), mxAlgebra(name = "EfStd", 
                                                                                                                        Ef/Vf, dimnames = list(selDVs, paste0(selDVs, "EfStd"))), 
                                  umxMatrix("expMeanGm", "Full", nrow = 1, ncol = nVar * 
                                              2, free = TRUE, values = obsMean, labels = paste0(selDVs, 
                                                                                                "_mean_m")), umxMatrix("expMeanGf", "Full", nrow = 1, 
                                                                                                                       ncol = nVar * 2, free = TRUE, values = obsMean, labels = paste0(selDVs, 
                                                                                                                                                                                       "_mean_f")), umxMatrix("expMeanGo", "Full", nrow = 1, 
                                                                                                                                                                                                              ncol = nVar * 2, free = TRUE, values = obsMean, labels = paste0(selDVs, 
                                                                                                                                                                                                                                                                              rep(c("_mean_m", "_mean_f"), each = nVar))), 
                                  mxAlgebra(name = "expCovMZm", rbind(cbind(Vm, Am + Cm), 
                                                                      cbind(Am + Cm, Vm))), mxAlgebra(name = "expCovDZm", 
                                                                                                      rbind(cbind(Vm, dzAr %x% Am + dzCr %x% Cm), cbind(dzAr %x% 
                                                                                                                                                          Am + dzCr %x% Cm, Vm))), mxAlgebra(name = "expCovMZf", 
                                                                                                                                                                                             rbind(cbind(Vf, Af + Cf), cbind(Af + Cf, Vf))), mxAlgebra(name = "expCovDZf", 
                                                                                                                                                                                                                                                       rbind(cbind(Vf, dzAr %x% Af + dzCr %x% Cf), cbind(dzAr %x% 
                                                                                                                                                                                                                                                                                                           Af + dzCr %x% Cf, Vf)))), mxModel("MZm", mxExpectationNormal("top.expCovMZm", 
                                                                                                                                                                                                                                                                                                                                                                        means = "top.expMeanGm", dimnames = selVars), mxFitFunctionML(), 
                                                                                                                                                                                                                                                                                                                                             mxData(mzmData, type = "raw")), mxModel("DZm", mxExpectationNormal("top.expCovDZm", 
                                                                                                                                                                                                                                                                                                                                                                                                                means = "top.expMeanGm", dimnames = selVars), mxFitFunctionML(), 
                                                                                                                                                                                                                                                                                                                                                                                     mxData(dzmData, type = "raw")), mxModel("MZf", mxExpectationNormal("top.expCovMZf", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                        means = "top.expMeanGf", dimnames = selVars), mxFitFunctionML(), 
                                                                                                                                                                                                                                                                                                                                                                                                                             mxData(mzfData, type = "raw")), mxModel("DZf", mxExpectationNormal("top.expCovDZf", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                means = "top.expMeanGf", dimnames = selVars), mxFitFunctionML(), 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                     mxData(dzfData, type = "raw")), mxFitFunctionMultigroup(c("MZf","DZf", "MZm", "DZm")))
  }
  if (sexlim == "Nonscalar") {
    if (A_or_C == "A") {
      if ("^Rc[fmo](_.*)$" %in% umxGetParameters(model)) {
        model = umxModify(model, regex = "^Rc[fmo](_.*)$", 
                          newlabels = "Rc\\1", autoRun = FALSE)
      }
    }
    else if (A_or_C == "C") {
      if ("^Ra[fmo](_.*)$" %in% umxGetParameters(model)) {
        model = umxModify(model, regex = "^Ra[fmo](_.*)$", 
                          newlabels = "Ra\\1", autoRun = FALSE)
      }
    }
  }
  else if (sexlim %in% c("Scalar", "Homogeneity")) {
    model = umxModify(model, regex = "^(R[ace])[f|m|o]_", 
                      newlabels = "\\1_", autoRun = FALSE)
    model = mxModel(model, mxModel(model$top, umxMatrix("Rao", 
                                                        "Stand", nrow = nVar, free = TRUE, values = 0.2, 
                                                        baseName = "Ra", lbound = -1, ubound = 1), umxMatrix("Rco", 
                                                                                                             "Stand", nrow = nVar, free = TRUE, values = 0.2, 
                                                                                                             baseName = "Rc", lbound = -1, ubound = 1)))
    if (sexlim == "Homogeneity") {
      model = umxModify(model, regex = "^a[fm]_", newlabels = "a_", 
                        autoRun = FALSE)
      model = umxModify(model, regex = "^c[fm]_", newlabels = "c_", 
                        autoRun = FALSE)
      model = umxModify(model, regex = "^e[fm]_", newlabels = "e_", 
                        autoRun = FALSE)
    }
  }
  model = omxAssignFirstParameters(model)
  model = as(model, "MxModelSexLim")
  model = xmu_safe_run_summary(model, autoRun = autoRun, tryHard = tryHard)
  invisible(model)
}


##Functions- misc

file_to_dataframe <- function(pattern, header, directory){
  setwd(directory)
  all_files <- list.files(pattern = pattern)
  all_data <- lapply(all_files, read.table,fill = TRUE, header = F)
  dataframe <- do.call(rbind, all_data)
  names(dataframe) <- header
  return(dataframes)
}

get_dist_matrix <- function(id_df, dist_df, id_var, dist_var, na_var){
  id_1 <- paste(id_var, 1, sep = "_")
  id_2 <- paste(id_var, 2, sep = "_")
  n <- dim(id_df)[1]
  m<- data.frame(matrix(ncol = n, nrow = n))
  colnames(m) <- id_df[[id_1]]
  rownames(m) <- id_df[[id_1]]
  dist_df$pairs_str <- paste(dist_df[[id_1]], dist_df[[id_2]], sep = "_")
  distvar_no <- which(colnames(dist_df) == dist_var)
  navar_no <- which(colnames(dist_df) == na_var)
  for(i in 1:length(m)){
    for(j in 1:length(m)){
      if(i != j){
        val <- as.numeric(dist_df[which((dist_df[[id_1]] == rownames(m)[i]) & (dist_df[[id_2]] == colnames(m)[j])), distvar_no])
        #val <- as.numeric( dist_df[[dist_var]][str_detect(dist_df$pairs_str, paste(rownames(m)[i])) & str_detect(dist_df$pairs_str, paste(colnames(m)[j]))]))
        m[i,j] = val
        if (is.na(val)){  
          na_val <- as.numeric(dist_df[which((dist_df[[id_1]] == rownames(m)[i]) & (dist_df[[id_2]] == colnames(m)[j])), navar_no])
          #val <- as.numeric( dist_df[[dist_var]][str_detect(dist_df$pairs_str, paste(rownames(m)[i])) & str_detect(dist_df$pairs_str, paste(colnames(m)[j]))]))
          m[i,j] = na_val
        }
      }else{m[i,j] = 0}
    }
  } 
  return(m)
}

#Doesn't actually show the plot when view == TRUE. :((((((
make_dscalar <- function(xii, df, match_val, map_val, save = FALSE, view = FALSE, colors = "Reds", zlim = NULL, color_mode = "sequential"){
  df_name <- deparse(substitute(df))
  if ("data.table" %in% class(df)) {
    df = as.data.frame(df)
  }
  if(!is.xifti(xii)){
    xii <- read_cifti(
      xii, brainstructures="all", 
      surfL_fname= NULL , surfR_fname= NULL,
      resamp_res= NULL)
  }
  if(!is.xifti(xii)) {
    stop(paste0("Please provided loaded xifti file or valid path to xifti file."), call. = FALSE)
  }
  if(!is.numeric(df[,map_val])){
    df[,map_val] <- as.factor(df[,map_val])
    factor_levels <- c(levels(df[,map_val]))
    df[,map_val] <- as.numeric(df[,map_val])
  }
  # Access CIFTI data ---------------------------------------------
  cortexL_mwall <- xii$meta$cortex$medial_wall_mask$left
  cortexR_mwall <- xii$meta$cortex$medial_wall_mask$right
  # subcortVol <- xii$data$subcort
  # subcortLabs <- xii$meta$subcort$labels
  # subcortMask <- xii$meta$subcort$mask
  surfL <- xii$surf$cortex_left
  surfR <- xii$surf$cortex_right
  
  xii_new <- as.xifti(
    cortexL=as.vector(df[match(xii$data$cortex_left, df[,match_val], NA), map_val]),   
    cortexL_mwall=cortexL_mwall,
    cortexR= as.vector(df[match(xii$data$cortex_right, df[,match_val], NA), map_val]) ,
    cortexR_mwall=cortexR_mwall,
    #subcortVol=subcortVol, subcortLabs=subcortLabs,
    #subcortMask=subcortMask,
    surfL=surfL, surfR=surfR
  )
  if (save == TRUE) {
    write_cifti(xii_new, paste(df_name, map_val, format(Sys.time(), "%m_%d_%y"), ".dscalar.nii", sep = "_"))
  }
  
  #this doesn't actually work, as in viewer doesn't show image :(
  if (view == TRUE) {
    if(is.null(zlim)) {
      z_min <- min(df[,map_val])
      z_max <- max(df[,map_val])
      zlim <- c(z_min, z_max)
    }
    colors <- make_color_pal(colors, color_mode = color_mode, zlim = zlim)
    view_xifti_surface(xii_new, colors = colors, legend_embed = TRUE, zlim = zlim, title = paste(df_name, map_val, 
                                                                                                 sep   =  "_"))
  }
  return(xii_new)
}

`%!in%` = Negate(`%in%`)

