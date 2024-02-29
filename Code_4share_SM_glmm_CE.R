library("nlme")
library("Matrix")
library("lme4")
library("readxl")
library("plyr")
library("dplyr")
library("stats")
library("car")
library("ggplot2")
library("xml2")
library("lmerTest")
library("tidyverse")
library("chron")
library('sjmisc')

rm(list=ls())

# Paths and settings

csvpath <- '' # <---- add your .csv path here !
csvfile <- paste(csvpath,'Fig2_Fig3.xlsx', sep='') # <---- add your .csv name here !

resultpath <- paste('', sep='') # <---- add your result path here !
dir.create(resultpath,showWarnings = FALSE,recursive = TRUE)

yourresultfilename <- ''; # <---- add your result filename here !

# settings
freq1 <- 'Delta1'
freq2 <- 'beta1'
nchan <- 2447

# comp1 = 1: CEWR, comp1 = 0: NE, comp1 = 2: CE
# For fig.2: comp1 <- 2 and comp2 <- 0;
# For fig.3: comp1 <- 2 and comp2 <- 1;

comp1 <- 2  # <---- change according to the model you want to run !
comp2 <- 0  # <---- change according to the model you want to run !

var2predict <- 'CE'

if (comp1 == 1){
  comp1_str = 'CEWR'
}else if (comp1 == 2){
  comp1_str = 'CE'
}else if (comp1 == 0){
  comp1_str = 'NE'
}

# reading excel file
allcol <- colnames(read_excel(csvfile))
allcoltype <- ifelse(allcol=='IDsub'| allcol=='duration'| allcol=='timeBeg','text','numeric')
data <- read_excel(csvfile1,col_types = allcoltype)
data <- data[c(which(data$P_H==1 & data$First_Unanimity_Episode==1)),]

#########################
namescol <- colnames(data[,-(1:4)])
namescol <- namescol[str_detect(namescol,'_sec_')]

# calculate log
for (c in namescol){
  data[paste(c,'_log',sep='')]<-log(data[,c])
  print(class(data[c]))
}

rm(namescol)
namescol <- colnames(data[,-(1:4)])

# reorganize data
dataCE2 <- data[which((data[,var2predict][[1]] == comp1)),]
dataCE0 <- data[which((data[,var2predict][[1]] == comp2)),]
dataCE <- rbind(dataCE2,dataCE0)
rm(dataCE2, dataCE0)


# factors
dataCE$CE <- as.factor(dataCE$CE)
dataCE$prov_spont <- as.factor(dataCE$prov_spont)


# detect existing timeframes
timeframelist <- namescol[str_detect(namescol,paste('Delta1_',sep='')) & str_detect(namescol,'log')]

# timeframe loop
for (t in timeframelist){
  
  # extracting frequency name
  t_freq1 <- paste(freq1,substr(t,nchar(freq1)+1,nchar(t)), sep = '')
  t_freq2 <- paste(freq2,substr(t,nchar(freq1)+1,nchar(t)), sep = '')
  timeframe <- substr(t,7,nchar(t))
  
  # setting formula
  formulastr <- paste(var2predict, '~ ', t_freq1,'+', t_freq2, '+ (1|IDsub) + (1|prov_spont)',sep='')
  main_stat_list = c(freq1,freq2)
  
  beta_val <- data.frame(matrix(ncol = length(main_stat_list), nrow = nchan))
  colnames(beta_val) <- lapply(main_stat_list,paste,"_beta",sep="")
  z_val <- data.frame(matrix(ncol = length(main_stat_list), nrow = nchan))
  colnames(z_val) <- lapply(main_stat_list,paste,"_z",sep="")
  beta_p_val <- data.frame(matrix(ncol = length(main_stat_list), nrow = nchan))
  colnames(beta_p_val) <- lapply(main_stat_list,paste,"_beta_P",sep="")
  
  glmform <- as.formula(formulastr)
  
  # voxel loop
  for(c in 1:nchan){
    
    # running model
    model <- glmer(glmform, dataCE[dataCE$chan==c,], family = binomial(link = 'logit'), control =  glmerControl(optimizer = 'bobyqa'), nAGQ = 0)
    sum_res <- summary(model)
    
    coef_model <- coef(sum_res)  
    beta_val[c,] <- coef_model[-1,1] # beta estimates
    z_val[c,] <- coef_model[-1,4] #
    beta_p_val[c,] <- coef_model[-1,5] # beta estimates p value 
    
    rm(model, sum_res, coef_model)
    
  }#end of channel loop
  
  # create a unique table
  stat_df <- cbind(beta_val,z_val,beta_p_val)
  
  save(stat_df, file = paste(resultpath, yourresultfilename, '_', timeframe,'_.RData', sep =''))
  rm(beta_val, z_val, beta_p_val, stat_df, t_freq1, t_freq2)
  
}#end of timeframe loop