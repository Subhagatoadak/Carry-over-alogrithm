start_time <- Sys.time()

library(data.table)
library(dplyr)

#View(t)
# normalize <-function(x)({
#   normalizing_factor   = sum(x)
#   normalized <- x/normalizing_factor
#   return (normalized)
#   
# })

lagpad <- function(x, k) {
  if (k>0) {
    return (c(rep(NA, k), x)[1 : length(x)] );
  }
  else {
    return (c(x[(-k+1) : length(x)], rep(NA, -k)));
  }
}

ztran <- function(x, na.rm = TRUE) {
  mns <- colMeans(x, na.rm = na.rm)
  sds <- apply(x, 2, sd, na.rm = na.rm)
  x <- sweep(x, 2, mns, "-")
  x <- sweep(x, 2, sds, "/")
  x
}



library(readxl)
med_data <- read_excel("~/med_data_1.xlsx")
#EQ1<- readRDS("COstcoEQ.RDS")
EQ <- readRDS("~/none/Models/EQ_RDS.RDS")
#EQ2 <- readRDS("Amz_EQ.RDS")
rm(EQ2)
EQ$ZIP=NULL
EQ=EQ[EQ$Date>"2019-01-05",]
#rm(EQ1)
#EQ$EQ<-NULL
#EQ$PROJ<-NULL
#EQ$Date <- as.Date(EQ$Date)
dates = as.Date(unique(med_data$Week))
str(dates)
library(dplyr)
EQ = EQ %>% group_by(TDLINXID) %>% mutate(EQ_z = scale(EQ_proj))
med_data$Week=as.Date(med_data$Week)
#EQ1 = merge(EQ,med_data,by.x = "Date",by.y = "Week")
EQ$EQ_proj=NULL
EQ[is.na(EQ)]<-0
#View(EQ[EQ$TDLINXID==2657,])

EQ = EQ[complete.cases(EQ),]
#SM <- read_excel("Multiplier.xlsx",sheet = "SRTRM_MULT")
#LM <- read_excel("Multiplier.xlsx", sheet = "LNGTRM_MULT")
final_all_variable = data.frame()
for(u in 2:(length(med_data))){
  
  med=med_data[[u]]
 #  # ###################dweibull#############################################
 #  t = seq(1:104)
 #  #t1 = rep(0,104)
 #  #t = c(t,t1)
 #  t = as.data.frame(t)
 #  colnames(t)[1]="X"
 #  start_parameter1 =1
 #  start_parameter2 = 1
 #  end_parameter1 = 5
 #  end_parameter2 = 70
 #  for (i in c(seq(0,0.9,0.1),seq(start_parameter1,end_parameter1))){
 #  #for (i in c(seq(start_parameter1,end_parameter1))){
 #    for(j in start_parameter2:end_parameter2){
 #      list =c()
 #      for(k in 1:length(t$X)){
 #        t2= dweibull(t$X[k],i,j)
 #        list = append(list,t2)
 #      }
 #      t[[paste0("Weibull_",i,j)]]=list
 #    }
 #  }
 # 
 #  t1 = as.data.frame(matrix(0,104,length(t)))
 #  colnames(t1) =colnames(t)
 #  t= rbind(t,t1)
 #  DT= t
 #  library(data.table)
 #  for (j in 1:ncol(DT)) set(DT, which(is.infinite(DT[[j]])), j, 0)
 #  library(dplyr)
 #  DT = DT %>% select_if(~ !any(is.na(.)))
 #  DT[,] <- lapply(DT[,], function(x) x/sum(x))
 # 
 #  # mult=DT$Weibull_0.11
 #  # mult = as.data.frame(mult)
 #  # mult$Week1=NULL
 #  # mult$Week2=NULL
 #  # for(i in 1:208){
 #  #   mult[[paste0("Week_",i)]] = lagpad(mult$mult,i-1)
 #  # }
 #  # mult$mult=NULL
 #  # mult[is.na(mult)]=0
 #  # mult[,] <- lapply(mult[,], function(x) x/sum(x))
 # 
 # 
 #  #med=runif(208)
 #  #t5 = as.data.frame(t(t(mult)*med))
 #  #rowSums(t5)
 #  DT1_weibull = as.data.frame(DT$X)
 # # DT1_SM_Weibull = as.data.frame(DT$X)
 #  #DT1_LM_Weibull= as.data.frame(DT$X)
 #  for(m in 2:length(DT)){
 #    mult=DT[[m]]
 #    mult = as.data.frame(mult)
 #    mult$Week1=NULL
 #    mult$Week2=NULL
 #    for(i in 1:208){
 #      mult[[paste0("Week_",i)]] = lagpad(mult$mult,i-1)
 #    }
 #    mult$mult=NULL
 #    mult[is.na(mult)]=0
 #    # mult[,] <- lapply(mult[,], function(x) x/sum(x))
 # 
 # 
 # 
 #    t5 = as.data.frame(t(t(mult)*med))
 #    t5_sum=rowSums(t5)
 #    name1 = colnames(DT)[m]
 #    DT1_weibull[[name1]]=t5_sum
 # 
 # 
 #    # t5_SM =t5*SM
 #    # t5_SM_sum = rowSums(t5_SM)
 #    # DT1_SM_Weibull[[paste0(name1,"_SM")]]=t5_SM_sum
 #    #
 #    # t5_LM =t5*LM
 #    # t5_LM_sum = rowSums(t5_LM)
 #    # DT1_LM_Weibull[[paste0(name1,"_LM")]]=t5_LM_sum
 # 
 #    DT1_Z_weibull = DT1_weibull
 #    # DT1_Z_SM_Weibull = ztran(DT1_SM_Weibull)
 #    # DT1_Z_LM_Weibull= ztran(DT1_LM_Weibull)
 # 
 #  }
  # 
  # #eq = runif(208)*100
  # 
  # # df = data.frame(Variable = character(),Correl = double(),Pval = double())
  # # for(i in 2:length(DT1_weibull)){
  # #   cor1 = cor.test(DT1_weibull[[i]],eq)
  # #   t1 = as.data.frame(t(c(colnames(DT1_weibull)[i],cor1$estimate,cor1$p.value)))
  # #   colnames(t1)=c("Variable","Correl","Pval")
  # #   df = rbind(df,t1)
  # # }
  # # 
  # # final_cor_weibull=df[order(df$Correl,decreasing = TRUE),]
  # #########################################################################
  # 
  # 
  
  # ###################dlognormal#############################################
  # t = seq(1:104)
  # #t1 = rep(0,104)
  # #t = c(t,t1)
  # t = as.data.frame(t)
  # colnames(t)[1]="X"
  # start_parameter1 =1
  # start_parameter2 = 10
  # end_parameter1 = 20
  # end_parameter2 = 20
  # 
  # for (i in seq(start_parameter1,end_parameter1)){
  #   for(j in start_parameter2:end_parameter2){
  #     list =c()
  #     for(k in 1:length(t$X)){
  #       t2= dlnorm(t$X[k],i,j)
  #       list = append(list,t2)
  #     }
  #     t[[paste0("lnorm_",i,j)]]=list
  #   }
  # }
  # 
  # 
  # t1 = as.data.frame(matrix(0,104,length(t)))
  # colnames(t1) =colnames(t)
  # t= rbind(t,t1)
  # 
  # DT= t
  # library(data.table)
  # for (j in 1:ncol(DT)) set(DT, which(is.infinite(DT[[j]])), j, 0)
  # library(dplyr)
  # DT = DT %>% select_if(~ !any(is.na(.)))
  # DT[,] <- lapply(DT[,], function(x) x/sum(x))
  # # mult=DT$Weibull_0.11
  # # mult = as.data.frame(mult)
  # # mult$Week1=NULL
  # # mult$Week2=NULL
  # # for(i in 1:208){
  # #   mult[[paste0("Week_",i)]] = lagpad(mult$mult,i-1)
  # # }
  # # mult$mult=NULL
  # # mult[is.na(mult)]=0
  # # mult[,] <- lapply(mult[,], function(x) x/sum(x))
  # 
  # 
  # #med=runif(208)
  # #t5 = as.data.frame(t(t(mult)*med))
  # #rowSums(t5)
  # DT1_lnorm = as.data.frame(DT$X)
  # # DT1_SM_lnorm = as.data.frame(DT$X)
  # # DT1_LM_lnorm= as.data.frame(DT$X)
  # for(m in 2:length(DT)){
  #   mult=DT[[m]]
  #   mult = as.data.frame(mult)
  #   mult$Week1=NULL
  #   mult$Week2=NULL
  #   for(i in 1:208){
  #     mult[[paste0("Week_",i)]] = lagpad(mult$mult,i-1)
  #   }
  #   mult$mult=NULL
  #   mult[is.na(mult)]=0
  #   #mult[,] <- lapply(mult[,], function(x) x/sum(x))
  #   
  #   
  #   
  #   t5 = as.data.frame(t(t(mult)*med))
  #   t5_sum=rowSums(t5)
  #   name1 = colnames(DT)[m]
  #   DT1_lnorm[[name1]]=t5_sum
  #   
  #   # 
  #   # t5_SM =t5*SM
  #   # t5_SM_sum = rowSums(t5_SM)
  #   # DT1_SM_lnorm[[name1]]=t5_SM_sum
  #   # 
  #   # t5_LM =t5*LM
  #   # t5_LM_sum = rowSums(t5_LM)
  #   # DT1_LM_lnorm[[name1]]=t5_LM_sum
  #   
  #   DT1_Z_lnorm = DT1_lnorm
  #   # DT1_Z_SM_lnrom = ztran(DT1_SM_lnorm)
  #   # DT1_Z_LM_lnorm= ztran(DT1_LM_lnorm)
  #   # 
  #   
  # }
  # 
  # #eq = runif(208)*100
  # 
  # # df = data.frame(Variable = character(),Correl = double(),Pval = double())
  # # for(i in 2:length(DT1_lnorm)){
  # #   cor1 = cor.test(DT1_lnorm[[i]],eq)
  # #   t1 = as.data.frame(t(c(colnames(DT1_lnorm)[i],cor1$estimate,cor1$p.value)))
  # #   colnames(t1)=c("Variable","Correl","Pval")
  # #   df = rbind(df,t1)
  # # }
  # # 
  # # final_cor_lnorm=df[order(df$Correl,decreasing = TRUE),]
  # ######################################################################### 
  # 
  # ###################dlogis##########################################
  # t = seq(1:104)
  # #t1 = rep(0,104)
  # #t = c(t,t1)
  # t = as.data.frame(t)
  # colnames(t)[1]="X"
  # start_parameter1 =1
  # start_parameter2 = 1
  # end_parameter1 = 16
  # end_parameter2 = 16
  # 
  # for (i in seq(start_parameter1,end_parameter1)){
  #   for(j in start_parameter2:end_parameter2){
  #     list =c()
  #     for(k in 1:length(t$X)){
  #       t2= dlogis(t$X[k],i,j)
  #       list = append(list,t2)
  #     }
  #     t[[paste0("logis_",i,j)]]=list
  #   }
  # }
  # 
  # t1 = as.data.frame(matrix(0,104,length(t)))
  # colnames(t1) =colnames(t)
  # t= rbind(t,t1)
  # DT= t
  # library(data.table)
  # for (j in 1:ncol(DT)) set(DT, which(is.infinite(DT[[j]])), j, 0)
  # library(dplyr)
  # DT = DT %>% select_if(~ !any(is.na(.)))
  # DT[,] <- lapply(DT[,], function(x) x/sum(x))
  # 
  # # mult=DT$Weibull_0.11
  # # mult = as.data.frame(mult)
  # # mult$Week1=NULL
  # # mult$Week2=NULL
  # # for(i in 1:208){
  # #   mult[[paste0("Week_",i)]] = lagpad(mult$mult,i-1)
  # # }
  # # mult$mult=NULL
  # # mult[is.na(mult)]=0
  # # mult[,] <- lapply(mult[,], function(x) x/sum(x))
  # 
  # 
  # #med=runif(208)
  # #t5 = as.data.frame(t(t(mult)*med))
  # #rowSums(t5)
  # DT1_logis = as.data.frame(DT$X)
  # #DT1_SM_logis = as.data.frame(DT$X)
  # #DT1_LM_logis= as.data.frame(DT$X)
  # for(m in 2:length(DT)){
  #   mult=DT[[m]]
  #   mult = as.data.frame(mult)
  #   mult$Week1=NULL
  #   mult$Week2=NULL
  #   for(i in 1:208){
  #     mult[[paste0("Week_",i)]] = lagpad(mult$mult,i-1)
  #   }
  #   mult$mult=NULL
  #   mult[is.na(mult)]=0
  #   # mult[,] <- lapply(mult[,], function(x) x/sum(x))
  #   
  #   
  #   
  #   t5 = as.data.frame(t(t(mult)*med))
  #   t5_sum=rowSums(t5)
  #   name1 = colnames(DT)[m]
  #   DT1_logis[[name1]]=t5_sum
  #   
  #   
  #   
  #   # t5_SM =t5*SM
  #   # t5_SM_sum = rowSums(t5_SM)
  #   # DT1_SM_logis[[name1]]=t5_SM_sum
  #   # 
  #   # t5_LM =t5*LM
  #   # t5_LM_sum = rowSums(t5_LM)
  #   # DT1_LM_logis[[name1]]=t5_LM_sum
  #   
  #   DT1_Z_logis = DT1_logis
  #   # DT1_Z_SM_logis = ztran(DT1_SM_logis)
  #   # DT1_Z_LM_logis= ztran(DT1_LM_logis)
  #   
  #   
  # }
  # 
  # #eq = runif(208)*100
  # # 
  # # df = data.frame(Variable = character(),Correl = double(),Pval = double())
  # # for(i in 2:length(DT1_logis)){
  # #   cor1 = cor.test(DT1_logis[[i]],eq)
  # #   t1 = as.data.frame(t(c(colnames(DT1_logis)[i],cor1$estimate,cor1$p.value)))
  # #   colnames(t1)=c("Variable","Correl","Pval")
  # #   df = rbind(df,t1)
  # # }
  # # 
  # # final_cor_logis=df[order(df$Correl,decreasing = TRUE),]
  # #########################################################################
  # 
  # ###################NegExp#############################################
  # 
  # t = seq(1:104)
  # #t1 = rep(0,104)
  # #t = c(t,t1)
  # t = as.data.frame(t)
  # colnames(t)[1]="X"
  # start_parameter1 =0.1
  # end_parameter1 = 1
  # 
  # for (i in seq(start_parameter1,end_parameter1,0.1)){
  #   list =c()
  #   for(k in 1:length(t$X)){
  #     t2= dexp(t$X[k],i)
  #     list = append(list,t2)
  #   }
  #   t[[paste0("exp_",i)]]=list
  # }
  # 
  # t1 = as.data.frame(matrix(0,104,length(t)))
  # colnames(t1) =colnames(t)
  # t= rbind(t,t1)
  # 
  # DT= t
  # #DT[,] <- lapply(DT[,], function(x) x/sum(x))
  # library(data.table)
  # for (j in 1:ncol(DT)) set(DT, which(is.infinite(DT[[j]])), j, 0)
  # library(dplyr)
  # DT = DT %>% select_if(~ !any(is.na(.)))
  # DT[,] <- lapply(DT[,], function(x) x/sum(x))
  # 
  # # mult=DT$Weibull_0.11
  # # mult = as.data.frame(mult)
  # # mult$Week1=NULL
  # # mult$Week2=NULL
  # # for(i in 1:208){
  # #   mult[[paste0("Week_",i)]] = lagpad(mult$mult,i-1)
  # # }
  # # mult$mult=NULL
  # # mult[is.na(mult)]=0
  # # mult[,] <- lapply(mult[,], function(x) x/sum(x))
  # 
  # 
  # #med=runif(208)
  # #t5 = as.data.frame(t(t(mult)*med))
  # #rowSums(t5)
  # DT1_exp = as.data.frame(DT$X)
  # # DT1_SM_exp = as.data.frame(DT$X)
  # # DT1_LM_exp= as.data.frame(DT$X)
  # for(m in 2:length(DT)){
  #   mult=DT[[m]]
  #   mult = as.data.frame(mult)
  #   mult$Week1=NULL
  #   mult$Week2=NULL
  #   for(i in 1:208){
  #     mult[[paste0("Week_",i)]] = lagpad(mult$mult,i-1)
  #   }
  #   mult$mult=NULL
  #   mult[is.na(mult)]=0
  #   
  #   
  #   
  #   
  #   t5 = as.data.frame(t(t(mult)*med))
  #   t5_sum=rowSums(t5)
  #   name1 = colnames(DT)[m]
  #   DT1_exp[[name1]]=t5_sum
  #   
  #   
  #   # t5_SM =t5*SM
  #   # t5_SM_sum = rowSums(t5_SM)
  #   # DT1_SM_exp[[name1]]=t5_SM_sum
  #   # 
  #   # t5_LM =t5*LM
  #   # t5_LM_sum = rowSums(t5_LM)
  #   # DT1_LM_exp[[name1]]=t5_LM_sum
  #   
  #   DT1_Z_exp = DT1_exp
  #   # DT1_Z_SM_exp = ztran(DT1_SM_exp)
  #   # DT1_Z_LM_exp= ztran(DT1_LM_exp)
  # }
  # 
  # #eq = runif(208)*100
  # # 
  # # df = data.frame(Variable = character(),Correl = double(),Pval = double())
  # # for(i in 2:length(DT1_exp)){
  # #   cor1 = cor.test(DT1_exp[[i]],eq)
  # #   t1 = as.data.frame(t(c(colnames(DT1_exp)[i],cor1$estimate,cor1$p.value)))
  # #   colnames(t1)=c("Variable","Correl","Pval")
  # #   df = rbind(df,t1)
  # # }
  # # 
  # # final_cor_exp=df[order(df$Correl,decreasing = TRUE),]
  # #########################################################################
  # 
  # ###################chisq#############################################

  t = seq(1:104)
  #t1 = rep(0,104)
  #t = c(t,t1)
  t = as.data.frame(t)
  colnames(t)[1]="X"
  start_parameter1 =1
  end_parameter1 = 70


  for (i in seq(start_parameter1,end_parameter1,1)){
    list =c()
    for(k in 1:length(t$X)){
      t2= dchisq(t$X[k],i)
      list = append(list,t2)
    }
    t[[paste0("chisq_",i)]]=list
  }


  t1 = as.data.frame(matrix(0,104,length(t)))
  colnames(t1) =colnames(t)
  t= rbind(t,t1)
  DT= t
  library(data.table)
  for (j in 1:ncol(DT)) set(DT, which(is.infinite(DT[[j]])), j, 0)
  library(dplyr)
  DT = DT %>% select_if(~ !any(is.na(.)))
  DT[,] <- lapply(DT[,], function(x) x/sum(x))
  # mult=DT$Weibull_0.11
  # mult = as.data.frame(mult)
  # mult$Week1=NULL
  # mult$Week2=NULL
  # for(i in 1:208){
  #   mult[[paste0("Week_",i)]] = lagpad(mult$mult,i-1)
  # }
  # mult$mult=NULL
  # mult[is.na(mult)]=0
  # mult[,] <- lapply(mult[,], function(x) x/sum(x))


  #med=runif(208)
  #t5 = as.data.frame(t(t(mult)*med))
  #rowSums(t5)
  DT1_chi = as.data.frame(DT$X)
  # DT1_SM_chi = as.data.frame(DT$X)
  # DT1_LM_chi= as.data.frame(DT$X)
  for(m in 2:length(DT)){
    mult=DT[[m]]
    mult = as.data.frame(mult)
    mult$Week1=NULL
    mult$Week2=NULL
    for(i in 1:208){
      mult[[paste0("Week_",i)]] = lagpad(mult$mult,i-1)
    }
    mult$mult=NULL
    mult[is.na(mult)]=0
    # mult[,] <- lapply(mult[,], function(x) x/sum(x))



    t5 = as.data.frame(t(t(mult)*med))
    t5_sum=rowSums(t5)
    name1 = colnames(DT)[m]
    DT1_chi[[name1]]=t5_sum



    # t5_SM =t5*SM
    # t5_SM_sum = rowSums(t5_SM)
    # DT1_SM_chi[[name1]]=t5_SM_sum
    #
    # t5_LM =t5*LM
    # t5_LM_sum = rowSums(t5_LM)
    # DT1_LM_chi[[name1]]=t5_LM_sum
    #
    DT1_Z_chi = DT1_chi
    # DT1_Z_SM_chi = ztran(DT1_SM_chi)
    # DT1_Z_LM_chi= ztran(DT1_LM_chi)


  }
  
  #eq = runif(208)*100
  # 
  # df = data.frame(Variable = character(),Correl = double(),Pval = double())
  # for(i in 2:length(DT1_chi)){
  #   cor1 = cor.test(DT1_chi[[i]],eq)
  #   t1 = as.data.frame(t(c(colnames(DT1_chi)[i],cor1$estimate,cor1$p.value)))
  #   colnames(t1)=c("Variable","Correl","Pval")
  #   df = rbind(df,t1)
  # }
  # 
  # final_cor_chi=df[order(df$Correl,decreasing = TRUE),]
  #########################################################################
  #DT1_Z_weibull$Week = NULL
  # DT1_Z_SM_Weibull$Week =NULL
  # DT1_Z_LM_Weibull$Week =NULL
  # 
  # DT1_Z_lnorm$Week = NULL 
  # DT1_Z_SM_lnorm$Week =NULL
  # DT1_Z_LM_lnorm$Week =NULL
  # 
  # 
  # DT1_Z_logis$Week = NULL 
  # DT1_Z_SM_logis$Week =NULL
  # DT1_Z_LM_logis$Week =NULL
  # 
  # 
  # DT1_Z_exp$Week = NULL 
  # DT1_Z_SM_exp$Week =NULL
  # DT1_Z_LM_exp$Week =NULL
  # 
  DT1_Z_chi$Week = NULL 
  # DT1_Z_SM_chi$Week =NULL
  # DT1_Z_LM_chi$Week =NULL
  #############################################
  #DT1_Z_weibull$`DT$X` = NULL 
  # DT1_Z_SM_Weibull$`DT$X` =NULL
  # DT1_Z_LM_Weibull$`DT$X` =NULL
  
  #DT1_Z_lnorm$`DT$X` = NULL 
  # DT1_Z_SM_lnrom$`DT$X`=NULL
  # DT1_Z_LM_lnorm$`DT$X`=NULL
  # 
  
  #DT1_Z_logis$`DT$X` = NULL 
  # DT1_Z_SM_logis$`DT$X` =NULL
  # DT1_Z_LM_logis$`DT$X` =NULL
  
  
  #DT1_Z_exp$`DT$X` = NULL
  # DT1_Z_SM_exp$`DT$X` =NULL
  # DT1_Z_LM_exp$`DT$X` =NULL
  
  DT1_Z_chi$`DT$X` = NULL 
  # DT1_Z_SM_chi$`DT$X` =NULL
  # DT1_Z_LM_chi$`DT$X` =NULL
  # 
  
  # DT1=cbind(DT1_Z_weibull,DT1_Z_lnorm,DT1_Z_logis,DT1_Z_exp,DT1_Z_chi)
  # DT1_SM=cbind(DT1_Z_SM_Weibull,DT1_Z_SM_lnrom,DT1_Z_SM_logis,DT1_Z_SM_exp,DT1_Z_SM_chi)
  # DT1_LM = cbind(DT1_Z_LM_Weibull,DT1_Z_LM_lnorm,DT1_Z_LM_logis,DT1_Z_LM_exp,DT1_Z_LM_chi)
  
  #DT1_Z_weibull$Week = dates
  #DT1_Z_lnorm$Week = dates
  #DT1_Z_logis$Week = dates
  #DT1_Z_exp$Week = dates
  DT1_Z_chi$Week = dates
  # DT1_Z_SM_Weibull$Week = dates
  # DT1_Z_SM_lnrom$Week = dates
  # DT1_Z_SM_logis$Week = dates
  # DT1_Z_SM_exp$Week = dates
  # DT1_Z_SM_chi$Week = dates
  # DT1_Z_LM_Weibull$Week =dates
  # DT1_Z_LM_lnorm$Week =dates
  # DT1_Z_LM_logis$Week =dates
  # DT1_Z_LM_exp$Week =dates
  # DT1_Z_LM_chi$Week =dates
  # 
  #adstock =c("DT1_Z_lnorm","DT1_Z_logis","DT1_Z_exp","DT1_Z_chi")
  adstock =c("DT1_Z_chi")
  #adstock =c("DT1_Z_weibull")
  # adstock_SM =c("DT1_Z_SM_Weibull","DT1_Z_SM_lnrom","DT1_Z_SM_logis","DT1_Z_SM_exp","DT1_Z_SM_chi")
  # #adstock_SM =c("DT1_Z_SM_logis","DT1_Z_SM_exp")
  # adstock_LM =c("DT1_Z_LM_Weibull","DT1_Z_LM_lnorm","DT1_Z_LM_logis","DT1_Z_LM_exp","DT1_Z_LM_chi")
  # #adstock_LM =c("DT1_Z_LM_logis","DT1_Z_LM_exp")
  # #t<-eval(as.name(dd[1]))
  # #keep(DT1,DT1_LM,DT1_SM,EQ,med_data,final_all_variable,sure = TRUE)
  # 
  memory.limit(size=9999999999999)
  #l=adstock[1]
 #rm(t,t1,t5,mult,DT1_weibull,DT)
  ###########################################################################
  adstock1 = data.frame(Variable=character(),Correl=numeric(),Pval=numeric())
  for(l in adstock){
    DT1= merge(EQ,eval(as.name(l)),by.x="Date",by.y = "Week")
    #DT1 = DT1_Z[DT1_Z$Date>"2019-01-05",]
    
    
    #DT1 = DT1_Z
    #DT1 = DT1 %>% select_if(~ !any(is.na(.)))
    #DT1 = DT1_Z[complete.cases(DT1),]
    #scale1 <- function(x) scale(x)[,1]
    #jj = colnames(DT1)[5:length(DT1)]
    #dd=DT1
    #DTy <- DT1 %>%group_by(TDLINXID) %>%mutate_at(jj, ~(scale(.) %>% as.vector))
    #setDT(dd)[, (jj) := lapply(.SD, scale1) , by = .(TDLINXID), .SDcols = jj] 
    #head(data_scale2)                 # Head of scaled data
    df = data.frame(Variable = character(),Correl = double(),Pval = double())
    for(i in 5:length(DT1)){
      DTm=DT1[,c(1,2,4,i)]
      mm = colnames(DTm)[4]
      DTm = DTm %>% group_by(TDLINXID) %>% mutate("{mm}_z" := scale(eval(as.name(colnames(DTm)[4]))))
      DTm[is.na(DTm)] <- 0
      cor1 =cor.test(DTm[[5]],DTm$EQ_z)
      t1 = as.data.frame(t(c(colnames(DTm)[5],cor1$estimate,cor1$p.value)))
      colnames(t1)=c("Variable","Correl","Pval")
      df = rbind(df,t1)
      write.csv(df,paste0("df",colnames(med_data)[u],".csv"))
    }
    df[df=="NaN"]=NA
    df = df[complete.cases(df),]
    df$Correl<-as.numeric(df$Correl)
    df$Pval<-as.numeric(df$Pval)
    final_cor=df
    #final_cor = df[df$Correl>0,]
    #final_cor=final_cor[order(final_cor$Correl,decreasing = TRUE),]
    DT1_best = final_cor[final_cor$Correl==max(final_cor$Correl),]
    adstock1 = rbind(adstock1,DT1_best)
  }
  
  # ###########################################################################
  # adstock1_SM = data.frame(Variable=character(),Correl=numeric(),Pval=numeric())
  # for(l in adstock_SM){
  #   DT1_Z_SM= merge(EQ,eval(as.name(l)),by.x="Date",by.y = "Week")
  #   DT1_SM=DT1_Z_SM
  #   DT1_SM = DT1_SM %>% select_if(~ !any(is.na(.)))
  #   df_SM = data.frame(Variable = character(),Correl = double(),Pval = double())
  #   for(i in 5:length(DT1_SM)){
  #     cor1 =cor.test(DT1_SM[[i]],DT1_SM$EQ_z)
  #     t1 = as.data.frame(t(c(colnames(DT1_SM)[i],cor1$estimate,cor1$p.value)))
  #     colnames(t1)=c("Variable","Correl","Pval")
  #     df_SM = rbind(df_SM,t1)
  #   }
  #   df_SM[df_SM=="NaN"]=NA
  #   df_SM = df_SM[complete.cases(df_SM),]
  #   df_SM$Correl<-as.numeric(df_SM$Correl)
  #   df_SM$Pval<-as.numeric(df_SM$Pval)
  #   final_cor_SM = df_SM[df_SM$Correl>0,]
  #   final_cor_SM=final_cor_SM[order(final_cor_SM$Correl,decreasing = TRUE),]
  #   DT1_best_SM = final_cor_SM[final_cor_SM$Correl==max(final_cor_SM$Correl),]
  #   adstock1_SM = rbind(adstock1_SM,DT1_best_SM)
  # }
  # 
  # ########################################################################################
  # adstock1_LM = data.frame(Variable=character(),Correl=numeric(),Pval=numeric())
  # for(l in adstock_LM){
  #   DT1_Z_LM = merge(EQ,eval(as.name(l)),by.x="Date",by.y = "Week")
  #   DT1_LM=DT1_Z_LM
  #   #keep(DT1,DT1_LM,DT1_SM,EQ,med_data,final_all_variable,sure = TRUE)
  #   DT1_LM = DT1_LM %>% select_if(~ !any(is.na(.)))
  #   df_LM = data.frame(Variable = character(),Correl = double(),Pval = double())
  #   for(i in 5:length(DT1_LM)){
  #     cor1 =cor.test(DT1_LM[[i]],DT1_LM$EQ_z)
  #     t1 = as.data.frame(t(c(colnames(DT1_LM)[i],cor1$estimate,cor1$p.value)))
  #     colnames(t1)=c("Variable","Correl","Pval")
  #     df_LM = rbind(df_LM,t1)
  #   }
  #   df_LM[df_LM=="NaN"]=NA
  #   df_LM = df_LM[complete.cases(df_LM),]
  #   df_LM$Correl<-as.numeric(df_LM$Correl)
  #   df_LM$Pval<-as.numeric(df_LM$Pval)
  #   final_cor_LM = df_LM[df_LM$Correl>0,]
  #   final_cor_LM=final_cor_LM[order(final_cor_LM$Correl,decreasing = TRUE),]
  #   DT1_best_LM = final_cor_LM[final_cor_LM$Correl==max(final_cor_LM$Correl),]
  #   adstock1_LM = rbind(adstock1_LM,DT1_best_LM)
  # }
  ##################################################################################
  adstock1=adstock1[order(adstock1$Correl,decreasing = TRUE),]
  # adstock1_SM=adstock1_SM[order(adstock1_SM$Correl,decreasing = TRUE),]
  # adstock1_LM=adstock1_LM[order(adstock1_LM$Correl,decreasing = TRUE),]
  # 
  # dt_temp = as.data.frame(t(as.data.frame(c(colnames(med_data)[u],adstock1$Variable[1],adstock1_SM$Variable[1],adstock1_LM$Variable[1]))))
  # rownames(dt_temp)=NULL

  adstock1$Vehicle = colnames(med_data)[u]

  final_all_variable = rbind(final_all_variable,adstock1)
  
}
write.csv(final_all_variable,"adstock_new_terations.csv")
end_time <- Sys.time()
write.csv(df,"final_df.csv")
