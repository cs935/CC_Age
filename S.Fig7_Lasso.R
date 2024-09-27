library(survival)
library(survminer)
library(tidyverse)
library(glmnet)
library(survivalROC)
library(timeROC)

cl<-read.csv("clinicalData.csv",header = T,row.names = 1) 
cl <- filter(cl,time>0) 
exp <- read.csv("expdataTumor.csv",header = T,check.names = F,row.names = 1)
genelist <- c("CD276","LAMA3","LDHA","PVR","CXCL2","CXCL5","CXCL8","CXCR1")
exp <- t(exp[genelist,])
exp <- as.data.frame(exp)
sameSample <- intersect(rownames(cl),rownames(exp))
cl <- cl[sameSample,]
exp <- exp[sameSample,]
x <- as.matrix(exp)
y <- as.matrix(cl)
lasso <- glmnet(x, y, family = "cox")
plot(lasso,xvar = "lambda",label = TRUE) 
fitCV <- cv.glmnet(x, y,family = "cox")
plot(fitCV) 
fitCV$lambda.min
coefficient <- coef(fitCV, s = "lambda.min")
coefficient
Active.index <- which(as.numeric(coefficient) != 0)
Active.coefficient <- as.numeric(coefficient)[Active.index]
sig_gene_mult_cox <- rownames(coefficient)[Active.index]
exp1 <- exp[,sig_gene_mult_cox] 
expr.lasso_cox <- cbind(exp1,cl)
multiCox <- coxph(Surv(time, status) ~ ., data =  expr.lasso_cox)
summary(multiCox)
riskScore=predict(multiCox,type="risk",newdata=expr.lasso_cox)
riskScore<-as.data.frame(riskScore)
riskScore$sample <- rownames(riskScore)
riskScore_cli <- cbind(riskScore,cl)
cutoff<- 60 
SROC= survivalROC(Stime = riskScore_cli$time, 
                  status = riskScore_cli$status,     
                  marker = riskScore_cli$riskScore,                      
                  predict.time = cutoff, method= "KM" ) 
cut.op= SROC$cut.values[which.max(SROC$TP-SROC$FP)]  
riskScore_cli$riskScore2 = ifelse(riskScore_cli$riskScore >=cut.op,'high','low')   
range(riskScore_cli$riskScore)
fit <- survfit(Surv(time, as.numeric(status)) ~ riskScore2, data=riskScore_cli)
lasso_KM <- ggsurvplot(fit, data = riskScore_cli,                             
                       pval = T,                             
                       risk.table = T,                             
                       surv.median.line = "hv",                          
                       legend.labs=c("High risk","Low risk"),                             
                       legend.title="RiskScore",                            
                       ylab="Cumulative survival")
lasso_KM
with(riskScore_cli,     
     ROC_riskscore <<- timeROC(T = time,                               
                               delta = status,                               
                               marker = riskScore,                               
                               cause = 1,                               
                               weighting = "marginal",                               
                               times = c(12,36,60),                               
                               ROC = TRUE,                               
                               iid = TRUE))
plot(ROC_riskscore, time = 12, col = "red", add = F,title = "")
plot(ROC_riskscore, time = 36, col = "blue", add = T)
plot(ROC_riskscore, time = 60, col = "purple", add = T)
legend("bottomright",c("1-Year","3-Year","5-Year"),col=c("red","blue","purple"),lty=1,lwd=2)
text(0.5,0.2,paste("1-Year AUC = ",round(ROC_riskscore$AUC[1],3)))
text(0.5,0.15,paste("3-Year AUC = ",round(ROC_riskscore$AUC[2],3)))
text(0.5,0.1,paste("5-Year AUC = ",round(ROC_riskscore$AUC[3],3)))
cl.ext <-read.csv("GSE44001_clinical.Data.csv",header = T,row.names = 1) 
cl.ext <- filter(cl.ext,time>0) 
exp.ext <- read.csv("GSE44001_expData.csv",header = T,row.names = 1)
exp.ext <- t(exp.ext[sig_gene_mult_cox,])
exp.ext <- as.data.frame(exp.ext)
expr.lasso_cox.ext <- cbind(cl.ext,exp.ext )
riskScore=predict(multiCox,type="risk",newdata=expr.lasso_cox.ext)
riskScore<-as.data.frame(riskScore)
riskScore$sample <- rownames(riskScore)
riskScore_cli <- cbind(riskScore,cl.ext)
cutoff<- 60
SROC= survivalROC(Stime = riskScore_cli$time, 
                  status = riskScore_cli$status,     
                  marker = riskScore_cli$riskScore,                      
                  predict.time = cutoff, method= "KM" ) 
cut.op= SROC$cut.values[which.max(SROC$TP-SROC$FP)]  
riskScore_cli$riskScore2 = ifelse(riskScore_cli$riskScore >=cut.op,'high','low')   
fit <- survfit(Surv(time, as.numeric(status)) ~ riskScore2, data=riskScore_cli)
lasso_KM <- ggsurvplot(fit, data = riskScore_cli,                             
                       pval = T,                             
                       risk.table = T,                             
                       surv.median.line = "hv",                          
                       legend.labs=c("High risk","Low risk"),                             
                       legend.title="RiskScore",                            
                       ylab="Cumulative survival")
lasso_KM
with(riskScore_cli,     
     ROC_riskscore <<- timeROC(T = time,                               
                               delta = status,                               
                               marker = riskScore,                               
                               cause = 1,                               
                               weighting = "marginal",                               
                               times = c(12,36,60),                               
                               ROC = TRUE,                               
                               iid = TRUE))
plot(ROC_riskscore, time = 12, col = "red", add = F,title = "")
plot(ROC_riskscore, time = 36, col = "blue", add = T)
plot(ROC_riskscore, time = 60, col = "purple", add = T)
legend("bottomright",c("1-Year","3-Year","5-Year"),col=c("red","blue","purple"),lty=1,lwd=2)
text(0.5,0.2,paste("1-Year AUC = ",round(ROC_riskscore$AUC[1],3)))
text(0.5,0.15,paste("3-Year AUC = ",round(ROC_riskscore$AUC[2],3)))
text(0.5,0.1,paste("5-Year AUC = ",round(ROC_riskscore$AUC[3],3)))