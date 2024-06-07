library(dplyr) # summary of concordance indices
library(tibble) # printing dataframe full
library(SurvMetrics) # BS, IBS
library(survminer) # visualizing survival analysis results
library(survival) # cox proportional hazardous model
library(caret) # stratified cross validation
library(randomForestSRC) # random survival forest
library(gbm) # gradient boosting model
library(pec) # calibration
library(survivalROC) # ROC curve
library(PerformanceAnalytics) # correlation matrix among continuous variables
library(survAUC) # Concordance index by Uno et al.
library(timeROC) # time dependent AUC comparison


## input and manipulate data
df <- read.csv("../codes/baseline_fibroscan_median_imputed.csv")
fcols <- c("sex", "dm", "antivirals", "cirrhosis", "hbeag")
indata <- df[df$Cohort==0, ]
indata[fcols] <- lapply(indata[fcols], factor)
print(dim(indata))
exdata <- df[(df$Cohort!=0) & (df$antivirals!=2), ]
exdata[fcols] <- lapply(exdata[fcols], factor)
print(dim(exdata))

## correlation matrix for continuous variables
pdf(file='correlation_maxtirx.pdf', width=9.7, height=9.7)
chart.Correlation(indata[, c('age', 'bilirubin', 'albumin', 'plt', 'ast', 'alt', 'loghbvdna', 
			     'Fibroscan', 'cap')], histogram=TRUE, pch=19)
dev.off()


## model training in cross validation setting of internal dataset
set.seed(20220902)
fold <- 5
cv_index <- createFolds(as.factor(indata$hcc), fold, returnTrain=TRUE)
c_idxs <- data.frame()
for (i in 1:length(cv_index)) {
	
    train_data <- indata[cv_index[[i]], ]
    test_data <- indata[-cv_index[[i]], ]

    ## cox proportional hazardous model
    cox <- coxph(Surv(time=hccfreesurvival, event=hcc) ~ age + sex + dm + antivirals + cirrhosis + 
                                                         bilirubin + albumin + plt + #ast + 
                                                         alt + loghbvdna + hbeag + log(Fibroscan+1) + cap,
  	    	  	 data=train_data, x=TRUE, y=TRUE)
    ### internal evaluation
    lp.new.cox <- predict(cox, newdata=test_data)
    surv.rsp.cox <- with(train_data, Surv(time=hccfreesurvival, event=hcc))
    surv.rsp.new.cox <- with(test_data, Surv(time=hccfreesurvival, event=hcc))


    ## random survival forest
    rsf <- rfsrc(Surv(time=hccfreesurvival, event=hcc) ~ age + sex + dm + antivirals + cirrhosis + 
                                                          bilirubin + albumin + plt + #ast + 
                                                          alt + loghbvdna + hbeag + log(Fibroscan+1) + cap, 
                  data=train_data, ntrees=1000,  seed=20220902)
    lp.new.rsf <- predict(rsf, newdata=test_data)
    
    ## mpageb scoring
    years <- c(3, 5, 6, 7)
    for (yr in years) {
        c_idxs <- rbind(c_idxs, c("mpageb", "in", yr, 
                                  UnoC(surv.rsp.cox, surv.rsp.new.cox, test_data$mpageb, time=365.24*yr)))
        c_idxs <- rbind(c_idxs, c("cox", "in", yr, 
                                  UnoC(surv.rsp.cox, surv.rsp.new.cox, lp.new.cox, time=365.24*yr)))
        c_idxs <- rbind(c_idxs, c("rsf", "in", yr, 
                                  UnoC(surv.rsp.cox, surv.rsp.new.cox, lp.new.rsf$predicted, time=365.24*yr)))
    }
    #cindex_in_rsf[[i]] <-sprintf("C-index in Uno et al. of RSF = %f", 
    #                             UnoC(surv.rsp.rsf, surv.rsp.new.rsf, lp.new.rsf$predicted, time=365.24*7))

    ## gradient boosting survival model
    #gbs <- gbm(Surv(time=hccfreesurvival, event=hcc) ~ age + sex + dm + antivirals + cirrhosis +
    #                                                   bilirubin + albumin + plt + #ast + 
    #                                                   alt + loghbvdna + hbeag + Fibroscan + cap,
    #           distribution='coxph', data=train_data, n.trees=2000, n.cores=15, verbose=0)
    #lp.new.gbs <- predict(gbs, newdata=test_data)
    #surv.rsp.gbs <- with(train_data, Surv(time=hccfreesurvival, event=hcc))
    #surv.rsp.new.gbs <- with(test_data, Surv(time=hccfreesurvival, event=hcc))
    #cindex_in_gbs[[i]] <- sprintf("C-index in Uno et al. of GBSM = %f", 
    #                              UnoC(surv.rsp.gbs, surv.rsp.new.gbs, lp.new.gbs, time=365.24*7))

       
}
colnames(c_idxs) <- c("models", "validation_type", "target_years", "C_indices")


## model training with internal dataset and evaluate the model with external dataset
### cox proportional hazardous model
cox <- coxph(Surv(time=hccfreesurvival, event=hcc) ~ age + sex + dm + antivirals + cirrhosis + 
                                                     bilirubin + albumin + plt + #ast + 
                                                     alt + loghbvdna + hbeag + log(Fibroscan+1) + cap,
    	  	 data=indata, x=TRUE, y=TRUE)
### external evaluation
lp.cox <- predict(cox)
lp.new.cox <- predict(cox, newdata=exdata)
surv.rsp <- with(indata, Surv(time=hccfreesurvival, event=hcc))
surv.rsp.new <- with(exdata, Surv(time=hccfreesurvival, event=hcc))


## random survival forest
rsf <- rfsrc(Surv(time=hccfreesurvival, event=hcc) ~ age + sex + dm + antivirals + cirrhosis + 
                                                      bilirubin + albumin + plt + #ast + 
                                                      alt + loghbvdna + hbeag + log(Fibroscan+1) + cap, 
              data=indata, ntrees=1000,  seed=20220902)
lp.rsf <- predict(rsf)
lp.new.rsf <- predict(rsf, newdata=exdata)
for (yr in years) {
    c_idxs <- rbind(c_idxs, c("rsf", "ex", yr, 
                              UnoC(surv.rsp, surv.rsp.new, lp.new.rsf$predicted, time=365.24*yr)))
    c_idxs <- rbind(c_idxs, c("mpageb", "ex", yr, 
                              UnoC(surv.rsp, surv.rsp.new, exdata$mpageb, time=365.24*yr)))
    c_idxs <- rbind(c_idxs, c("cox", "ex", yr, 
                              UnoC(surv.rsp, surv.rsp.new, lp.new.cox, time=365.24*yr)))
}

## printing out the summary concordance index
c_idxs[, 'C_indices'] = sapply(c_idxs[, 'C_indices'], as.numeric)
c_idx_bymodel_byvalidtype_byyear <- c_idxs %>%
  group_by(models, validation_type, target_years) %>%
  summarize(mean_cindices = mean(C_indices), 
            sd_cindices = sd(C_indices))

print(as_tibble(c_idx_bymodel_byvalidtype_byyear), n=30)

## calibration
calibrated3 <- calPlot(list("COX"=cox, "RSF"=rsf), time=365.24*3, data=exdata, type='risk', 
                      B=10, q=5, xlim=c(0, 0.1), ylim=c(0, 0.1), cex=0.8, percent=TRUE, col=c('blue', 'red'), 
                      lwd=c(1.5, 1.5), lty=c(2, 3))
calibrated5 <- calPlot(list("COX"=cox, "RSF"=rsf), time=365.24*5, data=exdata, type='risk', 
                      B=10, q=5, xlim=c(0, 0.4), ylim=c(0, 0.4), cex=0.8, percent=TRUE, col=c('blue', 'red'), 
                      lwd=c(1.5, 1.5), lty=c(2, 3))
calibrated6 <- calPlot(list("COX"=cox, "RSF"=rsf), time=365.24*6, data=exdata, type='risk', 
                      B=10, q=5, xlim=c(0, 0.4), ylim=c(0, 0.4), cex=0.8, percent=TRUE, col=c('blue', 'red'), 
                      lwd=c(1.5, 1.5), lty=c(2, 3))
calibrated7 <- calPlot(list("COX"=cox, "RSF"=rsf), time=365.24*7, data=exdata, type='risk', 
                      B=10, q=5, xlim=c(0, 0.4), ylim=c(0, 0.4), cex=0.8, percent=TRUE, col=c('blue', 'red'), 
                      lwd=c(1.5, 1.5), lty=c(2, 3))
#print(calibrated)
calplot_name <- sprintf("plot_calibration_exdata.pdf")
#pdf(file=calplot_name)
pdf(file=calplot_name)
par(mfrow=c(2,2))
plot(calibrated3)
mtext("3 years", side=3, line=1, cex=1, adj=-0.05)
plot(calibrated5)
mtext("5 years", side=3, line=1, cex=1, adj=-0.05)
plot(calibrated6)
mtext("6 years", side=3, line=1, cex=1, adj=-0.05)
plot(calibrated7)
mtext("7 years", side=3, line=1, cex=1, adj=-0.05)
dev.off()

## plot function of ROC and AUC

fun.survivalROC <- function(lp1, lp2, lp3, t) {
    res1 <- with(exdata,
                survivalROC(Stime        = hccfreesurvival,
                            status       = hcc,
                            marker       = get(lp1),
                            predict.time = t,
                            method       = "KM"))       # KM method without smoothing
    res2 <- with(exdata,
                survivalROC(Stime        = hccfreesurvival,
                            status       = hcc,
                            marker       = get(lp2),
                            predict.time = t,
                            method       = "KM"))       # KM method without smoothing
    res3 <- with(exdata,
                survivalROC(Stime        = hccfreesurvival,
                            status       = hcc,
                            marker       = get(lp3),
                            predict.time = t,
                            method       = "KM"))       # KM method without smoothing


    ## Plot ROCs
    with(res1, plot(TP ~ 1-FP, type="l", lty=1, lwd=1, col='black', xlab='1-Specificity', ylab='Sensitivity', 
                   main=sprintf("ROC, t = %.0f years", t/365.24)))
    with(res2, lines(TP ~ 1-FP, type='l', lty=2, lwd=1, col='blue'))
#                   main = sprintf("cox, t = %.0f, AUC = %.2f", t, AUC)))
    with(res3, lines(TP ~ 1-FP, type='l', lty=3, lwd=1, col='red')) 
#                   main = sprintf("rsf, t = %.0f, AUC = %.2f", t, AUC)))
    legend("bottomright", legend = c(sprintf('AUROC (mPAGE-B) = %.3f', res1$AUC), 
                                     sprintf('AUROC (Cox) = %.3f', res2$AUC),
                                     sprintf('AUROC (RSF) = %.3f', res3$AUC)), 
           lty = c(1, 2, 3), lwd=c(1, 1, 1), col = c('black', 'blue', 'red'), bty="n")
    abline(a = 0, b = 1, lty = 4)
}

## gradient boosting survival model
#gbs <- gbm(Surv(time=hccfreesurvival, event=hcc) ~ age + sex + dm + antivirals + cirrhosis +
#                                                   bilirubin + albumin + plt + #ast + 
#                                                   alt + loghbvdna + hbeag + Fibroscan + cap,
#           distribution='coxph', data=indata, n.trees=2000, n.cores=15, verbose=0)
#lp.new.gbs <- predict(gbs, newdata=exdata)
#surv.rsp.gbs <- with(indata, Surv(time=hccfreesurvival, event=hcc))
#surv.rsp.new.gbs <- with(exdata, Surv(time=hccfreesurvival, event=hcc))
#cindex_ex_gbs <- sprintf("C-index in Uno et al. of GBSM = %f", 
#                         UnoC(surv.rsp.gbs, surv.rsp.new.gbs, lp.new.gbs, time=365.24*7))

## computing Brier score of mPAGE-B, Cox and RSF at the 7th year.
mat_rsf <- predict(rsf, exdata)$survival
dis_time <- rsf$time.interest
surv_obj <- Surv(exdata$hccfreesurvival, exdata$hcc)
mat_cox <- predictSurvProb(cox, exdata, dis_time)
max_mpageb <- max(exdata$mpageb)

t_index <- 68
t_star <- 1096
brier_mpageb3 <- Brier(surv_obj, pre_sp=1-(exdata$mpageb/max_mpageb), t_star)
brier_cox3 <- Brier(surv_obj, pre_sp=mat_cox[, t_index], t_star)
brier_rsf3 <- Brier(surv_obj, pre_sp=mat_rsf[, t_index], t_star)

t_index <- 117
t_star <- 1764
brier_mpageb5 <- Brier(surv_obj, pre_sp=1-(exdata$mpageb/max_mpageb), t_star)
brier_cox5 <- Brier(surv_obj, pre_sp=mat_cox[, t_index], t_star)
brier_rsf5 <- Brier(surv_obj, pre_sp=mat_rsf[, t_index], t_star)

t_index <- 130
t_star <- 2181
brier_mpageb6 <- Brier(surv_obj, pre_sp=1-(exdata$mpageb/max_mpageb), t_star)
brier_cox6 <- Brier(surv_obj, pre_sp=mat_cox[, t_index], t_star)
brier_rsf6 <- Brier(surv_obj, pre_sp=mat_rsf[, t_index], t_star)

t_index <- 139
t_star <- 2544
brier_mpageb7 <- Brier(surv_obj, pre_sp=1-(exdata$mpageb/max_mpageb), t_star)
brier_cox7 <- Brier(surv_obj, pre_sp=mat_cox[, t_index], t_star)
brier_rsf7 <- Brier(surv_obj, pre_sp=mat_rsf[, t_index], t_star)

Brier_df <- data.frame('model'=c('mPAGE-B', 'Cox', 'RSF'), 
                       'BS_3yrs'=c(brier_mpageb3, brier_cox3, brier_rsf3),
                       'BS_5yrs'=c(brier_mpageb5, brier_cox5, brier_rsf5),
                       'BS_6yrs'=c(brier_mpageb6, brier_cox6, brier_rsf6),
                       'BS_7yrs'=c(brier_mpageb7, brier_cox7, brier_rsf7))
print(Brier_df)

## computing Integrated Brier score of Cox and BSF.
#ibs_mpageb <- IBS(surv_obj, sp_matrix=1-(exdata$mpageb/max_mpageb), dis_time)
ibs_cox <- IBS(surv_obj, sp_matrix=mat_cox, dis_time)
ibs_rsf <- IBS(surv_obj, sp_matrix=mat_rsf, dis_time)
ibs_df <- data.frame('model'=c('Cox', 'RSF'), 
                       'IBS'=c(ibs_cox, ibs_rsf))
print(ibs_df)


## printing ROC curve when 3, 5, 6, 7 years
exdata$lp.cox <- predict(cox, newdata=exdata, type='lp')
exdata$lp.rsf <- predict(rsf, newdata=exdata, type='lp')$predicted
#png(file="time_dependent_ROC.png", width=8.3, height=8.3)
pdf(file="time_dependent_ROC.pdf", width=8.3, height=8.3)
layout(matrix(1:4, nrow=2, ncol=2, byrow=TRUE))
res.survivalROC <- lapply(c(365.24*3, 365.24*5, 365.24*6, 365.24*7), function(t) {
                   fun.survivalROC(lp1='mpageb', lp2="lp.cox", lp3="lp.rsf", t)
})
dev.off()

## model comparison
roc.mpageb <- timeROC(T=exdata$hccfreesurvival, 
                      delta=exdata$hcc,
                      marker=exdata$mpageb,
                      cause=1,
                      weighting='marginal', 
                      time = c(365.24*3, 365.24*5, 365.24*6, 365.24*7),
                      #time = c(365.24*6, 365.24*7),
                      iid=TRUE)

roc.cox <- timeROC(T=exdata$hccfreesurvival, 
                   delta=exdata$hcc,
                   marker=exdata$lp.cox,
                   cause=1,
                   weighting='marginal', 
                   time = c(365.24*3, 365.24*5, 365.24*6, 365.24*7),
                   #time = c(365.24*6, 365.24*7),
                   iid=TRUE)

roc.rsf <- timeROC(T=exdata$hccfreesurvival, 
                   delta=exdata$hcc,
                   marker=exdata$lp.rsf,
                   cause=1,
                   weighting='marginal', 
                   time = c(365.24*3, 365.24*5, 365.24*6, 365.24*7),
                   #time = c(365.24*6, 365.24*7),
                   iid=TRUE)
print(compare(roc.mpageb, roc.cox, adjusted=TRUE))
print(compare(roc.mpageb, roc.rsf, adjusted=TRUE))


# Survival curve from HCC event of a patient assuming that the values of other variables 
# have the average values of validation data fro three values of Fibroscan = 0, 11, and 75.
 
dummy <- expand.grid(age=48.8, sex=as.factor(1), dm=as.factor(0), 
                     antivirals=as.factor(1), cirrhosis=as.factor(1), bilirubin=1.17,
                     albumin=4.12, plt=173, alt=121, loghbvdna=4.6, hbeag=as.factor(0), cap=214, 
                     Fibroscan=c(0, 11, 75))
print(dummy)
surv1 <- survfit(cox, newdata=dummy)
png(file='survival_probability_fibroscan.png', unit='px', width=720, height=480, bg='white')
plot(surv1, col=c('black', 'red', 'blue'), lty=1:3, xlab="Days", ylab="Survival probability from HCC", cex=1.2, cex.axis=1.2)
legend("bottomleft", col=c('black', 'red', 'blue'), 
       legend=paste(c('Fibroscan=0', 'Fibroscan=11', 'Fibroscan=75')), 
       lty=1:3, cex=1.2, text.font=4, bty='n')
#ggsurvplot(survfit(cox, data=exdata), palatte='#ffa500', ggtheme = theme_minimal())
dev.off()


png(file='hazard_ratio_indata.png', unit='px', width=960, height=640, bg='white')
ggforest(cox, data=indata, main='Hazard ratio of the Cox model trained by cohort with status of 0', fontsize=1)
dev.off()

png(file='hazard_ratio_exdata.png', unit='px', width=960, height=640, bg='white')
ggforest(cox, data=exdata, main='Hazard ratio of the Cox model validated by cohort with status of not equal to 0', fontsize=1)
dev.off()

