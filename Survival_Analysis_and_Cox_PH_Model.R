#Sadaf Saaber 
#Survial Analysis with Veteran Dataset in survival package 

############################
#Section 0
#DatasetDecription

##################################
#Description: Randomised trial of two treatment regimens for lung cancer. This is a standard survival analysis data set.
#Variable Descriptions 
#trt:	1=standard 2=test
#celltype:	1=squamous, 2=smallcell, 3=adeno, 4=large
#time:	survival time
#status:	censoring status
#karno:	Karnofsky performance score (100=good)
#diagtime:	months from diagnosis to randomisation
#age:	in years
#prior:	prior therapy 0=no, 10=yes
#Source: https://stat.ethz.ch/R-manual/R-devel/library/survival/html/veteran.html

#############################
#Section 1 
#installing packages
#############################
library(survival)
library(ggplot2)
library(dplyr)
library (psych)
library( Hmisc)

#############################
#Section 2
#loading and exploring the data 
#(including some linear regression analysis)
#############################

print(veteran) 

#one variable at a time using summary option 
summary(veteran$celltype)
summary(veteran$time)
summary (veteran$status)
summary(veteran$karno)
summary(veteran$diagtime)
summary(veteran$age)
summary(veteran$prior)

#also note that the following are numeric variables: time, status, karno, diagtime, age, prior 
#also note that the following is a factor variables: celltype 

#looking at the summary statistics by celltype 
describeBy(veteran, veteran$celltype)

#determing correlations between the variables 
#need to remove the celltype column to create a table of correlations between the numeric values
veteran2<-veteran[-c(2)]
cor(veteran2)


#creating x,y plots for the numeric variables 
xyplot(veteran$time~veteran$karno, data=veteran, xlab="time", ylab="karno", main= "time vs. karno")
xyplot(veteran$time~veteran$diagtime, data=veteran, xlab="time", ylab="diagtime", main= "time vs. diagntime")
xyplot(veteran$time~veteran$age, data=veteran, xlab="time", ylab="age", main= "time vs. age")
xyplot(veteran$karno~veteran$diagtime, data=veteran, xlab="karno", ylab="diagtime", main= "karno vs. diagtime")
xyplot(veteran$diagtime~veteran$age, data=veteran, xlab="diagtime", ylab="age", main= "diagtime vs age")

#calculating simple linear regressions for the numeric values mentioned above 
timekarno= lm(veteran$time~veteran$karno, data=veteran)
summary(timekarno)

timediagtime=lm(veteran$time~veteran$diagtime, data=veteran)
summary(timediagtime)

timeage=lm(veteran$time~veteran$age, data=veteran)
summary(timeage)

karnodiagtime=lm(veteran$karno~veteran$diagtime, data=veteran)
summary(karnodiagtime)

diagtimeage=lm(veteran$diagtime~veteran$age, data=veteran)
summary(diagtimeage)


#############################
#Section 3 
# Survival object and survfit 
#initial Kaplan-Meier survival estimates
#############################

#Creating a surival object linking the time and status variables together
km<-with(veteran, Surv(time, status))
head(km,100)

#surival analysis with the whole dataset for time being increiments of 10 
km_fit<-survfit(Surv(time, status)~1, data=veteran)
print(km_fit)
summary(km_fit, times = 400)
summary(km_fit, times = c(0,100,200,300,400,500,600,700,800,900,1000))

summary(km_fit, times=c(1, 10, 10*(1:10)))
plot (km_fit)

#surival analysis with the whole dataset for the time incriments being 20 
summary(km_fit, times=c(1, 20*(1:10)))
plot (km_fit)

#stratified by treatment group 
km_trt_fit<-survfit(Surv(time, status)~trt, data=veteran)
print(km_trt_fit)

summary(km_trt_fit, times = c(0,100,200,300,400,500,600,700,800,900,1000))
plot(km_trt_fit)

plot (km_trt_fit,lyt=c(1, 2), col=c("black", "grey"), 
      xlab = "survival time", ylab = "survival probabilities")
###########################
#section 4
#log rank test and the variations
#checking the PH assumption using graphical approaches 
##########################

#log rank test 
survdiff(Surv(time, status)~trt, data=veteran)

#potential variation of the log rank test which 
#(this is the equivalent to the Peto Peto modification of the Gehan-Wilcoxon test) 
survdiff(Surv(time, status)~trt, data=veteran, rho = 1)

#log rank test stratified by the variable of cell type 
survdiff(Surv(time, status)~trt+strata(celltype), data = veteran)

#log-log surival curves 
plot (km_trt_fit, fun="cloglog", xlab="time using logarithmic scale"
      , ylab = "log-log survival", main="log-log curves by treatment")

#log-log survival curves without having the time in the log scale
km_3=summary(km_trt_fit)
print (km_3)

km4=data.frame(km_3$strata, km_3$time, km_3$surv) 
names(km4)=c("clinic", "time", "surivival")

treatment1=km4[km4$clinic=="trt=1",]
treatment2=km4[km4$clinic=="trt=2",]

plot(treatment1$time, log(-log(treatment1$surivival)), xlab = "survival time"
     , ylab = "log-log surival", xlim = c(0, 1000), col="red", 
     lty="solid", main = "log-log curves by clinic")
par(new=T)
plot(treatment2$time, log (-log(treatment2$surivival)), axes=F, 
     xlab = "survival time",ylab = "", col="grey50", lty="dashed")
legend("bottomright", c("Treatment 1", "Treatment 2"),
       lty=c("solid", "dashed"), col=c("black", "grey50"))
par(new=F)

###########################
#Section 5
#Cox porportional Hazard Model
##############################
cox<-coxph(Surv(time, status)~trt+karno+diagtime+ age+prior, data=veteran)
summary(cox)
cox_fit<-survfit(cox)
plot(cox_fit)

#with other mefron, breslow and exact methods 
cox<-coxph(Surv(time, status)~trt+karno+diagtime+ age+prior, data=veteran, method = "efron")
summary(cox)

cox<-coxph(Surv(time, status)~trt+karno+diagtime+ age+prior, data=veteran, method = "breslow")
summary(cox)

cox<-coxph(Surv(time, status)~trt+karno+diagtime+ age+prior, data=veteran, method = "exact")
summary(cox)

#changing a factor variable back to numeric value 
is.factor(veteran$celltype)
veteran4<-as.numeric(veteran$celltype)
veteran$celltype2<-veteran4

#cox model with the factor value as numeric value with the different methods
cox2<-coxph(Surv(time, status)~trt+karno+diagtime+age+prior+celltype2, data=veteran)
summary(cox2)

cox2<-coxph(Surv(time, status)~trt+karno+diagtime+age+prior+celltype2, data=veteran, method = "efron")
summary(cox2)

cox2<-coxph(Surv(time, status)~trt+karno+diagtime+age+prior+celltype2, data=veteran, method = "breslow")
summary(cox2)

cox2<-coxph(Surv(time, status)~trt+karno+diagtime+age+prior+celltype2, data=veteran, method = "exact")
summary(cox2)

#############################
#Section 6 
#Stratified Cox Model for cell type by different methods 
#############################

coxstrat<-coxph(Surv(time, status)~trt+karno+diagtime+age+prior+strata(celltype) , data=veteran)
summary(coxstrat)

coxstrat<-coxph(Surv(time, status)~trt+karno+diagtime+age+prior+strata(celltype) , data=veteran, method="efron")
summary(coxstrat)

coxstrat<-coxph(Surv(time, status)~trt+karno+diagtime+age+prior+strata(celltype) , data=veteran, method = "breslow")
summary(coxstrat)

coxstrat<-coxph(Surv(time, status)~trt+karno+diagtime+age+prior+strata(celltype) , data=veteran, method = "exact")
summary(coxstrat)

