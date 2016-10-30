rm(list = ls())
setwd(getwd())
library(plyr)

data<-read.csv("bank.csv", header = TRUE, sep = ";" )

bank1<-data[,sapply(data,is.factor)]

sum(as.numeric(!complete.cases(data))) #check missing values

#convert all integer variables into factor

bank1$y<-as.numeric(bank1$y=='yes')
n_predictors=(ncol(bank1)-1)

# Data grouping and counts generation -----------------------------------------------------------
bank=ddply(bank1,as.formula(paste("~",paste(names(bank1)[1:n_predictors], collapse=" + "),collapse=" ")),summarise,Freq=sum(y))
sum(bank$Freq)

bank_y=bank[,(n_predictors+1):ncol(bank)]
bank_x=bank[,1:n_predictors]

##check factor levels
str(bank_x)

#generate design matrix, options= maineffects, twoways,threeways
main_effects_terms <-c()
twoway_terms <-c()
threeway_terms <-c()

#report basic statistics  
mean(bank_y)
sum(bank_y!=0)



# # design matrix generation ----------------------------------------------

main_effects_terms <-paste(names(bank_x), collapse=" + ")

design1=model.matrix(as.formula(paste("~",main_effects_terms,collapse=" ")) ,bank_x)
ncol(design1)
ncol(design1)-sum(apply(design1==0,2,all)) #report number of non-zero columns
design1=design1[,which(!apply(design1,2,FUN = function(x){t(x)%*%bank_y==0}))]
ncol(design1)

twoway<-combn(names(bank_x), 2, paste, collapse="*")
twoway_terms <-paste(twoway, collapse=" + ")

design2=model.matrix(as.formula(paste("~",twoway_terms,collapse=" ")) ,bank_x)
ncol(design2)
ncol(design2)-sum(apply(design2==0,2,all)) #report number of non-zero columns
design2=design2[,which(!apply(design2,2,FUN = function(x){t(x)%*%bank_y==0}))]
ncol(design2)

threeway<-combn(names(bank_x), 3, paste, collapse="*")
threeway_terms <-paste(threeway, collapse=" + ")

design3=model.matrix(as.formula(paste("~",threeway_terms,collapse=" ")) ,bank_x)
ncol(design3)
ncol(design3)-sum(apply(design3==0,2,all)) #report number of non-zero columns
design3=design3[,which(!apply(design3,2,FUN = function(x){t(x)%*%bank_y==0}))]
ncol(design3)


#make last column as observation
options='threeway';
design=design3
design=cbind(design, bank_y) 


#write.table(design,file=paste(options, "_design_matrix.dat", sep=""),sep=",", col.names = F, row.names = F)



# L1 selected features ----------------------------------------------------------------------

# design=design1
# index=read.table("L1_feature_main.txt")


index2<-read.table("L1_EBIC_feature_twoway.txt")[,1]
length(index2)
variables2<-colnames(design2)[index2]


index3<-read.table("L1_EBIC_feature_threeway.txt")[,1]
length(index3)
variables3<-colnames(design3)[index3]

intersect(index2,index3)

share_variables<-intersect(variables2,variables3)
length(share_variables)
twoway_only_variables<-setdiff(variables2, share_variables)
threeway_only_variables<-setdiff(variables3, share_variables)


share_variables
twoway_only_variables
threeway_only_variables

# index2<-read.table("L1_feature_twoway_start.txt")[,1]
# length(index2)
# twoway_start<-colnames(design2)[index2][1:16]
# 
# 
# index3<-read.table("L1_feature_threeway_start.txt")[,1]
# length(index3)
# threeway_start<-colnames(design3)[index3][1:16]


#rbind(twoway_start,threeway_start)


###MLE of shared variables
glmdata<-as.data.frame(cbind(bank_y,design2[,-1]))

colnames(glmdata)[1]<-'y'
model1<-glm(formula(paste('y',paste(variables3[-1], collapse=" + "), sep=" ~ ")),data=glmdata,family=poisson)

model2<-glm(formula(paste('y',paste(variables2[-1], collapse=" + "), sep=" ~ ")),data=glmdata,family=poisson)

model3<-glm(formula(paste('y',paste(share_variables[-1], collapse=" + "), sep=" ~ ")),data=glmdata,family=poisson)

summary(model1)
summary(model2)
summary(model3)

model4 <- update(model1, . ~ . - contacttelephone:poutcomenonexistent )
summary(model4)
model44 <- update(model1, . ~ . - monthoct:poutcomenonexistent  )
summary(model44)

model5<- update(model2, . ~ . - monthmar-monthoct-educationuniversity.degree)
summary(model5)

model6 <- update(model3, . ~ . - monthmar-monthoct-educationuniversity.degree)
summary(model6)


summary(model1)




glmdata3<-as.data.frame(cbind(bank_y,design3[,-1]))
library('corrplot') #package corrplot
corrplot(cor(glmdata3[,index2]),  method="number",tl.pos="d") #plot matrix
corrplot(cor(glmdata3[,index3]),  method="number",tl.pos="d") #plot matrix

corrplot(cor(glmdata3[,intersect(index2,index3)]),  method="number",tl.pos="d") #plot matrix


with(model1, cbind(res.deviance = deviance, df = df.residual,
               p = pchisq(deviance, df.residual, lower.tail=FALSE))) 

exp(0.16+0.57)
exp(0.93)
exp(-0.47)
exp(-0.88)
pchisq(sum(residuals(model1, type="pearson")^2), 18048)


colnames(design1)
 
