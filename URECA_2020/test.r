# Install package
#install.packages("ggplot2")

# Activate package
library(ggplot2)


# DATA TYPE

##Vector##
a<-c(1,4,7,-7,8,67)
b<-c("black", "lives", "matter")
c<-c(TRUE,TRUE,TRUE)

class(a)
class(b)
class(c)

b<-c("black","lives","matter","really","matter")
b1<-as.factor(b)

b[2]
b[c(2,5)]
c(a,b)

#Matrices
1:20
x<-matrix(1:20,nrow=5,ncol=4)
x
?matrix
y<-matrix(1:20,nrow=5,ncol=4,borrow=TRUE)
y

y[row,column]
y[1,]
y[,1]
y[1:3,1:2]

###Dataframe###
d<-c(1,2,3,4)
e<-c("red","white","red",NA)
f<-c(TRUE,FALSE,FALSE,TRUE)
df<-data.frame(d,e,f)
df
colnames(df)<-c("ID","color","passed")
df[3:5,]
df[c("ID","color")]
df$ID

###Import DATASET###
?read.csv
df<-read.csv("/home/seriouscomedian/code/GitHub/cosmology/URECA_2020/Davis.csv",header=TRUE,sep=",",na.strings="NA")
head(df)

df<-read.csv("/home/seriouscomedian/code/GitHub/cosmology/URECA_2020/Davis.csv",header=TRUE,sep="\t",na.strings="NA")
head(df)

setwd("/home/seriouscomedian/code/GitHub/cosmology/URECA_2020/")
df<-read.csv("Davis.txt",header=TRUE,sep="\t",na.strings="NA")
head(df)


colnames(df)[1]<-"RowID"
head(df)
View(df)
dim(df)
summary(df)
str(df)
df$sex
length(df$sex)
table(df$sex)
sort(df$weight, decreasing=TRUE)

#Quick visualisation
plot(df$weight~df$RowID)
plot(df$height~df$RowID)
boxplot(df$weight~df$sex)
plot(df$weight,df$repwt)
abline(0,1,col="red")
cor.test(df$weight,df$repwt)


#Subset a DATASET
girl<-df[df$sex %in% "F",]
head(girl)
dim(girl)
girl<-df[df$sex == "F",]
girl<-subset(df,df$sex %in% "F")

boy<-subset(df,df$sex %in% "M")
boy

?rbind
rbind()
cbind()

dfbis<-rbind(girl,boy)
dim(dfbis)
