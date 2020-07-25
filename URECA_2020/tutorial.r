#Install package
install.packages("ggplot2")

# Activate package
library(ggplot2)

# DATA TYPE

###Vector###
a<-c(1,4,67,-7,8,67)
b<-c("I", "have", "osteoporosis")

c<-c(TRUE,TRUE,FALSE)

class(a)
class(b)
class(c)

# Convert datatypes
b<-c("I", "have", "osteoporosis", "crippling", "osteoporosis")
b1<-as.factor(b)
b1
# 1 assigned to crippling, 2 to have, 3 to I, 4 to osteoporosis,
# and 4 to osteoporosis again

# call the second element of vector
b[2]
b[c(2,5)] # call the second and fifth elements, concatenated
c(a,b) # convert to character (concatenating numeric and character)

# Matrices
1:20 # generates a set of numbers (numeric) from 1 to 20
x<-matrix(1:20,nrow=5,ncol=4) # fill a matrix with said numbers
x
?matrix
y<-matrix(1:20,nrow=5,ncol=4, byrow=TRUE) # fill a matrix with said numbers
y

y[row,column]
y[1,2] #show me 1st row, 2nd column
y[1,] # show me entire row
y[,1] # show me entire column
y[1:3,1:2]  # show me 1st-3rd rows and 1st-2nd columns

###Dataframe###
# like a matrix but where elements are not necessarily numeric
d<-c(1,2,3,4)
e<-c("red","white","red",NA)
f<-c(TRUE,FALSE,FALSE,TRUE)
df<-data.frame(d,e,f)
df

colnames(df)<-c("ID","color","passed")

df[3:5,] # displays rows 3-5 of dataframe
df[c("ID","color")] # displays the columns with names "ID" and "color"
df$ID # displays every element of "ID" as a numeric

###IMPORT DATASET###
?read.csv    # the question mark is asking for documentation on a command
df<-read.csv("Davis.csv")
df
