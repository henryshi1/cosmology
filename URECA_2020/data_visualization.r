ls("package:datasets")
head(euro)
?euro

#Example 1: airquality
#New York air quality and climate measurements

airquality
?airquality()
head(airquality)
str(airquality)

airquality$Month<-as.factor(airquality$Month)
airquality$Day<-as.factor(airquality$Day)
str(airquality)

#First, name a histogram
hist(airquality$Temp)
library(ggplot2)

ggplot(data=airquality, aes(x=Temp)) + geom_histogram(fill="gray",color="black")
