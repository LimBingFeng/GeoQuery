#loading library
library(GEOquery)
library(dplyr)
library(TmCalculator)
library(ape)
library(ggpubr)

#download files from NCBI GEO
crc = getGEO("GSE92921")

#data cleaning
genes = crc$GSE92921_series_matrix.txt.gz@featureData@data
genesExprs = crc$GSE92921_series_matrix.txt.gz@assayData$exprs

msh6 = genes[8997,]
msh6 = msh6 %>% select(c("ID","Gene Title","Gene Symbol","RefSeq Transcript ID"))

msh2 = genes[18835,]
msh2 = msh2 %>% select(c("ID","Gene Title","Gene Symbol","RefSeq Transcript ID"))

msh6Exprs = genesExprs[8997,]
msh2Exprs = genesExprs[18835,]

#Checking for outlier
boxplot(msh6Exprs, main ="Boxplot of MSH6", col = "light blue") ## 3 outliers
boxplot(msh2Exprs, main ="Boxplot of MSH2", col = "orange") ## 1 outlier

#plotting graph
plot(msh6Exprs, main = "MSH6 Expression Values", type = 'l', xlab = "Samples", ylab = "MSH6 Expression Value", col = "light blue")
plot(msh2Exprs, main = "MSH2 Expression Values", type = 'l', xlab = "Samples", ylab = "MSH2 Expression Value", col = "orange")

#checking normality (msh6 not normally distributed, msh2 normally distributed)
hist(msh6Exprs, probability = T, main = "Normality of MSH6",xlab = "skewed data", ylim = c(0,0.25),col = "light blue")
lines(density(msh6Exprs),col = "red")

hist(msh2Exprs, probability = T, main = "Normality of MSH2",xlab="Approximately normally distributed data",ylim = c(0,0.001),col = "orange")
lines(density(msh2Exprs),col = "red")


qqnorm(msh6Exprs, main = "QQ plot of MSH6", pch = 19, col = "light blue")
qqline(msh6Exprs, col = "red")

qqnorm(msh2Exprs, main = "QQ plot of MSH2", pch = 19, col = "orange")
qqline(msh2Exprs, col = "red")


shapiro.test(msh6Exprs) ## low p-value
shapiro.test(msh2Exprs) ## high p-value


#calculate exprs mean and median
mean(msh6Exprs)
median(msh6Exprs)
summary(msh6Exprs)

mean(msh2Exprs)
median(msh2Exprs)
summary(msh2Exprs)

#download Genes Fasta file
msh6Fasta = read.GenBank(msh6$`RefSeq Transcript ID`, as.character = T)
msh2Fasta = read.GenBank(msh2$`RefSeq Transcript ID`, as.character = T)
msh6Fasta = msh6Fasta[[1]]
msh2Fasta = msh2Fasta[[1]]

#GC content
GC(msh6Fasta)
GC(msh2Fasta)
#comparing GC content
barplot(c(GC(msh6Fasta),GC(msh2Fasta)), main = "GC content", names.arg = c("MSH6","MSH2"), 
        col = c("light blue","orange"), ylab = "GC Content", xlab = "Genes", ylim=c(0,50))


#create pie chart
msh6Content = table(msh6Fasta)
msh2Content = table(msh2Fasta)
#MSH6 pie chart
percentage = round(msh6Content/sum(msh6Content)*100)
percentage = paste(percentage, "%", sep="")
lbls = paste(c("A","C","G","T"),percentage,sep = " ")
pie(msh6Content,labels = lbls, main = "Base Compostion of MSH6 Gene")
#MSH2 pie chart
percentage = round(msh2Content/sum(msh2Content)*100)
percentage = paste(percentage, "%", sep="")
lbls = paste(c("A","C","G","T"),percentage,sep = " ")
pie(msh2Content,labels = lbls, main = "Base Compostion of MSH2 Gene")

