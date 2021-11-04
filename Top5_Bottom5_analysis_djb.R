#Analysis of GS Pilot II data to determine necessary sample size

#setwd to source file location
setwd("/Users/danbarshis/dansstuff/Research/2021-2024_GlobalSearch/2021-08_Pilotv2Reproducability")

library(data.table)
library(Hmisc)
df <- read.csv("Poc_ED50s.csv")
df$meanED50<-rowMeans(df[,2:3]) #compute average ED50 between run 1 and run 2

#populate empty data frame for results from loop
ResOutput<-data.frame("numsamples"=numeric(), "numoverlapping"=numeric(), "pvalue"=numeric())

#run loop for every number of samples between 10 and 40
for (i in 10:nrow(df)) {
#run 1000 times for each sample number
  for (j in 1:1000){
#randomly sample i rows from data frame
    x<-df[sample(nrow(df),i),]
#compute ranks of 1st run and 2nd run
    x$Rank1<-rank(-x$ED50_R1)
    x$Rank2<-rank(-x$ED50_R2)
#compute top5 and bottom 5 from 1st run ranks
    dt1 <- data.table(x, key="Rank1")
    R1_top5<-head(dt1, n=5)
    R1_bot5<-tail(dt1, n=5)
#output number of overlapping samples in run 2 (i.e., the number of top ranks from run 1 that are in the bottom ranks of run 2)
    numoverlaptopbottom<-sum(R1_top5$Rank2>R1_bot5$Rank2)
    numoverlapbottomtop<-sum(R1_bot5$Rank2<R1_top5$Rank2)
#test whether average ED50s of top 5 are higher than average ED50s of bottom 5
    TestPvalue<-t.test(R1_top5$meanED50,R1_bot5$meanED50,alternative="g")$p.value
    ResOutput=rbind(ResOutput,data.frame("numsamples"=i, "numoverlapping"=numoverlaptopbottom+numoverlapbottomtop, "pvalue"=TestPvalue))
}
  }
sum(ResOutput$pvalue>0.05)
OutputSummary<-data.frame(
"numsamples"=10:40,
"Sigpvalue"=aggregate(pvalue~numsamples, function(x){sum(x<0.05)}, data=ResOutput)[[2]],
"meanoverlap"=aggregate(numoverlapping~numsamples, mean, data=ResOutput)[[2]],
"minoverlap"=aggregate(numoverlapping~numsamples, min, data=ResOutput)[[2]],
"maxoverlap"=aggregate(numoverlapping~numsamples, max, data=ResOutput)[[2]],
"medianoverlap"=aggregate(numoverlapping~numsamples, median, data=ResOutput)[[2]],
"stddevoverlap"=aggregate(numoverlapping~numsamples, sd, data=ResOutput)[[2]]
)
pdf(file="Pverrucosa_reproducibility.pdf",7,14)
par(mfrow=c(2,1))
plot(OutputSummary$numsamples,OutputSummary$Sigpvalue, ylab="Number of significant comparisons", xlab="Number of samples", main="Num Sig comparisons top5 vs. bottom5 ED50s")

errbar(OutputSummary$numsamples,OutputSummary$meanoverlap, OutputSummary$meanoverlap+OutputSummary$stddevoverlap, OutputSummary$meanoverlap-OutputSummary$stddevoverlap,ylab="Mean overlap +- 1 stddev", xlab="Number of samples")
title(main="Mean overlap top5 vs. bottom5 ED50s")
dev.off()

plot(OutputSummary$numsamples,OutputSummary$meanoverlap, type="l")
plot(OutputSummary$numsamples,OutputSummary$maxoverlap, type="l")
plot(OutputSummary$numsamples,OutputSummary$medianoverlap, type="l")
