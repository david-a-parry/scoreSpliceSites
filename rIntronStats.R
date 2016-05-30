#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("At least two arguments must be supplied (intron stats file and exon stats file).n", call.=FALSE)
} else if (length(args) == 2) {
  # default output file
  args[3] = "stats_out_dir"
}else if (length(args) > 3){
  stop("Too many arguments!.n", call.=FALSE)
}


dir.create(args[3])

types <- c("intron", "exon")
logcats  <- c("LENGTH", "REPEATS", "TOTAL_REPEAT_LENGTH", "MEAN_REPEAT_LENGTH")

for (i in 1:2){
  browser()
  type = types[i]
  data <-  read.table(args[i], header=T)
  u12 <- data[data$INTRON_TYPE == 'U12', ]
  u2 <- data[data$INTRON_TYPE == 'U2', ]

  sink(paste(args[3], paste(type, "Stats.txt", sep=""), sep="/"), split=T)  
  
  print("U12 %GC")
  print(summary(u12[, "PERCENT_GC"]))
  print("U2 %GC")
  print (summary(u2[, "PERCENT_GC"]) )
  print ( wilcox.test(u2[, "PERCENT_GC"], u12[, "PERCENT_GC"]) )
  png(paste(args[3], paste(type, "GCpercentBoxplot.png", sep=""), sep="/"), height=800, width=800 )
  boxplot(u2[, "PERCENT_GC"], u12[, "PERCENT_GC"], names = c("U2", "U12"), col=c( "brown1", "cyan" ))
  title(main = paste(type, "%GC vs Intron Type"), ylab = "%GC")
  for (l in logcats){
    print(paste("U12", l))
    print(summary(u12[, l]))
    print(paste("U2", l))
    print(summary(u2[, l]))
    # u2log <- log10(u2[, l])
    # u12log <- log10(u12[, l])
    # print ( t.test(u2log[is.finite(u2log)], u12log[is.finite(u12log)]) )
    print ( wilcox.test(u2[, l], u12[, l]) )
    png(paste(args[3], paste(type, l, "Boxplot.png", sep=""), sep="/"), height=800, width=800)
    boxplot(u2[, l], u12[, l], log='y', names = c("U2", "U12"), col=c( "brown1", "cyan" ))
    title(main = paste(type, l, "vs Intron Type"), ylab = paste(l))
    dev.off()
  }
  sink()
}
 

# intronStats <- read.table(args[1], header=T)
# exonStats <- read.table(args[2], header=T) 
# u12Introns <- intronStats[ intronStats$INTRON_TYPE == 'U12', ]
# u2Introns <- intronStats[ intronStats$INTRON_TYPE == 'U2', ]

# t.test(log10(u2Introns$LENGTH), log10(u12Introns$LENGTH))
# 
# png(paste(args[3], "intronLengthBoxplot.png", sep="/"))
# boxplot(u2Introns$LENGTH, u12Introns$LENGTH, log="y", names = c("U2", "U12"))
# title(main = "Intron Length vs Intron Type", ylab = "Intron Length (bp)")
# dev.off()
# 
# 
# t.test(u2Introns$PERCENT_GC, u12Introns$PERCENT_GC)
# 
# png(paste(args[3], "intronGCpercentBoxplot.png", sep="/"))
# boxplot(u2Introns$PERCENT_GC, u12Introns$PERCENT_GC, names = c("U2", "U12"))
# title(main = "Intron %GC vs Intron Type", ylab = "%GC")
# dev.off()
# 
# t.test(log10(u2Introns$REPEATS), log10(u12Introns$REPEATS))
# 
# png(paste(args[3], "intronRepeatNumberBoxplot.png", sep="/"))
# boxplot(u2Introns$REPEATS, u12Introns$REPEATS, names = c("U2", "U12"))
# title(main = "No. Repeats per Intron vs Intron Type", ylab = "No. Repeats")
# dev.off()
# 
# t.test(log10(u2Introns$TOTAL_REPEAT_LENGTH), log10(u12Introns$TOTAL_REPEAT_LENGTH))
# 
# png(paste(args[3], "intronTotalRepeatLengthBoxplot.png", sep="/"))
# boxplot(u2Introns$TOTAL_REPEAT_LENGTH, u12Introns$TOTAL_REPEAT_LENGTH, names = c("U2", "U12"))
# title(main = "Total Length of Repeats per Intron vs Intron Type", ylab = "Length (bp)")
# dev.off()
# 
# t.test(log10(u2Introns$MEAN_REPEAT_LENGTH), log10(u12Introns$MEAN_REPEAT_LENGTH))
# 
# png(paste(args[3], "intronMeanRepeatLengthBoxplot.png", sep="/"))
# boxplot(u2Introns$MEAN_REPEAT_LENGTH, u12Introns$MEAN_REPEAT_LENGTH, names = c("U2", "U12"))
# title(main = "Mean Length of Repeats per Intron vs Intron Type", ylab = "Length (bp)")
# dev.off()

