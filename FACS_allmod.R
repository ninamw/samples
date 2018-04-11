rm(list = ls())
Args=commandArgs(TRUE)
arg.data = Args[1]
arg.title = Args[2]
arg.nm = Args[3]

if (is.na(arg.data)== TRUE){arg.data <- "data.txt"}
if (is.na(arg.title)== TRUE){arg.title <- "96 Hr phenotype"}
if (is.na(arg.nm)== TRUE){arg.nm <- "apoptosis_all.pdf"}

plotFACS <- function(data, title, nm)
  
{
  library(lattice)
  sort = read.table(data, sep ='\t', header =T)
  pdf(nm, width=4.33, height=3.83)
  col.symbol=c(rgb(1,0,0,.6), rgb(0,0,1,.6))
  print(
    xyplot(Apoptotic~Cells, groups = Treatment,  data = sort,  pch=16, aspect =1, auto.key = TRUE, alpha = 0.5, ylab = '% apoptotic', xlab ="", scales = list(x=list(rot=45), alternating = FALSE, tck = c(1,0)),main='')
    )
  dev.off()
}

plotFACS(arg.data, arg.title, arg.nm)






 
