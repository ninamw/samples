#rm(list = ls())
#setwd('/media/ninadmw/sgate4/test/')
stats = read.table('project.stats', sep = '', header = T)

date <- Sys.Date()
attach(stats)
pdf(paste("TotalReads_",date,".pdf", sep=''));
print(plot(name, totalreads/1000000, main="Total Reads", 
            ylab="Total Reads in Millions ", pch=19))
 dev.off()


 pdf(paste("MappingRate_",date,".pdf", sep=''));
 print(plot(name, mappingrate, main="Mapping Rate", 
             ylab="% Mapping Rate ", pch=19,ylim=c(0, 100)))
 dev.off()
 
 pdf(paste("DupRate_",date,".pdf", sep=''));
 print(plot(name, duprate, main="Duplicate Reads", 
             ylab="% Duplicate Reads ", pch=19))
 dev.off()
 
 pdf(paste("MitoRate_",date,".pdf", sep=''));
 print(plot(name, mitorate, main="Mitochondrial Reads", 
             ylab="% mitochondrial Reads ", pch=19))
 dev.off()
 
 pdf(paste("FilteredReads_",date,".pdf", sep=''));
 print(plot(name, filteredreads/1000000, main="Filtered Reads", 
            ylab="Filtered Reads in Millions ", pch=19))
 dev.off()
 
 pdf(paste("FRiP_",date,".pdf", sep=''));
 print(plot(name, frip, main="Fraction of Reads in peaks", 
            ylab="FRiP score", pch=19))
 dev.off()
 

  
