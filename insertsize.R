#setwd('/media/ninadmw/sgate4/test/')
stats = read.table('Insert_count.txt', sep = '', header = T)
names(stats) <- c("count", "insert")
x <- colSums(stats)[1]

stats2 = stats[stats$insert < 601,]
y <- colSums(stats2)[1]

y/x*100
frac <- round(y/x*100, digits =2)
print(paste("covers ",frac," % of reads", sep=''))

library(ggplot2)
pdf("insertsize.pdf");
print(ggplot(stats2, aes (x = insert)) +
              geom_line(aes(y=count/1000))+
              labs(x="Insert size in bp", y='Read count in thousands', title= paste("Insert size distribution for ",frac, " % of reads", sep =''))+
              theme_bw() + 
              theme(                              
                      axis.text.x = element_text(face="bold", color="black", size=14),
                      axis.text.y = element_text(face="bold", color="black", size=14),
                      axis.title.x = element_text(face="bold", color="black", size=14),
                      axis.title.y = element_text(face="bold", color="black", size=14),
                      plot.title = element_text(face="bold", color = "black", size=16, hjust=0.5)
              ))
 dev.off()