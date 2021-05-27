args = commandArgs(trailingOnly=TRUE)
sample_name=args[1]

library(ggthemes)
library(ggplot2)

data=read.csv(paste(sample_name,"_custom_uniprot_hits.txt",sep=""),header=F)

nlim=round((head(data$V1,n = 1)+450),digits = -2)
p1<-ggplot(data=data, aes(x=reorder(data$V2,data$V1), y=data$V1))+
  geom_bar(stat="identity", fill="dark blue", width=.5)+
  coord_flip()+labs(x="UniProt Species",y="Number of Hits")+
  geom_text(aes(label=data$V1), position=position_dodge(width=0.3), vjust=0.25, hjust=-.10)+
  theme(axis.text=element_text(size=12))+ylim(0,nlim)+theme(axis.text.x=element_text(size=12,angle=0))+
  theme(axis.title=element_text(size=15,face="bold"))+ggtitle(paste(sample_name,"UniProt hits",sep=" "))+
  theme(plot.title = element_text(family="sans", colour = "black", size = rel(1.5)*1, face = "bold"))


# not working in docker
#ggsave(filename = paste(sample_name,"_custom_uniprot_hits.svg",sep=""),width = 15 ,height = 7)
#ggsave(filename = paste(sample_name,"_custom_uniprot_hits.pdf",sep=""),width = 15 ,height = 7)
pdf(paste(sample_name,"_custom_uniprot_hits.pdf",sep=""),width = 15 ,height = 7)
print(p1)
dev.off()
svg(paste(sample_name,"_custom_uniprot_hits.svg",sep=""),width = 15 ,height = 7)
print(p1)
dev.off()
