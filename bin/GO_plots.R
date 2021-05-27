args = commandArgs(trailingOnly=TRUE)
sample_name=args[1]

library(ggthemes)
library(ggplot2)

dataCC=read.delim("GO_cellular.txt", header = F, sep = "\t")
dataMF=read.delim("GO_molecular.txt", header = F, sep = "\t")
dataBP=read.delim("GO_biological.txt", header = F, sep = "\t")

#CC
nlim=round((head(dataCC$V1,n = 1)+150),digits = -2)
p1<-ggplot(data=dataCC, aes(x=reorder(dataCC$V2,dataCC$V1), y=dataCC$V1))+
  geom_bar(stat="identity", fill="green", width=.5)+
  coord_flip()+labs(x="Classification",y="Number of Sequences")+
  geom_text(aes(label=dataCC$V1), position=position_dodge(width=0.7), vjust=-0.0005, hjust=-.15)+
  theme(axis.text=element_text(size=10))+ylim(0,nlim)+theme(text = element_text(size = 15))+
  theme(axis.text.x=element_text(size=12,angle=0))+theme(axis.title=element_text(size=15,face="bold"))+
  ggtitle(paste(sample_name,"Cellular Componenet GOs",sep=" "))+
  theme(plot.title = element_text(family="sans", colour = "black", size = rel(1.1)*1, face = "bold"))

#ggsave(filename = paste(sample_name,"_Cellular_Component.svg",sep=""),width = 15 ,height = 7)
#ggsave(filename = paste(sample_name,"_Cellular_Component.pdf",sep=""),width = 15 ,height = 7)
pdf(paste(sample_name,"_Cellular_Component.pdf",sep=""),width = 15 ,height = 7)
print(p1)
dev.off()
svg(paste(sample_name,"_Cellular_Component.svg",sep=""),width = 15 ,height = 7)
print(p1)
dev.off()

#MF
nlim=round((head(dataMF$V1,n = 1)+150),digits = -2)
p2 <-ggplot(data=dataMF, aes(x=reorder(dataMF$V2,dataMF$V1), y=dataMF$V1))+
  geom_bar(stat="identity", fill="blue", width=.5)+
  coord_flip()+labs(x="Classification",y="Number of Sequences")+
  geom_text(aes(label=dataMF$V1), position=position_dodge(width=0.7), vjust=-0.0005, hjust=-.15)+
  theme(axis.text=element_text(size=10))+ylim(0,nlim)+theme(text = element_text(size = 15))+
  theme(axis.text.x=element_text(size=12,angle=0))+theme(axis.title=element_text(size=15,face="bold"))+
  ggtitle(paste(sample_name,"Molecular Function GOs",sep=" "))+
  theme(plot.title = element_text(family="sans", colour = "black", size = rel(1.1)*1, face = "bold"))

#ggsave(filename = paste(sample_name,"_Molecular_Function.svg",sep=""),width = 15 ,height = 7)
#ggsave(filename = paste(sample_name,"_Molecular_Function.pdf",sep=""),width = 15 ,height = 7)
pdf(paste(sample_name,"_Molecular_Function.pdf",sep=""),width = 15 ,height = 7)
print(p2)
dev.off()
svg(paste(sample_name,"_Molecular_Function.svg",sep=""),width = 15 ,height = 7)
print(p2)
dev.off()

#BP
nlim=round((head(dataBP$V1,n = 1)+150),digits = -2)
p3<-ggplot(data=dataBP, aes(x=reorder(dataBP$V2,dataBP$V1), y=dataBP$V1))+
  geom_bar(stat="identity", fill="red", width=.5)+
  coord_flip()+labs(x="Classification",y="Number of Sequences")+
  geom_text(aes(label=dataBP$V1), position=position_dodge(width=0.7), vjust=-0.0005, hjust=-.15)+
  theme(axis.text=element_text(size=10))+ylim(0,nlim)+theme(text = element_text(size = 15))+
  theme(axis.text.x=element_text(size=12,angle=0))+theme(axis.title=element_text(size=15,face="bold"))+
  ggtitle(paste(sample_name,"Biological Processes GOs",sep=" "))+
  theme(plot.title = element_text(family="sans", colour = "black", size = rel(1.1)*1, face = "bold"))

#ggsave(filename = paste(sample_name,"_Biological_Processes.svg",sep=""),width = 15 ,height = 7)
#ggsave(filename = paste(sample_name,"_Biological_Processes.pdf",sep=""),width = 15 ,height = 7)
pdf(paste(sample_name,"_Biological_Processes.pdf",sep=""),width = 15 ,height = 7)
print(p3)
dev.off()
svg(paste(sample_name,"_Biological_Processes.svg",sep=""),width = 15 ,height = 7)
print(p3)
dev.off()
