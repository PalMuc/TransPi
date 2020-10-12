args = commandArgs(trailingOnly=TRUE)
sample_name=args[1]
comp_table=args[2]
transpi_table=args[3]

library(plotly)
library(reshape2)

# comparison table
csv=read.csv(comp_table, header=TRUE, sep="\t")

csv <- data.frame(lapply(csv, function(x) {gsub("Complete", "3", x)}))
csv <- data.frame(lapply(csv, function(x) {gsub("Duplicated", "2", x)}))
csv <- data.frame(lapply(csv, function(x) {gsub("Fragmented", "1", x)}))
csv <- data.frame(lapply(csv, function(x) {gsub("Missing", "0", x)}))
csv
c=melt(csv,id.vars = 'Busco.ID')
dec=c(0,.25,.25,.50,.50,.75,.75,1)
my_colors <- c("#081D58","#081D58", "#2280B8","#2280B8", "#99D6B9", "#99D6B9","#f8f9fc","#f8f9fc")
colz <- setNames(data.frame(dec, my_colors), NULL)
fig <- plot_ly(c,x=~variable, y=~Busco.ID, z=~value, colorscale=colz, reversescale=T, type = "heatmap",
              colorbar=list(tickmode='array', tickvals=c(.35,1.1,1.87,2.60), thickness=30,
              ticktext= c("Missing","Fragmented","Duplicated","Complete"), len=0.4))
fig <- fig %>% layout(xaxis=list(title="", showline = TRUE, mirror = TRUE),
              yaxis=list(title="BUSCO ID", tickmode="auto", nticks=length(csv$Busco.ID),
              tickfont=list(size=8), showline = TRUE, mirror = TRUE))

orca(fig, paste(sample_name,"_all_missing_BUSCO.png",sep=""))
orca(fig, paste(sample_name,"_all_missing_BUSCO.pdf",sep=""))

# TransPi table
csv=read.csv(transpi_table, header=TRUE, sep="\t")

csv <- data.frame(lapply(csv, function(x) {gsub("Complete", "3", x)}))
csv <- data.frame(lapply(csv, function(x) {gsub("Duplicated", "2", x)}))
csv <- data.frame(lapply(csv, function(x) {gsub("Fragmented", "1", x)}))
csv <- data.frame(lapply(csv, function(x) {gsub("Missing", "0", x)}))
csv
c=melt(csv,id.vars = 'Busco.ID')
dec=c(0,.25,.25,.50,.50,.75,.75,1)
my_colors <- c("#081D58","#081D58", "#2280B8","#2280B8", "#99D6B9", "#99D6B9","#f8f9fc","#f8f9fc")
colz <- setNames(data.frame(dec, my_colors), NULL)
fig <- plot_ly(c,x=~variable, y=~Busco.ID, z=~value, colorscale=colz, reversescale=T, type = "heatmap",
              colorbar=list(tickmode='array', tickvals=c(.35,1.1,1.87,2.60), thickness=30,
              ticktext= c("Missing","Fragmented","Duplicated","Complete"), len=0.4))
fig <- fig %>% layout(xaxis=list(title="", showline = TRUE, mirror = TRUE),
              yaxis=list(title="BUSCO ID", tickmode="auto", nticks=length(csv$Busco.ID),
              tickfont=list(size=8), showline = TRUE, mirror = TRUE))

orca(fig, paste(sample_name,"_TransPi_missing_BUSCO.png",sep=""))
orca(fig, paste(sample_name,"_TransPi_missing_BUSCO.pdf",sep=""))
