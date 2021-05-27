## rxdmag5x_edger_pair.R
## R:edgeR stats -----------------------------------------------
## cluster stub program, PBS shell caller sets
#. p1 = 'CO'
#. p2 = 'BN'
#. indat = 'rsamdmag5x_hsx.tocounts' 

library(edgeR)

GSNAM <- "dmag5xau13evc_hssam"   # caller should set, other opts?
MINCOUNT <- 50
qval <- FDR <- 0.05;

tab<- read.delim(indat,header=T,comment.char = "#") 

rownames(tab) <- tab[,1]; tab<- tab[,-1]
colnames(tab) <- sub("HS_CO_X_r4","HS_CO_X_r3",colnames(tab))  #? still need
tcol<- colnames(tab)
src<- gsub("_.*","",tcol)  # all HS or ND
treat<- gsub("^.._(\\w\\w)_.*","\\1",tcol) # BN .. TR
clone<- gsub(".*_(\\w+)_.*","\\1",tcol) # all X
repl <- gsub(".*_r","",tcol)  ## all 1,2,3 but for change HS_CO_X_r4 to HS_CO_X_r3
tgroup<- gsub("_r[0-9]*$","",tcol)

#............................................................
# Pairwise 2-group Treatment or Clone, using exactTest(dct)
# Treatment only, collapse clones.
# treat: CO x other only

targets <- data.frame(cbind(treat))
rownames(targets) <- tcol
design <- model.matrix( ~ treat, data = targets)  #? need for pairwise?

# all alts
sum( rx <- (apply(tab,1,max) >= MINCOUNT) )  
dc <- DGEList(counts=round(tab[rx,]), group=treat)
dc <- calcNormFactors(dc)
dct <- estimateCommonDisp(dc, design) 

#OFF# p1 <- "CO"; p2 <- "BN" # caller sets
#OFF# for ( p2 in unique(treat) ) { # NO loop, caller does
#OFF#   if( p1 == p2 ) next

tpnam= paste(p2,p1,sep="x")  

dct.xt <- exactTest(dct,pair=c(p1,p2)) 
nsig <- sum(dct.xt$table$p.value < FDR) 
dctOut<- topTags(dct.xt,n=nsig)$table
dctOut<- dctOut[ dctOut[,"FDR"] <= FDR, ];
dct.genesig <- rownames(dctOut)
outf<- paste(GSNAM,tpnam,"edger.txt",sep=".")
write.table( dctOut,outf,quote=F,sep="\t",row.names=T)
system(paste("perl -pi -e \'s/(\\.\\d\\d\\d\\d)\\d+/$1/g;\'",outf)) # skip?

pdf(paste(GSNAM,tpnam,"fc.pdf",sep="."))
plotSmear(dct, pair=c(p1,p2), cex=0.26, de.tags = dct.genesig, 
    main = paste(GSNAM,"FC plot for",tpnam) )
abline(h = c(-2, 2), col = "dodgerblue", lwd = 2)
dev.off()
 
#OFF# }

#............................................................
#. ##! /bin/bash
#. ### qsub -q batch edgerbatch.sh
#. #PBS -N derstat
#. #PBS -l vmem=48gb,nodes=1:ppn=8,walltime=20:55:00
#. #PBS -o derstat.$$.out
#. #PBS -e derstat.$$.err
#. #PBS -V
#. 
#. module add R
#. module add bioconductor
#. 
#. ncpu=8
#. workd=$HOME/scratch/chrs/daphmag/rnas/XXXX
#. indat=$workd/XXX.tocounts
#. rtempl=$workd/rxdmag5x_edger_pair.R
#. odir=de2xt
#. 
#. cd $workd/
#. mkdir $odir
#. cd $odir
#. 
#. for treat in BN BX CA CR FI PA TR; do {
#.   R_script="`basename $rtempl .R`.$treat.R"
#.   echo "p1 = 'CO'" > $R_script
#.   echo "p2 = '$treat'" >> $R_script
#.   echo "indat = '$indat'" >>  $R_script 
#.   cat $rtempl >> $R_script
#.   R --vanilla < $R_script 2>log.$treat &
#. } done
#. 
#. wait
#. 
###.............

