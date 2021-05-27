#!/bin/bash
# estmapstat.sh *.gmap.out.gz 
# summarize est mapping rates from GMAP output files

gmapset=$*

#no# if [ "X" = "X$evigene" ]; then evigene=/bio/bio-grid/mb/evigene; fi
#no# if [ ! -d $evigene ]; then echo "ERR: need evigene=/path/to/evigene param"; exit -1; fi

noh=0
for goz in $gmapset; do {
 if [ $goz = "stdin" ]; then
   nn="stdin"
   CAT="cat -"
 else
   ena=`basename $goz .gz | sed 's/.gmap.*\.out.*//;'`
   nn=`echo $ena | sed 's/\..*//;'`
   CAT="gunzip -c $goz "
 fi 
 
 $CAT | egrep '^\>|^  Path |^Paths| Percent identity:| Coverage:' | env nohd=$noh nam=$nn perl -ne\
'if(/^>(\S+)/){ $id=$1; $keep=(/>cgb:P_/)?0:1; $nt++ if($keep); $ischi=0; }
elsif(/^Paths .(\d+)/) { $np=$1; $ischi=(/chimera/)?1:0; $nchi++ if($ischi); }
elsif(/ Path (\d+):/){ $p=$1; } elsif(/ Coverage:/  and ($p == 1 or ($ischi and $p == 2)) and $keep)
{ ($cov, $qlen)= m/\s+([\d\.]+) .query length: (\d+)/; $scov+=$cov; }
elsif(/ Percent / and ($p == 1) and $keep)
{ @v= m/: (\d+)\S* \((\d+) matches, (\d+) mismatches, (\d+) indels/; $n++;
for $i (0..3) { $sv[$i]+= $v[$i]; } $ib=$v[3]/($v[1]||1); $sv[4]+=$ib; }
BEGIN{ $nohd=$ENV{nohd}||0; $nam=$ENV{nam} ||"noname";
unless($nohd) { printf "%-16s\t","Group"; print join("\t",
  qw(Nin Nmap %mapped %cover %ident match  mismat indel), "indel/bp","chimer\n"); } }
END {  printf "%-16s\t$nt\t$n\t%.2f\t%.2f",$nam, (100*$n/$nt), $scov/$n;
$n||=1; $sv[4]= $n*$sv[3]/$sv[1]; $sv[5]= $nchi; @fm=("%.2f","%.0f","%.2f","%.2f","%.4f","%.3f");
for $i (0..5){ printf "\t$fm[$i]",$sv[$i]/$n;} print "\n"; }'  

  noh=1
} done

