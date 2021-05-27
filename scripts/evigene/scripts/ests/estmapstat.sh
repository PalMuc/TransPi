#!/bin/tcsh
# estmapstat.sh *.gmap.out.gz 

set gmapset="$*"
set evigene=/bio/bio-grid/mb/evigene
set noh=0

foreach goz ( $gmapset )
 if ( $goz == "stdin" ) then
   set nn="stdin"
   set CAT="cat -"
# debug this
   # echo "DEBUG..."
   # $CAT | egrep '^\>|^  Path | Percent identity:| Coverage:' | head -20
   # echo "DONE..."
   # exit;

 else
 set ena=`basename $goz .gz | sed 's/.gmap.*\.out.*//;'`
 set nn=`echo $ena | sed 's/.s[1-9]r.//; s/.mars11//; s/.kfish2a//; s/Assembly/a/; s/reads/er/;'`
 set CAT="gunzip -c $goz "
 endif
 
## add /^Paths/ for "*** Possible chimera" counts?
## Paths (2): *** Possible chimera with breakpoint at 175..178

 $CAT | egrep '^\>|^  Path |^Paths| Percent identity:| Coverage:' | env nohd=$noh nam=$nn perl -ne\
'if(/^>(\S+)/){ $id=$1; $keep=(/>cgb:P_/)?0:1; $nt++ if($keep); $ischi=0; }\
elsif(/^Paths .(\d+)/) { $np=$1; $ischi=(/chimera/)?1:0; $nchi++ if($ischi); } \
elsif(/ Path (\d+):/){ $p=$1; } elsif(/ Coverage:/  and ($p == 1 or ($ischi and $p == 2)) and $keep) \
{ ($cov, $qlen)= m/\s+([\d\.]+) .query length: (\d+)/; $scov+=$cov; } \
elsif(/ Percent / and ($p == 1) and $keep) \
{ @v= m/: (\d+)\S* \((\d+) matches, (\d+) mismatches, (\d+) indels/; $n++; \
for $i (0..3) { $sv[$i]+= $v[$i]; } $ib=$v[3]/($v[1]||1); $sv[4]+=$ib; } \
BEGIN{ $nohd=$ENV{nohd}||0; $nam=$ENV{nam} ||"noname"; \
unless($nohd) { printf "%-16s\t","Group"; print \
join("\t",qw(Nin Nmap %mapped %cover %ident match  mismat indel), "indel/bp","chimer\n"); } }\
END {  printf "%-16s\t$nt\t$n\t%.2f\t%.2f",$nam, (100*$n/$nt), $scov/$n; \
$n||=1; $sv[4]= $n*$sv[3]/$sv[1]; $sv[5]= $nchi; @fm=("%.2f","%.0f","%.2f","%.2f","%.4f","%.3f");\
for $i (0..5){ printf "\t$fm[$i]",$sv[$i]/$n;} print "\n"; }'  

  set noh=1
end

