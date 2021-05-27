#!/bin/bash
# run_evgomcltab.sh : should be part of run_evgomcl

date=20171110
idprefix=CIF7b_G
# idprefix=ARPCIF_G # default: "ARP9_G"

#myspp for brief table
myspecies=dpx17evg

speciesmap=bemtab=Whitefly,dapmaevg14=Daphnia_magna,dpx17evg=Daphnia_pulex,drosmel16nc=Fruitfly,tribcas16nc=Beetle,guppync1=Guppy,zebrafish3=Zebrafish

export clade=CrustInsectFish

## mcl : used? fix path
export MCL=/bio/bio-grid/mb/mcl9
export evigene=/bio/bio-grid/mb/evigene

env myspecies=$myspecies xml=0 date=$date clade=$clade speciesmap="$speciesmap" \
 $evigene/scripts/omcl/orthomcl_tabulate.pl -debug -idprefix=$idprefix -omclpath ./ -namepath ../names

#______________________
# arp8daph.spplist
#   13901 bemtab
#   29127 dapmaevg14
#   33040 dpx17evg
#   13916 drosmel16nc
#   22914 guppync1
#   12861 tribcas16nc
#   26247 zebrafish3
#..........
# env sppmap=FISH11 myspecies=Fundulus orthoutisokay=1 xml=0 date=20131225 clade=Fish \
#  goodname='mayzebr|human|kfish2|platyfish' poorname='stickleback|medaka|tetraodon' \
#  sppindex="0,2,3,4,5,7,8,9,10" tmin=7 $evigene/scripts/omcl/orthomcl_tabulate.pl \
#   -debug -idprefix="FISH${pt}_G" -gtag="FISH$pt" -omclpath ./ -namepath ../names/

