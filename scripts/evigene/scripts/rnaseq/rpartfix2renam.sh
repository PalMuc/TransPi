#!/bin/bash

mv aphidpe_SRR075803.fa2 aphidpe_SRR075803.fa2.tmp
mv aphidpe_SRR075803.fa2.pairs aphidpe_SRR075803.fa2
cat aphidpe_SRR075803.fa2.unpair >> aphidpe_SRR075803.fa1

mv aphidpe_SRR075802.fa2 aphidpe_SRR075802.fa2.tmp
mv aphidpe_SRR075802.fa2.pairs aphidpe_SRR075802.fa2
cat aphidpe_SRR075802.fa2.unpair >> aphidpe_SRR075802.fa1

# for nam in (aphidpe_SRR071347)
mv aphidpe_SRR071347.bam.fa2 aphidpe_SRR071347.bam.fa2.tmp
mv aphidpe_SRR071347.bam.fa2.pairs aphidpe_SRR071347.fa2
mv aphidpe_SRR071347.bam.fa2.unpair aphidpe_SRR071347.fa1

