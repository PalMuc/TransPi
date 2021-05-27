# gg-velcall.sh

workd=$HOME/scratch/chrs/aphid2
cd $workd/rnas/$part_dir
$workd/rnas/gg-velrun.sh > log.vel$$ 2>&1 &
# wait
