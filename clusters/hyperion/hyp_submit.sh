#!/bin/bash
# Submit the Hyperion chain of campaign dots: hyp_submit.sh <id>... (records in qdw_campaign/records, cp2k/chain.sh
# written). Shape: 240 ranks below 1,600 atoms, 360 from 1,600 (Ivan's HgTe 1,299-atom jobs: 240 ranks over any free
# cores); memory per rank from the pilot (1,300 atoms GPW on 112 ranks: ~500 GB in all, largest rank 4.6 GB), x 1.5,
# ~ n^2, at least 4G: mem_per_cpu = max(4, ceil(1.5 * 520 (n/1300)^2 / ranks)) GB. Written to cp2k/.hyp_shape (every
# later piece uses it), then the first piece is submitted. Prints "<job> <id> <shape>".
C=/scratch/iinfante/qdw_campaign
for id in "$@"; do
    d=$C/records/$id/cp2k
    [ -s $d/chain.sh ] || { echo "$id: no cp2k/chain.sh"; continue; }
    [ -f $d/.done ] && { echo "$id: done"; continue; }
    na=$(head -1 $C/records/$id/start.xyz)
    shape=$(python3 -c "
import math; n=$na; r=240 if n < 1600 else 360
print(f'-n {r} --mem-per-cpu={max(4, math.ceil(1.5 * 520 * (n / 1300) ** 2 / r))}G')")
    echo "$shape" > $d/.hyp_shape
    j=$(cd $d && sbatch --parsable -J hyp-$id $shape $C/bin/hyp_chain.sh) && echo "$j $id $shape"
done
