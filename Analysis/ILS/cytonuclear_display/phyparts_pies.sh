#!/bin/bash

date;pwd

display=$(shuf -i 100-200 -n 1)
export DISPLAY=:${display}
Xvfb :${display} -screen 0 1024x768x16 > /dev/null 2>&1 &
echo "export DISPLAY=:${display}" > ~/.xvfb

#module load ete3
#python3 ./phypartspiecharts.py ../data/Prunus_speices_tree/Prunus_astral_species.br.rt.tre Prunus_phyparts 170 --svg_name Prunus_PhypartsPiecharts2.svg --to_csv
python  /data/xiongtao/scripts/phypartspiecharts.py /data/xiongtao/project/Rosaceae/Rosaceae_cytonuclear/ags353/results/03-treeshrink_more100/05-sptree/rosa_ags353_treeshrink_sp_rt.tre phyparts 353 --svg_name PhypartsPiecharts.svg --to_csv

date