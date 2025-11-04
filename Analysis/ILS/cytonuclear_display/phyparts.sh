    
#phyparts分析
/home/xiongtao/data/tree/rosa_tree/rosa_best_gene_tree

java -jar /data/xiongtao/software/phyparts/target/phyparts-0.0.1-SNAPSHOT-jar-with-dependencies.jar -a 1 -v -d ../../04-genetrees_reroot -m ../../05-sptree/rosa_ags353_treeshrink_sp_rt.tre -o ./phyparts 

# remember to change the path of phyparts_pies.sh
bash /data/xiongtao/scripts/phyparts_pies.sh