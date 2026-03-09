    
#phyparts分析
#-Xmx50g 设置最大内存50G,防止处理大批量基因树时内存溢出
java -Xmx600g -jar /data/xiongtao/software/phyparts/target/phyparts-0.0.1-SNAPSHOT-jar-with-dependencies.jar -a 1 -v -d rosa_orthofinder_MO_treeshrink_genetrees.tre -m rosa_orthofinder_MO_treeshrink_sp_rt.tre -o ./phyparts 

# remember to change the path of phyparts_pies.sh
bash /data/xiongtao/scripts/phyparts_pies.sh