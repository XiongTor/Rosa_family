# !/bin/bash
# Author: Tao Xiong
# Date: 2026-01-16
# Description: Keep only one outgroup in gene trees
# ==== 主体代码开始 ====
echo "Zelkova_schneideriana" >../old_z.txt;
echo "Elaeagnus_angustifolia" >../old_e.txt;
echo "Morus_indica" >../old_m.txt;
echo "Outgroup" >../new.txt;

for name in ./*.tre;do
  tt=$(basename $name .tre);   
  z=$(nw_labels -I $name|sort|uniq|grep -c "Zelkova_schneideriana");   
  e=$(nw_labels -I $name|sort|uniq|grep -c "Elaeagnus_angustifolia");   
  m=$(nw_labels -I $name|sort|uniq|grep -c "Morus_indica");   
  if [[ $z -eq 1 ]];then     
    rm=$(nw_labels -I "$name" | grep -E "Elaeagnus_angustifolia|Morus_indica" | paste -sd "," -);     
    echo $name "z";
    if [[ -n $rm ]];then
      pxrmt -t $name -n $rm>./${tt}_oneoutg.tre;     
      pxrlt -t ./${tt}_oneoutg.tre -c ../old_z.txt -n ../new.txt >./${tt}_oneoutg_final.tre;  
    else
      pxrlt -t ./${name} -c ../old_z.txt -n ../new.txt >./${tt}_oneoutg_final.tre
    fi   
  elif [[ $e -eq 1 ]];then     
    pxrmt -t $name -n Morus_indica>./${tt}_oneoutg.tre;     
    echo $name "e";     
    pxrlt -t ./${tt}_oneoutg.tre -c ../old_e.txt -n ../new.txt >./${tt}_oneoutg_final.tre;   
  else     
    echo $name "m";     
    pxrlt -t $name -c ../old_m.txt -n ../new.txt >./${tt}_oneoutg_final.tre;   
  fi; 
 done

#检查是否最终所有的树中都有Outgroup
 for tree in ./*_final.tre;do
   tt=$(grep -c "Outgroup" $tree);
   if [[ $tt -eq 0 ]];then
     echo $tree "has no outgroup!";
   fi
 done