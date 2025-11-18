for name in trimal/*.fasta;do
  pxlssq -s $name -m;
done|datamash mean 1 median 1

for name in ./*.fasta;do
  seqkit fx2tab -l -n $name|cut -f2|uniq;
done|datamash mean 1 median 1