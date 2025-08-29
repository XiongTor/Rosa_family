for name in trimal/*.fasta;do
  pxlssq -s $name -m;
done|datamash mean 1 median 1

for name in trimal/*.fata;do
  seqkit fx2tab -l -n $name
done|datamash mean 2 median 2
