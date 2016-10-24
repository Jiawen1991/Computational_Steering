for j in {4..5}
do
for i in {1..70}
do
  ./stencil2d $i $j 
done
echo -e "\n"
done
