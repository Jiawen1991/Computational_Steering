for j in {4..10}
do
for i in {1..4}
do
  ./stencil2d $i $j 
done
echo -e "\n"
done
