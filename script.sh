for j in {3..4}
do
for i in {1..20}
do
  ./stencil2d $i $j 
done
echo -e "\n"
done
