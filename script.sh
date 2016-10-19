for j in {20..21}
do
for i in {1..20}
do
  ./stencil2d $i $j 
done
echo -e "\n"
done
