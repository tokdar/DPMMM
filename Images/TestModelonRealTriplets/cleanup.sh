START=1
              END=4
EXT=".pdf"
for ((i=START;i<=END;i++));
do
file="../../Figures/TestModelonRealTriplets/Triplet_$i/Summary_Triplet$i"
file+=".pdf"	
echo "Copying $file"
cp $file ./
  echo "Converting Summary_Triplet$i$EXT "
file="Summary_Triplet$i$EXT"
pdftoppm -rx 300 -ry 300 -png "$file" "$file"
done
awk -F, '{print("mv \"" $1 "\" \"" $2 "\"")}' rename.csv | bash -