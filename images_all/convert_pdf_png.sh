START=1
END=32
EXT=".pdf"
for ((i=START;i<=END;i++));
do
	file="../Figures_All/Triplet_$i/Summary_Triplet$i"
	file+=".pdf"	
	echo "Copying $file"
	cp $file ./
	echo "Converting Summary_Triplet$i$EXT "
	file="Summary_Triplet$i$EXT"
	pdftoppm -rx 300 -ry 300 -png "$file" "$file"
done
#for file in *.pdf
	#do pdftoppm -rx 300 -ry 300 -png "$file" "$file"
#done
