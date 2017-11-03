# usage: bash myscript.bash Spp_ID fasta_file_name
# example: bash Rename_fasta_genes.bash test Aedes-aegypti-Liverpool_PEPTIDES.fa

# This script will rename all genes in a fasta file according to the following convention: 
# >Spp_ID@0
# >Spp_ID@1
# >Spp_ID@2
# ...
#
# 1st agrument should be the species identifier
# 2nd argument should be the fasta file ending in ".fa"

# save arguments
Spp_ID=$1
FA_file=$2

FA_file_MOD=$(echo $FA_file | sed 's/.fa//g' | awk '{print $1"_MOD.fa"}')
FA_file_KEY=$(echo $FA_file | sed 's/.fa//g' | awk '{print $1"_gene_names.txt"}')

echo $FA_file_MOD $FA_file_KEY

echo -e "Creating mod fasta file using fa file $FA_file and spp ID $Spp_ID\n..."
<<<<<<< HEAD
awk -v var="$Spp_ID" '/^>/{print ">"var”|” ++i-1; next}{print}' $FA_file > $FA_file_MOD
echo -e "Complete"

echo -e "Creating key file named $FA_file_KEY\n..."
awk -v var="$Spp_ID" '/^>/{print ">"var”|” ++i-1; next}{print}' $FA_file | grep '>' > temp.txt
=======
awk -v var="$Spp_ID" '/^>/{print ">"var"|" ++i-1; next}{print}' $FA_file > $FA_file_MOD
echo -e "Complete"

echo -e "Creating key file named $FA_file_KEY\n..."
awk -v var="$Spp_ID" '/^>/{print ">"var"|" ++i-1; next}{print}' $FA_file | grep '>' > temp.txt
>>>>>>> 37dca9b4e0911f77309dc7e5924a81eada18d0e4
grep '>' $FA_file | paste temp.txt - > $FA_file_KEY
rm temp.txt
echo -e "Complete"
