# do krona visualisation
#  runKrona.sh $infile $mothurtax $outputfile

inputfile=$1
isMothur=$2

# preprocess if taxonomy from MOTHUR
echo "starting"
head $inputfile

if [[ $isMothur == "Y" ]]
then
  awk 'BEGIN{FS="\t";OFS="\t"}{print $2,$3}' $inputfile > tempinput
  head tempinput
  sed -i 's/;/\t/g' tempinput
  head tempinput
 
else
  cp $inputfile tempinput
  
fi

head tempinput

# run Krona
echo $PATH
ktImportText tempinput


