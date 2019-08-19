for ((  i = 0 ;  i < 8;  i++  ))
do
  echo "Scanning spark # $i..."
  ./linescan -f ../output/Mar21/cru3d_grid -i $i -n 26 || break
done
