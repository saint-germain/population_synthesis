NAME=$1
CRE=createdfiles.txt
TOT=totalsims.txt
FIN=finished.txt
RES=finalresults.csv
QDI=qdiag.txt
LOG=log.txt
mkdir $NAME
cd $NAME
unzip "../"$NAME".zip" > $LOG
ls results/res* | xargs cat > temp.txt
cat results/header.txt temp.txt > $RES
rm temp.txt
cp results/*.o* $NAME.o
ls -lt results/rescortos* | awk '{print $9" " $5}' | cut -c18- | sort -n | sed s/".txt"//g > $CRE
ls results/info* | xargs -n1 | cut -c13- | sed s/".txt"// > $TOT
ls $NAME.o | xargs cat | grep -o -P  ".{0,31}run" | sort -n | awk '{print $1}' | uniq > $FIN
python ../quality.py $CRE $TOT $FIN > $QDI
python ../reduce.py $RES >> $QDI
cat $QDI
echo "***************************************"
zip $NAME"_results.zip" $RES terrestrial.txt giant.txt $QDI header.txt parameters.in >> $LOG
echo "Zip logs in "$LOG
echo "***************************************"
echo "Quality diagnostics output in "$QDI
echo "***************************************"
echo "Processed results in "$NAME"_results.zip"
echo "***************************************"
