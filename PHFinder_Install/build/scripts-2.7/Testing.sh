#! /bin/bash

mkdir ./Test
HP_Detect.py -f -a $1 -o "./Test/Extraction_Alignment"
for i in {15,20,25,30}
do
	for a in {0.2,0.3,0.4,0.5}
	do
		for b in {30,40,50,60}
		do
			HP_Detect.py -d -r $i -sr $a -q $b -s $2 -e $3 -o "./Test/MR$i-SR$a-AQ$b"
			echo "MR$i-SR$a-AQ$b" >> ./Test/Folderlist.txt
		done
	done
done

cp List_HP.csv ./Test
cd Test
touch Results.txt
echo "Parameters	HP_AB1	HP_Samples	False_Positives	Analysed_Samples	Analysed_AB1	Open_AB1" >> Results.txt
while read line
do
	echo $line
	cd $line
	COUNTER=$((COUNTER))
	while read sentence
	do
		echo $sentence
		if cat Final.txt | grep -q "$sentence"; then
			COUNTER=$((COUNTER + 1))
			echo $sentence >> temporal.txt
		fi
	done < ../List_HP.csv
	Dsamples=$(cat temporal.txt | cut -f1 -d" " | sort | uniq | wc -l) #Samples where HP was detected
	rm temporal.txt
	Asamples=$(cat log.txt | grep "from a total" | cut -d" " -f1,6 | awk '{print $2 - $1}') #samples analysed
	declare -i possibilities=$(cat log.txt | grep "Possible" | cut -d " " -f3) #Number of posssible HP
	Aab1=$(cat log.txt | grep "Analysed chromatograms" | cut -d " " -f3)
	falsepositives=$(( possibilities - $COUNTER ))
	AB1=$(cat Final.txt | grep ab1 | grep -v "no detection" | cut -f2 | sort | uniq | wc -l) #Numer of AB1 that would need checking
	cd ..
	echo "$line	$COUNTER	$Dsamples	$falsepositives	$Asamples	$Aab1	$AB1" >> Results.txt
	COUNTER=$((0))
done < Folderlist.txt
