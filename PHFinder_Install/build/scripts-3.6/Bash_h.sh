#! /bin/bash

if [[ $1 == 0 ]] ; then
	bowtie2-build -q $2 $3
	list=$(find -type f -name "*.fastq")
	echo "Alinging sequences..."
	for i in $list
		do
			out=$(echo $i | sed -E 's/(\w+\/.+\.)\w+/\1sam/' )
			bowtie2 --local --quiet -q $i -x $3 -S $out
		done
fi

if [[ $1 == 1 ]] ; then
	allabi=$( find -type f -name "*.ab1" | sed -e 's/\.\///g' | sed -e 's/\//\t/g')
	sample=$( find -type f -name "*.ab1" | sed -e 's/\.\///g' | sed -e 's/\//\t/g' | cut -f1 | sort -n | uniq)
	numbersamples=$(echo "$sample" | wc -l)
	cd ./$2
	samplesdetect=$(grep -E "\w+" Detection.txt | sed -e 's/\//\t/g' | cut -f1 | sort | uniq)
	data="$(cat ./NoDetection.txt | sed -e 's/\//\t/g' )"
	echo "Sample ID	AB1 file	Position	Main ratio	Average quality" > Final.txt
	for i in $samplesdetect
		do
			grep -E "\w+" Detection.txt | sed -e 's/\//\t/g' | awk -v id="$i" '$1 ==id' | sort	-k3n -k4n >> Final.txt
			grep -E "\w+" NoDetection.txt | sed -e 's/\//\t/g' | awk -v id="$i" '$1 ==id' |grep "no detection" >> Final.txt
			echo "######" >> Final.txt
		done
	COUNTER=$((COUNTER))
	printf "Samples in which no chromatograms were analysed:\n" >> log.txt
	for i in $sample
		do
			info=$(echo "$allabi" | grep "\b$i\b" | wc -l)
			noinfo=$(echo "$data" | grep "\b$i\b" | grep -E "Error|quality" | wc -l)
			if [[ "$info" = "$noinfo" ]]; then
				COUNTER=$((COUNTER + 1))
				echo "$i" >> log.txt
			fi
		done
	printf "\n$COUNTER from a total of $numbersamples" >> log.txt
fi
