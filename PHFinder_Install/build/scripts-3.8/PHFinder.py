#!/usr/bin/python3

"""PHFinder: Assisted detection of point heteroplasmy in Sanger
sequencing chromatograms.
Copyright (C) 2020 Marcos Suarez

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation (version 3 of the License)

This program is distributed in the hope that it will be useful,
but without any warranty; without even the implied warranty of
merchantability or fitness for a particular purpose. See the
GNU General Public License for more details.

Version=1.0
Marine Evolution and Conservation
University of Groningen
GitHub: https://github.com/MSuarezMenendez/PHFinder
Contact: marcos.sume@gmail.com
"""


import argparse
from argparse import ArgumentParser
import sys
import re
import os
from os import listdir
import subprocess
from collections import defaultdict
import datetime
import warnings
from Bio.Seq import Seq
import numpy
import Bio
from Bio import SeqIO

parser = ArgumentParser()
flag = parser.add_argument_group('Arguments')
flag.add_argument("-f", action="store_true", dest="Fastq", help="Extract fastq\
 from AB1 files")
flag.add_argument("-a", action="store", dest="Align", default=False, help="Align fastq\
 to reference sequence (provide fasta file)")
flag.add_argument("-d", action="store_true", dest="Detection", help="Detect\
 possible heteroplasmies")
flag.add_argument("-t", action="store_true", dest="Test",default=False, help="Test run of parameters")
flag.add_argument("-o", action="store", dest="Output", help="Name\
 of output directory")
flag.add_argument("-r", action="store", type=int, default=20, dest="Ratio", help="Ratio\
 threshold of heteroplasmy (Default: 20)")
flag.add_argument("-q", action="store", type=int, default=40, dest="Qthreshold", help="Quality\
 threshold (Default: 40)")
flag.add_argument("-s", action="store", type=int, default=1, dest="Start", help="Position\
 where the analysis starts (Default: 1)")
flag.add_argument("-e", action="store", type=int, dest="End", help="Position\
 where the analysis ends (Required when '-d' is provided)")
flag.add_argument("-sr", action="store", type=float, default=0.4, dest="SRatio", help="Secondary ratios (Default:0.4)")

if len(sys.argv) < 2: #If no arguments are provided, help is printed
	sys.stderr.write("PHfinder: Assisted detection of point heteroplasmy in Sanger sequencing chromatograms (AB1 files).\n\
Marine Evolution and Conservation\
\nUniversity of Groningen\nGitHub: https://github.com/MSuarezMenendez/PHFinder\nContact: marcos.sume@gmail.com\n\n")
	parser.print_help()
	sys.exit(1)
args = parser.parse_args()

if args.Detection is True and (args.End is None):
	parser.error("-d requires -e.")
if args.Test is True and (args.Align is None):
	parser.error("-t requires -a.")
if args.Test is True and (args.End is None):
	parser.error("-t requires -e.")
if args.Test is None and (args.Output is None):
	parser.error("argument -o is required")
if args.Test is True:
	if os.path.isfile("List_HP.csv") is not True:
		parser.error("Input file List_HP.csv is required")
#Test run
if args.Test is True:
	subprocess.call(['Testing.sh' , args.Align, str(args.Start), str(args.End)])

#Output files
Directory = args.Output
Detection = "/Detection.txt"
Log = "/log.txt"
NoDetection = "/NoDetection.txt"

#Parameters
start = args.Start #position where the analysis starts
end = args.End # position where the analysis ends
Threshold = args.Qthreshold #Quality threshold
MainRatio = args.Ratio # Heteroplasmy ratio threshold
#Second peaks surrounding the heteroplasmy position need to be lower than the
#specified percentage of the main second peak
Peak1b = args.SRatio #first peak before the heteroplasmy position
Peak2b = args.SRatio
Peak3b = args.SRatio
Peak1a = args.SRatio #first peak after the heteroplasmy position
Peak2a = args.SRatio
Peak3a = args.SRatio

#Directory where output will be in
try:
	os.mkdir(Directory)
except:
	print("Output directory already exists")
	quit()

nohp = open(Directory + NoDetection, "a") #Output file for no detections/errors
log = open(Directory + Log, "a")
now = datetime.datetime.now()
print ("Started %s" % now)
log.write("Started: %s\n" % now)

#Function that gets the paths to files with specific extension
def pathtofiles(Allfiles, arg1):
	count = 0
	for root, dirs, files in os.walk("./"):
		for f in files:
			path = os.path.relpath(os.path.join(root, f), "./")
			if path.endswith("%s" % arg1):
				Allfiles.append(path)
				count += 1
#Setting counters for fastq extraction
error= 0
extracted = 0

#Fastq extraction from ab1 files
if args.Fastq is True:
	#Extracting fastq
	print ("Extracting fastq from ab1 files...")
	log.write("Fastq extraction\n")
	Allfiles = list() #Necessary for the next function
	pathtofiles(Allfiles, arg1='.ab1')
	Numberab1 = len(Allfiles) #Total number of AB1 files
	#Creates a list with the output of the conversion from ab1 to fastq
	fastq = [ext.replace(".ab1", ".fastq") for ext in Allfiles] #Change extension
	for filename, output in zip(Allfiles, fastq): #Loops through both lists at the same time
		try:
			SeqIO.convert(filename, "abi", output, "fastq-sanger") #Convert ab1 to fastq
			extracted += 1
		except: #If there is an error converting the ab1 files:
			nohp.write("%s\tError extracting\n" % filename)
			os.remove(output) #Removes the empty fastq file
			error += 1
			pass
	log.write("%d fastq extracted out of %d AB1 files\n" % (extracted, Numberab1))
#This conditional ensures all log files include extraction errors
if args.Fastq is False:
	Listaab1 = list()
	pathtofiles(Listaab1, arg1='.ab1')
	Fastqfiles = [ext.replace(".ab1", ".fastq") for ext in Listaab1]
	for fastq, ab1 in zip(Fastqfiles, Listaab1):
		try:
			fh = open(fastq)
		except:
			nohp.write("%s\tError extracting\n" % ab1)
			error += 1
			pass

#Calls section of bash script that aligns the fastq to the provided reference
if args.Align is not False:
	if os.path.isfile(args.Align):
		refname = args.Align.split(".")
		refname = refname[0]
		try:
			subprocess.call(['Bash_h.sh', "0", args.Align, refname])
		except:
			subprocess.call(['./Bash_h.sh', "0", args.Align, refname])
	else:
		print("Error: Fasta file for alignment not found in directory")


Fastq = list()
pathtofiles(Fastq, arg1='.fastq') #path to all the generated fastq
Abifiles = [ext.replace(".fastq", ".ab1") for ext in Fastq]
Sam = [ext.replace(".fastq", ".sam") for ext in Fastq]

if args.Detection is True:
	log.write("Detection of heteroplasies\n")
	hp = open(Directory + Detection, "a")
	log.write("Parameters used:\nAnalysis starts in position: %d\n" % start)
	log.write("Analysis ends in position: %d\n" % end)
	log.write("Main ratio:%d\nQuality threshold: %d\nPercentage surrounding \
	peaks:\n" % (MainRatio, Threshold))
	log.write("Previous\tPosterior\n1: %.2f\t1: %.2f\n2: %.2f\t2: %.2f\n3: \
	%.2f\t3: %.2f\n" % (Peak1b,Peak1a,Peak2b,Peak2a,Peak3b,Peak3a))
	Total = len(Abifiles)
	print ("Files to be analysed: %d" % Total)
	print ("Analysing...")
	#Counters used along the script
	analysed = 0
	possibilities = 0
	lowquality = 0
	Averagequality = list()
	#Function that gets the value of the second highest peak
	def peakhight(peak):
		global Lowerpeak
		Mainpeak = peak.split("\t")
		Mainpeak = [int(x) for x in Mainpeak]
		Mainpeak.sort()
		Lowerpeak = Mainpeak[2]

	for chroma, align in zip(Abifiles, Sam): #Loop over all the ab1
		start = args.Start # reset every loop
		try:
			counter = 0
			traceinfo = list()
			#Extracting ab1 file information
			record = SeqIO.read(chroma, "abi")
			record.annotations.keys()
			record.annotations['abif_raw'].keys()
			trace = defaultdict(list)
			channels = ['DATA9', 'DATA10', 'DATA11', 'DATA12', 'PLOC2', 'PBAS2', 'PCON2']
			for c in channels:
				trace[c] = record.annotations['abif_raw'][c]
			#Peak height of each nucleotide
			NucleotideG = trace['DATA9']
			NucleotideA = trace['DATA10']
			NucleotideT = trace['DATA11']
			NucleotideC = trace['DATA12']
			Infopoint= trace['PLOC2'] #Right trace data points (used in base calling)
			Quality = trace['PCON2'] #Quality ascii characters
		   #Base called sequence
			Basecalls = trace['PBAS2'] #Basecalled sequence
			Length = len(Basecalls)
			#Alingment information
			with open(align, "r") as sam:
				content = sam.readlines()
				sam.close()
			info = content[3]
			info = info.split("\t")
			cigar = info[5]
			cigar = re.sub("(\d+)S\w+", r"\1", cigar)
			if info[3] == "1": #Ab1 contains whole reference sequence
				inicio = 0
			else:
				nohp.write("%s\tsequence not complete\n" % chroma)
				error += 1
				continue

			#Extracting trace data from exact positions to analyse
			if info[1] == "0": #Ab1 from forward strand
				Reverse = 0
				extra = int(cigar) #Position in ab1 before start of reference
				Location = 0
				for i, (A, T, G, C) in enumerate(zip(NucleotideA, NucleotideT, NucleotideG, NucleotideC)):
					for position in Infopoint:
						if position == i:
							Location += 1
							if Location not in range(extra + 1):
								traceinfo.append("%d\t%d\t%d\t%d" % (A, T, G, C))

			if info[1] == "16": #Ab1 from reverse strand
				extra = int(cigar)
				Reverse = 1
				Location = 0
				for i, (A, T, G, C) in enumerate(zip(NucleotideA, NucleotideT, NucleotideG, NucleotideC)):
					for position in Infopoint:
						if position == i:
							Location += 1
							if Location in range(Length - (extra - 1)):
								traceinfo.append("%d\t%d\t%d\t%d" % (A, T, G, C))
				traceinfo = traceinfo[::-1]
			if len(traceinfo) < (end - start): #Ab1 does not have complete sequence
				nohp.write("%s\tsequence not complete\n" % chroma)
				error += 1
				continue

			#Mean quality of file
			with warnings.catch_warnings():
				warnings.simplefilter("ignore", category=RuntimeWarning) #avoids warning when quality is 0
				qualitylist = list()
				if Reverse == 0:
					for i, (symbol) in enumerate(Quality):
						if i in range(extra - 1 + start,(extra + end - inicio)):#quality of the region analysed
							value = ord(symbol) #Transforms ascii into numerical values
							qualitylist.append(value)
					Meanquality = numpy.mean(qualitylist)
					if Meanquality >= Threshold:
						Averagequality.append(Meanquality)
				if Reverse == 1:
					for i, (symbol) in enumerate(Quality):
						if i in range(0, Length - (extra + start - 1)):#quality of the region analysed
							value = ord(symbol)
							qualitylist.append(value)
							final = end - inicio
						qualitylist1 = qualitylist[::-1][:final]
						Meanquality = numpy.mean(qualitylist)
					if Meanquality >= Threshold:
						Averagequality.append(Meanquality)
			#Filtering possibilities
			if Meanquality >= Threshold:
				analysed += 1
				for i, (Peak) in enumerate(traceinfo):
					if i in range(start, end - inicio):
						try:
							#Trace data surrounding studied postion
							priorpeak = traceinfo[i-1]
							priorpeak1 = traceinfo[i-2]
							priorpeak2 = traceinfo[i-3]
							posteriorpeak = traceinfo[i+1]
							posteriorpeak1 = traceinfo[i+2]
							posteriorpeak2 = traceinfo[i+3]
							#Main peak
							Mainpeak = Peak.split("\t")
							Mainpeak = [int(x) for x in Mainpeak]
							Mainpeak.sort()
							Highest = Mainpeak[3]
							Lower = Mainpeak[2]
							Ratio = (Lower*100)/Highest
							#prior1
							peakhight(priorpeak)
							Lower1 = Lowerpeak
							#prior2
							peakhight(priorpeak1)
							Lower12 = Lowerpeak
							#prior3
							peakhight(priorpeak2)
							Lower13 = Lowerpeak
							#posterior1
							peakhight(posteriorpeak)
							Lower2 = Lowerpeak
							#posterior2
							peakhight(posteriorpeak1)
							Lower21 = Lowerpeak
							#posterior2
							peakhight(posteriorpeak2)
							Lower22 = Lowerpeak
							#Conditions that have to be met to consider a double peak an heteroplasmy
							if Ratio >= MainRatio and Lower1 <= (Lower * Peak1b) and Lower12 <= (Lower * Peak2b)\
							and Lower13 <= (Lower * Peak3b)	and Lower2 <= (Lower * Peak1a)\
							and Lower21 <= (Lower * Peak2a) and Lower22 <= (Lower * Peak3a) :
								position = i + 1 + inicio
								counter += 1
								possibilities += 1
								#print "funciona"
								if Reverse == 0:
									hp.write("%s\t%d\t%d\t%d\n" % (chroma, position, Ratio, Meanquality))
								if Reverse == 1:
									hp.write("%s\t%d\t%d\t%d\tReverse\n" % (chroma, position, Ratio, Meanquality))
						except:
							pass
					else:
						continue
				if counter == 0: #If no possibilities found in the whole file
					#print "File %s no heteroplasmy detected" % filename
					nohp.write("%s\tno detection\n" % chroma)
			else:
				nohp.write("%s\tLow quality\t%d\n" % (chroma, Meanquality))
				lowquality += 1
		except Exception as e:
			pass
			nohp.write("%s\tError\t(%s)\n" % (chroma, e))
			error += 1

	Averagequality = numpy.nanmean(Averagequality)
	now = datetime.datetime.now()
	print ("Finished %s" % now)
	log.write("Finished: %s\n" % now)
	log.write("Number of chromatograms to be analysed: %d\n" % Total)
	log.write("Analysed chromatograms: %d\n" % analysed)
	log.write("Average quality analysed chromatograms: %.1f\n" % Averagequality)
	log.write("Quality lower than threshold in %d chromatograms\n" % lowquality)
	log.write("Possible heteroplasmies: %d\n" % possibilities)
	log.write("Not analysed chromatograms (errors): %d\n" % error)
	log.close()
	hp.close()
	nohp.close()
	try:
		subprocess.call(['Bash_h.sh', "1", Directory])
	except:
		subprocess.call(['./Bash_h.sh', "1", Directory])
