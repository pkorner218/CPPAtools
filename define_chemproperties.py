import sys
import time
import argparse
import os
import peptides
import logging
import numpy as np

from lib import basictools
from lib import allparsers
from multiprocessing import Process, Queue

def exec_function(input_key, sequence, header, qout):

	peptide = peptides.Peptide(sequence)

	func_dict = {
	"aliphatic_index":peptide.aliphatic_index(),
	"boman":peptide.boman(),
	"charge":peptide.charge(),
	"hydrophobic_moment":peptide.hydrophobic_moment(),
	"hydrophobicity":peptide.hydrophobicity(),
	"instability_index":peptide.instability_index(),
	"isoelectric_point":peptide.isoelectric_point(),
	"molecular_weight":peptide.molecular_weight(),
	"mz": peptide.mz()
	}

	qout.put({input_key:func_dict[input_key]})


##################

def getprocesses(sequence, props, header, qout):

	processes  = [Process(target=exec_function, args=(prop, sequence, header, qout)) for prop in props]

	for p in processes:
		p.start()
	for p in processes:
		p.join()

	properties = [qout.get() for p in processes]

	indict = {}

	for prop in properties:
		key = list(prop.keys())[0]
		value = list(prop.values())[0]
		indict[key] = str(value)
        
	return indict


##################


def loop_seqs(fastadict, WTdict, args):

	qout = Queue()
	props = get_chemical_properties(args)  #basictools.chemicalproperties

	if "/" in args.outfilename:
		chemoutname = args.outfilename
	else:
		cwd = os.getcwd()
		chempath = "".join([cwd,"/chem/"])

		outfilename = args.outfilename.split("/")[-1]
		chemoutname = "".join([chempath,outfilename,"_chem_properties.txt"])

	infotext = "Results written to " + chemoutname
	logging.info(infotext)

	outfile = open(chemoutname, "w")
	outfile.write("ID\tAApos\tGene\tchemproperty\tfullOOF\tfullWT\tafterOOF\tafterWT\tperc10OOF\tperc10WT\tcondition\tOriginalID\n")

	for header, sequence in fastadict.items():
		condition = "aberrant"
 
		if header.split("|")[-1] != "WT":
			if header.split("|")[-2] != "WT":

				originalheader = header

				if "Ctrl_" in header:
					condition = "ctrl"

				if "exon" in header:
					condition = "mutation"
					mcalevel = header.split("_")[-1]
					pos = int(header.split("|")[-1].split("_")[-2])
					header  = "|".join(header.split("|")[:-2]) + "|CDS|CDS"

				if "_pos" in header:
					condition = "_".join([condition,header.split("_pos")[1].split("_")[0]])
					mcalevel =  header.split("_")[-1]
					pos = int(header.split("|")[-1].split("_")[-2].replace("Ctrl", ""))

					if mcalevel == "mRNA":
						pos = int(pos/3)

				if "_FILEpos" in header:
					condition = "FILEpos"
					mcalevel = header.split("|")[-1].split("_")[-1]
					if "." in header.split("|")[-1].split("_")[-2]:
						pos = int(header.split("|")[-1].split("_")[-2].split(".")[0])
					else:
						pos = int(header.split("|")[-1].split("_")[-2])

					if mcalevel == "mRNA":
						pos = int(pos/3)

				ID = header.split("|")[0][1:]
				Gene = header.split("|")[1] 

				if "." in header.split("|")[:-1]:
					FShead = "|".join(header.split("|")[:-1].split(".")[0])
				else:
					FShead = "|".join(header.split("|")[:-1])

				if pos < len(sequence) and pos < len(WTdict[FShead][0]): ### otherwiseFS happened at last AA
					event = "_".join(header.split("|")[-1].split("_")[0:2])
					base = "|".join(header.split("|")[:-1])

					WTseq = WTdict[base][0] #["sequence"][0]

					before = sequence[:pos]
					after = sequence[pos:]
					WTbefore = WTseq[:pos]
					WTafter = WTseq[pos:]
 
					totalAB = getprocesses(sequence, props, header, qout)
					totalWT = getprocesses(WTseq, props, header, qout)

					afterAB = getprocesses(after, props, header, qout)
					afterWT = getprocesses(WTafter, props, header, qout)

					WTperc10 = int((len(WTseq) / 100) *  90)
					ABperc10 = int((len(sequence) / 100) *  90)
					WTperc10seq = WTseq[WTperc10:]
					ABperc10seq = sequence[ABperc10:]

					perc10AB = getprocesses(ABperc10seq, props, header, qout)
					perc10WT = getprocesses(WTperc10seq, props, header, qout)

					for prop in props:
						outstring = "\t".join([ID,str(pos),Gene,prop,totalWT[prop],totalAB[prop],afterAB[prop],afterWT[prop],perc10AB[prop],perc10WT[prop],condition,originalheader[1:]])
						outfile.write(outstring)
						outfile.write("\n")

###     ID      pos     Gene    chemproperty    fullOOF fullWT  afterOOF        afterWT        condition      originalID
#       ENST    x       Xx      charge  -0.2    0.2     -0.6    0.1        aberrant          ENST|Gene|FS //...
#       ENST    x       Xx      hydrophobic  -12    2     -16    -2        ctrl          ENST|Gene|WT //...


def get_chemical_properties(args):

	if args.properties == "all":
		chemicalproperties = ["aliphatic_index","boman","charge","hydrophobic_moment","hydrophobicity","instability_index","isoelectric_point","molecular_weight","mz"]
	else:
		chemicalproperties = list(args.properties.split(","))

	return chemicalproperties

def main():

	logging.basicConfig(level = logging.INFO , format='{asctime} - {levelname} - {message}',filemode='w', style = "{", datefmt = "[%d-%m-%Y] [%H:%M] ")
	logging.info("start run")

	starttime = time.time()

	parser = argparse.ArgumentParser()
	allparsers.parser_get_chem(parser)
	args = parser.parse_args()

	fastadict = basictools.get_fastadict(args.input)
	WTdict = basictools.obtain_basedicts(fastadict)

	loop_seqs(fastadict, WTdict, args)

	endtime = time.time() - starttime
	endtime = round((endtime/60),2)

	infotext = "run finished in " + str(endtime) +  " m"
	logging.info(infotext)

main()
