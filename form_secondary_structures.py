import tempfile
import sys
import argparse
import os
import logging
import time
import subprocess

from lib import allparsers
from lib import basictools

def header2dict(header, cond):

	ID = header.split("|")[0]
	Gene = header.split("|")[1]

	if cond == "WT":
		return {"header":header, "ID":ID, "Gene":Gene, "FSdir": None, "FSpos": None, "mca": None, "WT": True}
	else:

		FSdir = header.split("|")[2]
		descr = header.split("_")[-2]

		condition = "aberrant"

		if "Ctrl_" in header:
			condition = "ctrl"

		if "exon" in header:
			condition = "mutation"

		if "_pos" in header:
			condition = "_".join([condition,header.split("_pos")[1].split("_")[0]])

		if "_FILEpos" in header:
			condition = "FILEpos"
######################################

		if "|" in descr:
			FSpos = descr.split("|")[0]
		else:
			FSpos = header.split("|")[-1].split("_")[-2]

		mca = header.split("|")[-1].split("_")[-1]
		if mca == "mRNA":
			FSpos = int(FSpos)
			FSpos = str(int(FSpos/ 3))

		return {"header":header, "ID":ID, "Gene":Gene, "FSdir": FSdir, "FSpos": FSpos, "mca": mca, "WT": False, "condition":condition}

def read_secstrucfile(headerdict, filelines):

	struc_dict = {}
	struc_dict[headerdict["header"]] = {}
	struc_dict[headerdict["header"]] = {"struc":[], "conf":[]}

	filelines = list(filelines.split("\n"))

	i = 0

	for line in filelines:
		line = line.strip()

		if not line.startswith("#") and not line == "":
			pos = int(line.split(" ")[0])
			AA = line.split(" ")[1]

			struc = line.split(" ")[2]
			struc_dict[headerdict["header"]]["struc"].append(struc)

			Coil_confidence = float(line.split(" ")[5])
			Helix_confidence = float(line.split(" ")[7])
			TurnE_confidence = float(line.split(" ")[9])
			struc_dict[headerdict["header"]]["conf"].append([Coil_confidence,Helix_confidence,TurnE_confidence])

		i = i + 1

	return struc_dict


def matchconfidence(confidencelist, letterlist):

	matches = {"C":[], "H":[], "E":[]}

	for i in range(0,len(letterlist)):
		matches["C"].append(confidencelist[i][0])
		matches["H"].append(confidencelist[i][1])
		matches["E"].append(confidencelist[i][2])

	fracC = sum(matches["C"])/len(matches["C"])
	fracH = sum(matches["H"])/len(matches["H"])
	fracE = sum(matches["E"])/len(matches["E"])

	return fracC, fracH, fracE

def compare(compositiondict, WTdict, args):

	headerdict = header2dict(list(compositiondict.keys())[0], "")
	WTheaderdict = header2dict(list(WTdict.keys())[0], "WT")

	ab_structseq = compositiondict[list(compositiondict.keys())[0]]["struc"]
	ab_confseq = compositiondict[list(compositiondict.keys())[0]]["conf"]

	WT_structseq = WTdict[list(WTdict.keys())[0]]["struc"]
	WT_confseq = WTdict[list(WTdict.keys())[0]]["conf"]

	if len(ab_structseq) > int(headerdict["FSpos"]) and len(WT_structseq) > int(headerdict["FSpos"]):

		AB_conffrac = matchconfidence(ab_confseq, ab_structseq)
		WT_conffrac = matchconfidence(WT_confseq, WT_structseq)
		diffC = AB_conffrac[0] - WT_conffrac[0]
		diffH = AB_conffrac[1] - WT_conffrac[1]
		diffE = AB_conffrac[2] - WT_conffrac[2]

		aft_ab_structseq = ab_structseq[int(headerdict["FSpos"]):]
		aft_WT_structseq = WT_structseq[int(headerdict["FSpos"]):]
		aft_ab_confseq = ab_confseq[int(headerdict["FSpos"]):]
		aft_WT_confseq = WT_confseq[int(headerdict["FSpos"]):]

		aft_AB_conffrac = matchconfidence(aft_ab_confseq, aft_ab_structseq)
		aft_WT_conffrac = matchconfidence(aft_WT_confseq, aft_WT_structseq)
		aft_diffC = aft_AB_conffrac[0] - aft_WT_conffrac[0]
		aft_diffH = aft_AB_conffrac[1] - aft_WT_conffrac[1]
		aft_diffE = aft_AB_conffrac[2] - aft_WT_conffrac[2]

		if args.outputformat == "diff":
			return [str(diffC), str(diffH), str(diffE) ,str(aft_diffC), str(aft_diffH), str(aft_diffE)]
		else:
			return [str(AB_conffrac[0]),str(WT_conffrac[0]),str(WT_conffrac[1]),str(AB_conffrac[1]),str(WT_conffrac[2]),str(AB_conffrac[2]),str(aft_AB_conffrac[0]),str(aft_WT_conffrac[0]),str(aft_WT_conffrac[1]),str(aft_AB_conffrac[1]),str(aft_WT_conffrac[2]),str(aft_AB_conffrac[2])]


def loop_seqs(fastadict, args):

	outfilename = args.outfilename
	outfile = open(outfilename, "w")

	if args.outputformat == "diff":
		outfile.write("ID\tAApos\tGene\tdiff_Coil\tdiff_Helix\tdiff_E_turn\tafter_diff_Coil\tafter_diff_Helix\tafter_diff_E_turn\tOriginalID\tcondition\n")
	else:
		outfile.write("ID\tAApos\tGene\tAB_Coil\tWT_Coil\tAB_Helix\tWT_Helix\tAB_Eturn\tWT_Eturn\tOriginalID\tcondition\n")

	for header, seq in fastadict.items():

		if header.split("|")[-1] == "WT":
			if header.split("|")[-2] != "WT":

				WTheaderdict = header2dict(header, "WT")
				WTdict = {}
				base = "|".join(header.split("|")[:-1])

				WTdict[base] = seq

				temp = tempfile.NamedTemporaryFile(delete=False)
				basictools.write_fasta_to_temp(temp, WTdict, "")
				temp.seek(0)
				outstring = "python3 s4pred/run_model.py " + temp.name + " -t ss2"
				filelines = subprocess.getoutput(outstring)

				temp.close()

				WTcompositiondict = read_secstrucfile(WTheaderdict, filelines)
		else:
			base = "|".join(header.split("|")[:-1])

			headerdict = header2dict(header, "")

			ABdict = {}
			ABdict[header] = seq

			temp2 = tempfile.NamedTemporaryFile(delete=False)
			basictools.write_fasta_to_temp(temp2, ABdict, "")
			temp2.seek(0)
			outstring2 = "python3 s4pred/run_model.py " + temp2.name + " -t ss2"
			filelines2 = subprocess.getoutput(outstring2)

			temp2.close()

			compositiondict = read_secstrucfile(headerdict, filelines2)

			fracdiffs = compare(compositiondict, WTcompositiondict, args)
			headerdict = header2dict(list(compositiondict.keys())[0], "")

#			print(headerdict["condition"])

			if fracdiffs != None:  # write out.

				if args.outputformat == "diff":
					outline = "\t".join([headerdict["ID"][1:], headerdict["FSpos"],headerdict["Gene"],fracdiffs[0], fracdiffs[1], fracdiffs[2] ,fracdiffs[3], fracdiffs[4], fracdiffs[5], headerdict["header"][1:], headerdict["condition"]])
				else:
					outline = "\t".join([headerdict["ID"][1:], headerdict["FSpos"],headerdict["Gene"],fracdiffs[0], fracdiffs[1], fracdiffs[2] ,fracdiffs[3], fracdiffs[4], fracdiffs[5],fracdiffs[6], fracdiffs[7], fracdiffs[7] ,fracdiffs[9], fracdiffs[10], fracdiffs[11], headerdict["header"][1:], headerdict["condition"]])
				outfile.write(outline)
				outfile.write("\n")
	outfile.close()




def main():

	logging.basicConfig(level = logging.INFO , format='{asctime} - {levelname} - {message}',filemode='w', style = "{", datefmt = "[%d-%m-%Y] [%H:%M] ")
	logging.info("start run")

	if not os.path.exists("./s4pred/run_model.py"):
		logging.warning("s4pred not installed ! stop run")
		sys.exit()

	starttime = time.time()

	parser = argparse.ArgumentParser()
	allparsers.parser_get_secstruc(parser)
	args = parser.parse_args()
	fastadict = basictools.get_fastadict(args.input)

	loop_seqs(fastadict, args)

main()



