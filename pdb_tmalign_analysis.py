#!/usr/bin/env python3
import os
import math
import argparse
import tempfile
import subprocess
import statistics
import logging
import time

from collections import defaultdict
from scipy import stats

from lib import allparsers


def header2dict(header, cond):

	ID = header.split("|")[0]
	Gene = header.split("|")[1]

	if cond == "WT":
		return {"header":header,"ID":ID,"Gene":Gene,"WT":True}

	parts = header.split("|")[-1].split("_")

	if len(parts) > 1:
		FSdir = parts[1]
	else:
		FSdir = "NA"

	condition = "aberrant"

	if "Ctrl_" in header:
		condition = "ctrl"

	if "exon" in header:
		condition = "mutation"

	if "_pos" in header:
		condition = "_".join([condition,header.split("_pos")[1].split("_")[0]])

	if "_FILEpos" in header:
		condition = "FILEpos"

	if "|" in header.split("_")[-2]:
		FSpos = header.split("_")[-2].split("|")[-1].split("_")[-2]
	else:
		FSpos = header.split("|")[-1].split("_")[-2]

	mca = header.split("_")[-1]

	if mca == "mRNA":
		FSpos = int(int(FSpos) / 3)
	else:
		FSpos = int(FSpos)

	return { "header":header, "ID":ID, "Gene":Gene, "FSdir":FSdir, "FSpos":FSpos, "mca":mca, "WT":False, "condition":condition}

def listallfiles(inputdir):

	maindict = {}

	for filename in os.listdir(inputdir):

		if filename.endswith("WT.pdb"):
			tmpdict = header2dict(filename[:-4], "WT")

		elif filename.endswith(".pdb"):
			tmpdict = header2dict(filename[:-4], "")

		else:
			continue

		tmpdict["filename"] = filename
		tmpdict["filepath"] = os.path.join(inputdir, filename)

		maindict[filename] = tmpdict

	return maindict

def createpairs(headerdict):

	pairs = []
	wt = {}

	for f in headerdict:

		if headerdict[f]["WT"]:
			wt[headerdict[f]["ID"]] = f

	for f in headerdict:

		if headerdict[f]["WT"]:
			continue

		if headerdict[f]["ID"] in wt:
			pairs.append((wt[headerdict[f]["ID"]], f))

	return pairs

###############################################3

def residue_confidences(pdbfile):

	out = defaultdict(list)

	with open(pdbfile) as f:

		for line in f:

			if not line.startswith("ATOM"):
				continue

			res = int(line[22:26])
			b = float(line[60:66])

			out[res].append(b)

	return {k:sum(v)/len(v) for k,v in out.items()}

def conf_summary(confdict):

	if len(confdict) == 0:
		return [None,None,None,None]

	vals = list(confdict.values())

	return [statistics.mean(vals), statistics.median(vals), min(vals), max(vals)]

def conf_region(confdict,start,end):

	tmp = {}

	for k,v in confdict.items():
		if start <= k <= end:
			tmp[k] = v

	return conf_summary(tmp)

def pdb_length(pdbfile):

	res = set()

	with open(pdbfile) as f:
		for line in f:
			if line.startswith("ATOM"):
				res.add(int(line[22:26]))

	return len(res)

def splitpdb(infile,outfile,start,end):

	with open(outfile,"w") as out:

		with open(infile) as f:

			for line in f:

				if not line.startswith("ATOM"):
					continue

				res = int(line[22:26])

				if start <= res <= end:
					out.write(line)

####################################################


def run_tmalign(pdb1,pdb2):

	try:

		result = subprocess.run( ["TMalign",pdb1,pdb2], capture_output=True, text=True)

	except Exception:
		return [None]*5

	tm1=None
	tm2=None
	rmsd=None
	alnlen=None

	for line in result.stdout.split("\n"):
		if line.startswith("Aligned length="):
			tmp=line.split(",")
			alnlen=int(tmp[0].split("=")[1].strip())
			rmsd=float(tmp[1].split("=")[1].strip())

		if "TM-score=" in line:
			score=float(line.split("=")[1].split()[0])

			if tm1 is None:
				tm1=score
			elif tm2 is None:
				tm2=score

#			if tm1 is None:
#				print("FAILED FULL ALIGNMENT")
#				print(wt["filepath"])
#				print(oof["filepath"])

	return tm1,tm2,rmsd,alnlen,result.returncode

#######################################


def ci95(values):

	if len(values) < 2:
		return None,None

	sem = stats.sem(values)
	mean = statistics.mean(values)
	ci = stats.t.ppf(0.975,len(values)-1)*sem

	return mean-ci, mean+ci



#################################################

def main():

	logging.basicConfig(level = logging.INFO , format='{asctime} - {levelname} - {message}',filemode='w', style = "{", datefmt = "[%d-%m-%Y] [%H:%M] ")
	starttime = time.time()
	logging.info("started pdb TM-align run")

	parser=argparse.ArgumentParser()
	allparsers.parser_get_pdb_scores(parser)
	args = parser.parse_args()

	headerdict=listallfiles(args.inputdir)
	pairs=createpairs(headerdict)

	results=args.outfilename+"_TM_results.tsv"

	rows=[]

	with open(results,"w") as out:

		header=[ "ID","Gene","condition","FSpos", "WT_length","OOF_length", "WT_after_length","OOF_after_length", "after_length_ratio", "TM_full_WTnorm","TM_full_OOFnorm", "TM_after_WTnorm","TM_after_OOFnorm", "TM_after_aligned_fraction", "RMSD_full","RMSD_after", "OOF_after_conf_mean", "delta_conf_after"]
		out.write("\t".join(header)+"\n")

		for wtf,ooff in pairs:

			wt=headerdict[wtf]
			oof=headerdict[ooff]

			fs=oof["FSpos"]

			tmf1,tmf2,rmsdf,aln,ret = run_tmalign( wt["filepath"], oof["filepath"])

			wtconf=residue_confidences(wt["filepath"])
 
			tmf1,tmf2,rmsdf,_,_=run_tmalign(wt["filepath"],oof["filepath"])

			oofconf=residue_confidences(oof["filepath"])

			wt_after_len=max(0,pdb_length(wt["filepath"])-fs)
			oof_after_len=max(0,pdb_length(oof["filepath"])-fs)

			with tempfile.TemporaryDirectory() as tmp:

				wta=os.path.join(tmp,"wta.pdb")
				oofa=os.path.join(tmp,"oofa.pdb")

				splitpdb(wt["filepath"],wta,fs,999999)
				splitpdb(oof["filepath"],oofa,fs,999999)

				tma1=tma2=rmsda=aln=None

				if pdb_length(wta)>=10 and pdb_length(oofa)>=10:
					tma1,tma2,rmsda,aln,_=run_tmalign(wta,oofa)

			wt_after_conf=conf_region(wtconf,fs,999999)
			oof_after_conf=conf_region(oofconf,fs,999999)

			ratio=None
			if wt_after_len>0:
				ratio=oof_after_len/wt_after_len

			alnfrac=None
			if aln and min(wt_after_len,oof_after_len)>0:
				alnfrac=aln/min(wt_after_len,oof_after_len)

			delta_conf=None
			if wt_after_conf[0] is not None and oof_after_conf[0] is not None:
				delta_conf=oof_after_conf[0]-wt_after_conf[0]

			row={"ID":oof["ID"], "Gene":oof["Gene"], "condition":oof["condition"], "TM_after_WTnorm":tma1, "TM_full_WTnorm":tmf1}
			rows.append(row)

			out.write("\t".join(map(str,[oof["ID"],oof["Gene"],oof["condition"],fs, pdb_length(wt["filepath"]), pdb_length(oof["filepath"]), wt_after_len,oof_after_len, ratio, tmf1,tmf2, tma1,tma2, alnfrac, rmsdf,rmsda, oof_after_conf[0], delta_conf ]))+"\n")

	statsfile=args.outfilename+"_TM_statistics.tsv"

	with open(statsfile,"w") as out:
		out.write("condition\tmetric\tn\tmean\tmedian\tstd\tCI95_low\tCI95_high\n")

		for metric in ["TM_after_WTnorm","TM_full_WTnorm"]:
			conds=sorted(set(r["condition"] for r in rows))

			for cond in conds:
				vals=[float(r[metric]) for r in rows if r["condition"]==cond and r[metric] is not None]

				if len(vals)==0:
					continue

				lo,hi=ci95(vals)
				std=statistics.stdev(vals) if len(vals)>1 else 0

				out.write(f"{cond}\t{metric}\t{len(vals)}\t{statistics.mean(vals)}\t{statistics.median(vals)}\t{std}\t{lo}\t{hi}\n")


	endtime = time.time() #
	endtime = endtime - starttime #
	endtime = round((endtime/60),2)
	printstring =  "pdb TM-align finished in " + str(endtime) + " m"
	logging.info(printstring )

if __name__=="__main__":
	main()
