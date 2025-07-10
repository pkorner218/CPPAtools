import sys
import time
import argparse
import os
import re
import logging
import tempfile
import ast
from collections import Counter
import subprocess

from lib import basictools
from lib import allparsers


##############

def getempty_file(args, temp):

	infile = open(args.inputfastafile, "r")

	lines = infile.readlines()

	for line in lines:    
		line = line.strip()

		if line.startswith(">"):
			header = line
			condition = header.split("|")[-1].split("_")[0]

			if condition != "WT":
				ID = header.split("|")[0][1:]
				Gene = header.split("|")[1]
				abpos = str(getposfromheader(header)[1])

				outstring = "\t".join([header[1:],ID,Gene,abpos])

				temp.write(str(outstring).encode("utf-8"))
				temp.write(b"\n")

###############



def mappos(abpos, start, end):

	if start < abpos < end: #### abpos is in middle of domain
		part = "spanning"
	elif abpos > end: ### abpos is behind the domain # domain is "before" abpos
		part = "before"
	elif abpos <= start: ### abpos is before the domain # domain is "after" abpos
		part = "after"
	elif abpos == end:
		part = "truncated" ### abpos is endpoint of the domain # it potentially truncates domain

	return part


###############

#
def loop_fill(args,temp):

	outfile = open(args.outfilename, "w")
	outfile.write("OriginalID\tID\tGene\taberrantpos\tbefore\tbeforeWT\tspanning\tspanningWT\ttruncated\ttruncatedWT\tafter\tafterWT\n") 

	temp.seek(0)
	infilename = open(temp.name, "r")
	lines = infilename.readlines()

	for line in lines:
		line = line.strip()

		abpos = line.split("\t")[3]
		origID = line.split("\t")[0]
		WTID = origID.replace(origID.split("|")[-1],"WT")

		p = subprocess.run(["grep", origID, args.inputdomainfile], capture_output = True)
		px = p.stdout.decode("utf-8")
		pxl = px.split("\n")[:-1]

		pw = subprocess.run(["grep", WTID, args.inputdomainfile], capture_output = True)
		pxw = pw.stdout.decode("utf-8")
		pxwl = pxw.split("\n")[:-1]

		emptyscores = "0\t0\t0\t0\t0\t0\t0\t0"

		if pxl != []:
			for domainline in pxl:
				start = domainline.split("\t")[6]
				end = domainline.split("\t")[7]
				part = mappos(int(abpos), int(start), int(end))
				emptyscores  = level_counter(part, emptyscores, "Ab")

		if pxwl != []:
			for domainline in pxwl:
				start = domainline.split("\t")[6]
				end = domainline.split("\t")[7]
				part = mappos(int(abpos), int(start), int(end))

				emptyscores = level_counter(part,emptyscores,  "WT")

		outstring = "\t".join([line,emptyscores])
		outfile.write(outstring)
		outfile.write("\n")        

	outfile.close()






def level_counter(part, line, condition):

	partdict = {"before":0,"spanning":2,"truncated":4,"after":6}

	if condition != "WT":
		partindex = int(partdict[part])
	else:
		partindex = int(partdict[part]) + 1

	xline = line.split("\t")
	xline[partindex] = str(int(xline[partindex]) + 1)
	xline2 = "\t".join(xline)

	return xline2


######################

def getposfromheader(header):

	condition = "aberrant"

	if "_Ctrl" in header:
		condition = "ctrl"

	elif "_pos" in header:
		condition = "_".join([condition,header.split("_pos")[1].split("_")[0]])
	elif "_Filepos" in header:
		condition = "_".join([condition,header.split("_Filepos_")[1].split("_")[0]])

	if "exon" in header: ### if its mt based
		mcalevel =  header.split("_")[-1]
		pos = int(header.split("|")[-1].split("_")[-2].replace("Ctrl", ""))
	else:
		mcalevel =  header.split("_")[-1]
		pos = int(header.split("|")[-1].split("_")[-2].replace("Ctrl", ""))
                        
	if mcalevel == "mRNA":
		pos = int(pos/3)

	return header, pos, mcalevel, condition

##############


def main():


	parser = argparse.ArgumentParser()
	allparsers.parser_get_domNr(parser)
	args = parser.parse_args()
	print(args)
	print("")

	logging.basicConfig(level = logging.INFO , format='{asctime} - {levelname} - {message}',filemode='w', style = "{", datefmt = "[%d-%m-%Y] [%H:%M] ")
	logging.info("start run")

	starttime = time.time()

	temp = tempfile.NamedTemporaryFile(delete=False)

	getempty_file(args, temp)
	loop_fill(args,temp)

	temp.close()
 
	endtime = time.time() - starttime
	endtime = round((endtime/60),2)
	infotext2 = "run finished in " + str(endtime) +  " m"
	logging.info(infotext2)


	print("")

main()
