import time
import re
import argparse
import itertools
import sys
import os
import math
import random
import tempfile
import logging
import json

from lib import basictools
from lib import allparsers
from lib import aberranttranslation
from lib import call_transvar

def is_vcf_mode(args):
	return str(args.translation_error).replace("-", "_") == "vcf_file"

def exec_function(input_key, *params):

	func_dict = {
	"frameshift":aberranttranslation.frameshift,
	"repeat":aberranttranslation.repeat,
	"skip_deletion":aberranttranslation.skip,
	"truncation":aberranttranslation.truncate,
	"alternative_start":aberranttranslation.altstart,
	"insertion":aberranttranslation.insertion,
	"reversion":aberranttranslation.reversion,
	}

	func = func_dict.get(input_key)

	return func(*params)


def validate_TSV_REF(REFdict,infilename, reftype):

	ID_REFdict = {}
	ID_tr_REFdict = {}
	ID_gene_REFdict = {}

	trREFids = []
	geneREFids = []

	validTSVids = {}
	nonvalidTSVids = []

	if reftype == "given": # check that given ref file has unique headers

		duplicated_headers = basictools.duplicated_seqs(REFdict.keys()) # find duplicated headers !!!

		if duplicated_headers != []:
			print("Fatal ! headers ", duplicated_headers, "were found more than once in your Reference file")
			print(len(duplicated_headers), " duplicated headers out of all headers ", len(REFdict.keys()))
			sys.exit(1)

		duplicated_sequences = basictools.duplicated_seqs(REFdict.values())

		if duplicated_sequences != []:
			print("Warning !", len(duplicated_sequences), " duplicated sequences out of all sequences ", len(REFdict.values()))

	for header, seq in REFdict.items():

		## check that header has right format   transcr | gene
		if len(header.split("|")) > 1:
			ID_tr_REFdict[header.split("|")[0][1:]] = {header:REFdict[header]}
			ID_gene_REFdict[header.split("|")[1]] = {header:REFdict[header]}
		else:
			print("The provided header", header," in reference fasta file do not follow the described format")

		trREFids.append(header.split("|")[0][1:])
		geneREFids.append(header.split("|")[1])

	infile = open(infilename, "r")
	lines = infile.readlines()

	for line in lines:
		line = line.strip()
		ID = line.split("\t")[0]
		pos = line.split("\t")[1]

		if ID in trREFids:
			header = list(ID_tr_REFdict[ID].keys())[0]

			if header not in validTSVids:
				validTSVids[header] = [[int(pos),int(pos)]]
			else:
				validTSVids[header].append([int(pos),int(pos)])
			ID_REFdict[ID] = ID_tr_REFdict[ID]

		elif ID in geneREFids:
			header = list(ID_gene_REFdict[ID].keys())[0]

			if header not in validTSVids:
				validTSVids[header] = [[int(pos),int(pos)]]
			else:
				validTSVids[header].append([int(pos),int(pos)])
			ID_REFdict[ID] = ID_gene_REFdict[ID]
		else:
			nonvalidTSVids.append(ID)

	validcount = len(validTSVids.keys())

	if len(nonvalidTSVids) > 0:
		print("! Warning", str(len(nonvalidTSVids)), " given IDs do not match between your Input tsv file ")

	return validTSVids, ID_REFdict


def find_all_Regex_matches(pattern, string):

	pat = re.compile(pattern)
	pos = 0
	out2 = []

	while (match := pat.search(string, pos)) is not None:
		pos = match.start() + 1
		out2.append(match)
	return out2

def find_aberrant_positions(ID, mainseq, args ):

	position_string_type = args.position_string_type
	position_string = args.position_string

	ab_start_poss = []
	ab_end_poss = []

	finalposs = []

	metapos = {}

	wtlength = len(basictools.returnAAseq(args.mRNA_codon_aminoacid, mainseq)) ## ab pos can not be after CDS

	if args.mRNA_codon_aminoacid == "mRNA":
		wtlength = wtlength * 3

	if position_string_type == "codon":
		ab_start_poss = [i for i, x in enumerate(mainseq) if x == position_string]
		ab_end_poss = ab_start_poss

	elif position_string_type == "aminoacid":
		for ab_pos in re.finditer(re.escape(position_string), mainseq):
			if args.force_frame:
				if basictools.check_inframe(ab_pos.start()):
					ab_start_poss.append(ab_pos.start())
					ab_end_poss.append(ab_pos.start())
			else:
				ab_start_poss.append(ab_pos.start())
				ab_end_poss.append(ab_pos.start())

	elif position_string_type  == "regex":

		for ab_pos in find_all_Regex_matches(position_string, mainseq):
			if args.force_frame:
				if basictools.check_inframe(ab_pos.start()):
					ab_start_poss.append(ab_pos.start())
					ab_end_poss.append(ab_pos.end())
			else:
				ab_start_poss.append(ab_pos.start())
				ab_end_poss.append(ab_pos.end())

	elif position_string_type  == "sequence":

		for ab_pos in re.finditer(re.escape(position_string), mainseq):
			if args.force_frame:
				if basictools.check_inframe(ab_pos.start()):
					ab_start_poss.append(ab_pos.start())
					ab_end_poss.append(ab_pos.end())
			else:
				print(ab_pos)
				ab_start_poss.append(ab_pos.start())
				ab_end_poss.append(ab_pos.end())

	else: ##  position_string_type is rigid position based

		if position_string_type == "startATGpos": # start codon based
			ab_start_poss = [int(args.position_string)]
			ab_end_poss = ab_start_poss

		if position_string_type == "STOPpos": # stop codon based
			ab_start_poss = [wtlength - int(args.position_string)]
			ab_end_poss = ab_start_poss

	if args.NRafter != None or args.NRbefore != None: ### if nrs around are given calculate spacearound
		for i in range(len(ab_start_poss)):

			if spacearound(ab_start_poss[i],ab_end_poss[i],args.NRbefore,args.NRafter,len(mainseq)):
				if ab_end_poss[i] < wtlength:
					finalposs.append([ab_start_poss[i],ab_end_poss[i]])
					metapos[ab_start_poss[i]] = True
				else:
					metapos[ab_start_poss[i]] = False
			else:
				metapos[ab_start_poss[i]] = False

	else:
		for i in range(len(ab_start_poss)):
			if ab_end_poss[i] < wtlength:
				finalposs.append([ab_start_poss[i],ab_end_poss[i]])
				metapos[ab_start_poss[i]] = True

	return finalposs, metapos#, regexheaderseq


###########################

def spacearound(potentialposstart,potentialposend, NRbefore, NRafter, length):

	if (int(potentialposstart) > NRbefore) and (int(potentialposend) < (length - NRafter)):
		return True
	else:
		return False

############

def aberrant_header_names(header, abpos, mainseq, newseq, args):

	########## include insertion fs direction
	########## nil is not taken into account. header shows position of the event whether included or not
	########## all events in same sequence ....

	position = abpos[0]

	if (args.position_string_type == "sequence") and (len(args.position_string) > 10):
			positionstr = "pos" + args.position_string[:3] + "." + str(len(args.position_string)-6) + "." + args.position_string[len(args.position_string)-3:] + "_" + args.mRNA_codon_aminoacid

	elif args.position_string_type == "regex":
		positionstr = "pos" + str(mainseq[abpos[0]:abpos[1]]) + "_" + str(position) + "_" + args.mRNA_codon_aminoacid

	elif args.position_string_type == "FILEpos":
		positionstr = "FILEpos_" + str(position) + "_" + args.mRNA_codon_aminoacid

	elif "pos" not in args.position_string_type:
		positionstr = "pos" + args.position_string + "_" + str(position) + "_" + args.mRNA_codon_aminoacid
	else:
		positionstr = args.position_string + "_" + args.mRNA_codon_aminoacid

	if hasattr(args, "aberrantsequence") and args.aberrantsequence is not None:
		if len(args.aberrantsequence) > 10:
			aberstr = args.aberrantsequence[:3] + "." + str(len(args.aberrantsequence)-6) + "." + args.aberrantsequence[len(args.aberrantsequence)-3:]
		else:
			aberstr = args.aberrantsequence
##############

	newheaders = {}

	if args.codon_aware: ## if you wanted codon awareness it adds the codon to the header
		codon = header.split("|")[-1]
		header = "|".join(header.split("|")[:-1])

		positionstr = "".join(["pos",args.position_string,"codon",codon,"_",str(position),"_",args.mRNA_codon_aminoacid])

	if args.translation_error == "frameshift":
		fsrelabels = {"+1":"p1","-1":"m1"}

		for fsdir, fsseq in newseq.items():

			translation_name = args.translation_error

			if args.mRNA_codon_aminoacid == "mRNA":
				if len(basictools.returnAAseq(args.mRNA_codon_aminoacid, fsseq))*3 == int(position):
					translation_name = translation_name + "Truncation"
			else:
				if (len(basictools.returnAAseq(args.mRNA_codon_aminoacid, fsseq)) == int(position)):
					translation_name = translation_name + "Truncation"

			fslabel = fsdir.split(":")[0]

			headindex = "_".join([translation_name,fslabel,args.position_string_type,positionstr])
			newheader = "|".join([header,headindex])
			newheaders[newheader] = newseq[fsdir]
			if not args.position_string_type == "FILEpos":
				headindex = "_".join([translation_name,fslabel,args.position_string_type,positionstr])
			else:
				posstr = args.position_string_type + str(abpos[0])
				headindex = "_".join([translation_name,fslabel,args.position_string_type,args.position_string_type, str(abpos[0]), args.mRNA_codon_aminoacid ])

			newheader = "|".join([header,headindex])
			newheaders[newheader] = newseq[fsdir]

	if args.translation_error == "repeat":

		repeatindex = args.translation_error + "_" + args.aberrantsequence + "x" + str(args.number)
		headindex = "_".join([repeatindex,args.position_string_type,positionstr])
		# newhead >NP_001005484.2|OR4F5|NC_000001.11|repeat_TGGx15_sequence_posTGG_30_mRNA

	if (args.translation_error == "skip_deletion") or (args.translation_error == "reversion"):
		endpos = int(abpos[0]) + int(args.length)
		headindex = "_".join([args.translation_error,args.length,args.position_string_type,positionstr])
		# newhead >NP_001005484.2|OR4F5|NC_000001.11|skip_deletion_15_sequence_posTGG_30_mRNA

	if (args.translation_error == "truncation") or (args.translation_error == "alternative_start") :
		headindex = "_".join([args.translation_error, args.position_string_type,positionstr])
		# newhead >NP_001005484.2|OR4F5|NC_000001.11|truncation_sequence_posTGG_30_mRNA
		# newhead >NP_001005484.2|OR4F5|NC_000001.11|alternative_start_sequence_posTGG_30_mRNA

	if args.translation_error == "insertion":
		if (args.mRNA_codon_aminoacid != "aminoacid") and (len(args.aberrantsequence) % 3 != 0):
			positionstr = positionstr.replace("_","fs_")

		headindex = "_".join([args.translation_error, aberstr,args.position_string_type,positionstr])

	if args.translation_error != "frameshift": # not frameshift
		newheaders["|".join([header,headindex])] = newseq

	if args.remove_truncation:
		TRremovs = []
		for key in newheaders.keys():
			if "frameshiftTruncation" in key:
				TRremovs.append(key)
		for keyname in TRremovs:
			del newheaders[keyname]

	# FILEpos headers now include the numeric position directly, e.g.
	# frameshift_+1_FILEpos_FILEpos_123_mRNA. Do not delete FILEpos headers here.
	return newheaders


#################################

def aberrant_translation_file(temp, temp2, args, infastadict, posdict):

	if posdict != {}: # the positions were already given, but the fastafile is not in the right format yet
		fastadict = {}

		for ID, subdict in infastadict.items():
			seq = list(infastadict[ID].values())[0]
			WTheader = list(infastadict[ID].keys())[0]
			fastadict[WTheader] = seq
	else:
		fastadict = infastadict

	for header, rnaseq in fastadict.items():

		wt_already = False

		codonseq = basictools.dna_to_codon(rnaseq)
		codon2AA = basictools.codon_AA_dictionary()
		AAseq = basictools.codon_to_AA(codon2AA, codonseq, False)[0]
		mRNA_codon_aminoacid_dict = {"mRNA":rnaseq, "codon":codonseq, "aminoacid":AAseq}
		mainseq = mRNA_codon_aminoacid_dict[args.mRNA_codon_aminoacid]
		originallength = len(mainseq)

		ID = header.split("|")[0][1:]
		Gene = header.split("|")[1]

		if posdict == {}:
			ab_poss, metapos = find_aberrant_positions(ID,mainseq, args) # positions to be found based on input
		else:

			ab_poss = posdict[header] # positions given with input file
			metapos = {}

			for abpos in ab_poss:
				metapos[abpos[0]] = True

		if not args.aminoacid_outlevel:
			outlevel = ""
		else:
			outlevel = args.mRNA_codon_aminoacid

		if not wt_already:
			WTdict = {"|".join([header,"WT"]):mainseq}

			wt_already = True
			WTout = WTdict

			if len(WTout.keys()) > 1:
				basictools.write_fasta_to_temp(temp, WTout, outlevel)
			else:
				basictools.write_fasta_to_temp(temp, WTdict, outlevel)

		for abpos in ab_poss: # for each pair (startpos,endpos)
			if int(abpos[0]) > 1: # and int(abpos[1] < len(codonseq)): # pos has to be between start and stop

				if args.codon_aware: ## if you wanted codon awareness. get codonseq abpos position and at it to the header
					ca_codonpos = codonseq[abpos[0]]
					cheader = header + "|" + str(ca_codonpos)
					newseq = exec_function(args.translation_error, mainseq, abpos, args, codonseq) # for FS is a dict {+1:seq, -1:seq} # for all others its the new sequence 
					aberrantdict = aberrant_header_names(cheader, abpos, mainseq, newseq, args)
				else:
					newseq = exec_function(args.translation_error, mainseq, abpos, args, codonseq)
					aberrantdict = aberrant_header_names(header, abpos, mainseq, newseq, args)

				basictools.write_fasta_to_temp(temp, aberrantdict, outlevel)
				basictools.write_meta_to_temp(temp2, "|".join([header,"WT"]), metapos[abpos[0]], abpos[0], aberrantdict)

##########################################

def aberrant_translation_seq(args, temp):

	if args.mRNA_codon_aminoacid == "aminoacid":
		mainseq = args.input
	else:
		rnaseq = args.input
		codonseq = basictools.dna_to_codon(rnaseq)
		codon2AA = basictools.codon_AA_dictionary()
		AAseq = basictools.codon_to_AA(codon2AA, codonseq, False)[0]

		if args.mRNA_codon_aminoacid == "mRNA":
			mainseq = rnaseq
		else:
			mainseq = codonseq

	header = ">peptideseq_len" + str(len(mainseq))
	WTdict = {"|".join([header,"WT"]):mainseq}

	if not args.aminoacid_outlevel:
		outlevel = ""
	else:
		outlevel = args.mRNA_codon_aminoacid

	WTout = prepwriteout(args, temp, WTdict)

	if len(WTout.keys()) > 1:
		basictools.write_fasta_to_temp(temp, WTout, outlevel)
	else:
		basictools.write_fasta_to_temp(temp, WTdict, outlevel)

	ab_poss, metapos = find_aberrant_positions("",mainseq, args)

	for abpos in ab_poss: # for each pair (startpos,endpos)

		if args.codon_aware:
			ca_codonpos = codonseq[abpos[0]]
			cheader = header + "|" + str(ca_codonpos)
			newseq = exec_function(args.translation_error, mainseq, abpos, args, codonseq) # for FS is a dict {+1:seq, -1:seq}
			aberrantdict = aberrant_header_names(cheader, abpos, mainseq, newseq, args, regexheaderseq)
		else:
			newseq = exec_function(args.translation_error, mainseq, abpos, args, codonseq)
			aberrantdict = aberrant_header_names(header, abpos, mainseq, newseq, args, regexheaderseq)

		basictools.write_fasta_to_temp(temp, aberrantdict, outlevel)

##########################################

def filter_AB_seqlength(args, maindict):

	keywords = ["WT","Ctrl"]
	outdict = {}

	for header, AAseq in maindict.items():
		if "WT" not in header.split("|")[3] and "Ctrl" not in header.split("|")[3]:
			pos = int(header.split("|")[-1].split("_")[4])

			if args.mRNA_codon_aminoacid == "mRNA":
				pos = math.ceil(pos / 3)

			if len(AAseq[pos:]) > args.sequencelength:
				outdict[header] = AAseq
		else:
			outdict[header] = AAseq

	return outdict


def select_n_ctrls(controlnumber,maindict):

	iterations = 1
	n_controldict = {}
	already = []

	while (iterations < controlnumber+1):

		roundstart  = (len(n_controldict.keys()))

		key = random.choices(list(maindict.keys()),k=1)[0]

		if key.split("|")[-1] != "WT" and "Ctrl_" not in key and key not in already:
			already.append(key)
			mainseq = maindict[key]

			n_controldict[key] = mainseq

			iterations = iterations + 1

	return n_controldict

def scramble_controls(n_controldict):

	controldict = {}

	for abheader, seq in n_controldict.items():
		abpos = int(abheader.split("|")[-1].split("_")[-2])
		mca = abheader.split("|")[-1].split("_")[-1]

		if mca == "mRNA":
			abpos = abpos*3

		abseq = seq[abpos:]
		scrambledseq = list(abseq)
		random.shuffle(scrambledseq)
		scrambledseq  = "".join(scrambledseq)
		ctrlscrambledseq = seq[:abpos] + scrambledseq

		ctrlheader = abheader.split("|")[:-1]
		ctrlheader = "|".join(ctrlheader)  + "|" + "ScrambledCtrl_" + "_".join(abheader.split("|")[-1].split("_")[1:])

		controldict[ctrlheader] = ctrlscrambledseq

	return controldict


def reverse_controls(n_controldict):

	controldict = {}

	for abheader, seq in n_controldict.items():

		abpos = int(abheader.split("|")[-1].split("_")[-2])
		mca = abheader.split("|")[-1].split("_")[-1]

		if mca == "mRNA":
			abpos = abpos*3

		abseq = seq[abpos:]

		revseq = abseq[::-1]
		ctrlscrambledseq = seq[:abpos] + revseq

		ctrlheader = abheader.split("|")[:-1]
		ctrlheader = "|".join(ctrlheader)  + "|" + "ReversedCtrl_" + "_".join(abheader.split("|")[-1].split("_")[1:])

		controldict[ctrlheader] = ctrlscrambledseq

	return controldict




def runcontrols2(args, maindict):

	controldict = {}
	controlnumber = int(args.controlnumber)
	n_controldict = select_n_ctrls(controlnumber,maindict)

	if args.controltype == "all":#
		controldict.update(reverse_controls(n_controldict))
		controldict.update(scramble_controls(n_controldict))

	elif args.controltype == "reverse":
		controldict.update(reverse_controls(n_controldict))

	elif args.controltype == "scramble":
		controldict.update(scramble_controls(n_controldict))

	return controldict

######
def second_iteration(premaindictcopy, args):

	iterate_dict = {}
	WTonly = []

	for header, seq in premaindictcopy.items():
		baseheader = "|".join(header.split("|")[:-1])
		condition = header.split("|")[-1].split("_")[0]

		if baseheader not in iterate_dict:
			iterate_dict[baseheader] = {}
			if condition in iterate_dict[baseheader]:
				iterate_dict[baseheader][condition] = iterate_dict[baseheader][condition] + 1 
			else:
				iterate_dict[baseheader][condition] = 1
		else:
			if condition in iterate_dict[baseheader]:
				iterate_dict[baseheader][condition] = iterate_dict[baseheader][condition] + 1 
			else:
				iterate_dict[baseheader][condition] = 1
        
	for baseheader in iterate_dict.keys():
		if len(iterate_dict[baseheader]) == 1: 
			WTname = baseheader + "|WT"
			WTonly.append(WTname)

	for singlekey in WTonly:
		premaindictcopy.pop(singlekey, None)

	return premaindictcopy

def prepwriteout(args, temp, controldict):

	outdict = {}

	if args.sequencelength not in (False, None):
		premaindict = filter_AB_seqlength(args, basictools.get_fastadict(temp.name))
	else:
		premaindict = basictools.get_fastadict(temp.name)


	if controldict != {}:
		premaindict.update(controldict)
		controldict2 = runcontrols2(args, premaindict)
		premaindict.update(controldict2) ## adds controls to temp fastafile

	if args.translation_error != 'mutation' and not is_vcf_mode(args):

		seconditerationstring =  "Removing WT headers without associated aberrant header"
		logging.info(seconditerationstring)
		maindict = second_iteration(premaindict, args)
	else:
		maindict = premaindict

	if args.trypsin != False:
		trp_outdict = TRYPSIN(maindict)
		outdict = Trypsin_header(trp_outdict)

	if args.trimm != False:
		outdict = trimm_seqs(args, maindict)

	if outdict == {}:
		outdict = maindict

	return outdict

def runcontrols1(temp, args, REFdict):

	controlfastafile = basictools.get_fastadict(temp.name) ##
	controldict = {}
	controlnumber = int(args.controlnumber)

	iterations = 1

	while (iterations < controlnumber+1):
		roundstart  = (len(controldict.keys()))

		key = random.choices(list(controlfastafile.keys()),k=1)[0]

		if key.split("|")[-1] != "WT":

			newkey = key.split("|")[:-1]
			newkey = "|".join(newkey) + "|WT"
		else:
			newkey = key

		mainseq = controlfastafile[newkey]

		if "|".join(key.split("|")[:-1]) in REFdict:
			rnaseq = REFdict["|".join(key.split("|")[:-1])]
		else:
			rnaseq = REFdict[key.split("|")[0][1:]]["|".join(key.split("|")[:-1])]

		codonseq = basictools.dna_to_codon(rnaseq)
		codon2AA = basictools.codon_AA_dictionary()
		AAseq = basictools.codon_to_AA(codon2AA, codonseq, False)[0]

		if args.mRNA_codon_aminoacid == "mRNA":
			mainseq = rnaseq
		else:
			mainseq = codonseq

		randrangestart = 1
		randrangeend = len(codonseq)-2

		if args.NRbefore != None:
			randrangestart = args.NRbefore + 1

		if args.NRafter != None:
			randrangeend = (len(codonseq)-2) -  args.NRbefore

		if randrangestart < randrangeend:
			randABpos = random.randint(randrangestart,randrangeend) ## ensures that event is in CDS
			randABrange = (randABpos, randABpos+1)

		randABdict = exec_function(args.translation_error, mainseq, randABrange, args, codonseq)

		if isinstance(randABdict, dict):
			for head in randABdict.keys():
				if ":" in head:
					chead = head.replace(":","_")
					chead1 = chead.split("_")[0]
					chead2 = chead.split("_")[1]
				else:
					chead = head
					chead1 = chead
					chead2 = "None"

				randABseq = randABdict[head]
				randheader = "|".join(key.split("|")[:-1]) + "|" + str(args.translation_error) + "Ctrl" + "_" + str(chead1) + "_" + str(args.mRNA_codon_aminoacid) + "_posCtrl_" + str(randABpos) + "_" + str(chead2)

				randABAAseq = basictools.returnAAseq(args.mRNA_codon_aminoacid, randABseq) #####

				if args.sequencelength != None:
					if len(randABAAseq[randABpos:]) > args.sequencelength:
						controldict[randheader] = randABAAseq

				elif randABAAseq != "":
					controldict[randheader] = randABAAseq
		else:
 
			randABAAseq = basictools.returnAAseq(args.mRNA_codon_aminoacid, randABdict)
			randABheader = "|".join(key.split("|")[:-1]) + "|" + str(args.translation_error) + "Ctrl" + str(randABpos)

			if args.sequencelength != None:
				if len(randABAAseq[randABpos:]) > args.sequencelength:
					controldict[randABheader] = randABAAseq

			elif randABAAseq != "":
				controldict[randABheader] = randABAAseq

		if len(controldict.keys()) > roundstart:
			iterations = iterations + 1

	if len(controldict.keys()) > int(args.controlnumber):
		allcontrolkeys = list(controldict.keys())
		diff = int(len(allcontrolkeys) - controlnumber)
		remcontrolkeys = allcontrolkeys[:diff]

		for key in remcontrolkeys:
			controldict.pop(key, None)

	return controldict

#####################


def mkdir_meta(metapathname):

	try:
		os.mkdir(metapathname)
		logging.info(f" '{metapathname}' created successfully.")
	except FileExistsError:
		logging.warning(f" '{metapathname}' already exists.Metadata was overwritten.")
	except PermissionError:
		logging.info(f"Permission denied: Unable to create '{metapathname}'.")
		logging.warning("Pipeline can not run without permissions. Please assure permission policy is correct.")
		sys.exit(1)
	except Exception as e:
		logging.info(f"An unknown error occurred: {e}")
		sys.exit(1)


def runmeta(temp, temp2, args):

	controlfastafile = basictools.get_fastadict(temp.name)

	outfilename = args.outfilename.split("/")[-1]

	if "/" in args.outfilename:
		metaoutname = args.outfilename.split(".")[0] + "_meta.txt"
		logfileoutname = args.outfilename.split(".")[0] + "_logfile.txt"
	else:
		cwd = os.getcwd()
		metapath = cwd + "/aberr_fasta/metadata/" #######
		mkdir_meta(metapath)
		metaoutname = metapath + outfilename.split(".")[0] + "_meta.txt"

		logfileoutname = metapath + outfilename.split(".")[0] + "_logfile.txt"

	logfile = open(logfileoutname, "w")

	metaoutfile = open(metaoutname, "w")
	metaoutfile.write("ID\tGene\tlength_WT\tlength_aberrant\tdescription\taberrant_position\tspacepassed\tUTR3\tminimimumlength_passed\tPeptide_created\n")

	infile = open(temp2.name, "r")
	lines = infile.readlines()

	UTR3counts = 0
	minlenpassedcounts = 0
	totalpeptidescounts = 0
	peptideincludedcounts = 0

	for line in lines:
		totalpeptidescounts = totalpeptidescounts + 1 

		UTR3 = False
		minlenpassed = False
		peptideincluded = False

		line = line.strip()
		listline = line.split("\t")

		pos = listline[2]
		spacepassed = listline[3]

		WTheader = listline[0]
		WTseq = controlfastafile[WTheader]

		ABheader = listline[1]
		ABdescr = ABheader.split("|")[3]
		ABseq = controlfastafile[ABheader]

		ID = WTheader.split("|")[0][1:]
		Gene = WTheader.split("|")[1]

		if args.mRNA_codon_aminoacid == "mRNA":
			newpos = math.ceil(int(pos) / 3)
			lenABseq = len(ABseq) * 3
			lenWTseq = len(WTseq) * 3
		else:
			newpos =  int(pos)
			lenABseq = len(ABseq)
			lenWTseq = len(WTseq)

		if len(ABseq) > len(WTseq): ## based on codon sequence (CDS) length
			UTR3 = True
			UTR3counts = UTR3counts + 1

		if len(ABseq[newpos:]) > args.sequencelength:
			minlenpassed = True
			minlenpassedcounts = minlenpassedcounts + 1

		if spacepassed and minlenpassed:
			peptideincluded = True
			peptideincludedcounts = peptideincludedcounts + 1

		outstring = "\t".join([ID,Gene,str(lenWTseq),str(lenABseq), ABdescr, str(pos), str(spacepassed), str(UTR3), str(minlenpassed), str(peptideincluded)])

		metaoutfile.write(outstring)
		metaoutfile.write("\n")



	logfile.write("Arguments:\n")
	logfile.write(json.dumps(vars(args), indent=4))
	logfile.write("\n")

	logfile.write(f"#total_peptides = {totalpeptidescounts}\n")
	logfile.write(f"#peptides_included = {peptideincludedcounts}\n")
	logfile.write(f"#peptides_passing_minimum_length = {minlenpassedcounts}\n")
	logfile.write(f"#peptides_longer_than_WT = {UTR3counts}\n")

	logfile.close()

	metaoutfile.close()

################################



def TRYPSIN(maindict):

	cutAA = ["R","K"]
	outdict = {}

	for header, seq in maindict.items():
		cut_sites = [0]

		if header.split("|")[-1] == "WT":
			for i in range(len(seq)-1):
				if seq[i] in cutAA and seq[i+1] != "P":
					cut_sites.append(i)

			if cut_sites[-1] != len(seq):
				cut_sites.append(len(seq))

			if len(cut_sites) > 2:
				for j in range(0, len(cut_sites) - 1):
					selected = seq[cut_sites[j]: cut_sites[j + 1]]

					if 5 <= len(selected) <= 55:
						outheader = header + "|" + str(cut_sites[j]) + "_" + str(cut_sites[j + 1])
						outdict[outheader] = selected
			else:
				if 5 <= len(seq) <= 55:
					outheader = header + "|fulllength"
					outdict[outheader] = seq

		else: ######### FS ones

			for i in range(len(seq)-1):
				if seq[i] in cutAA and seq[i+1] != "P":
					cut_sites.append(i)

			if cut_sites[-1] != len(seq):
				cut_sites.append(len(seq))

			if len(cut_sites) > 2:
				for j in range(0, len(cut_sites) - 1):
					selected = seq[cut_sites[j]: cut_sites[j + 1]]

					if 5 <= len(selected) <= 55:
						if selected not in outdict.values():
							Abpart = header.split("|")[-1]
							Abpos = int(Abpart.split("_")[4])

							if cut_sites[j] < Abpos < cut_sites[j+1]: 
								outheader = header + "|" + str(cut_sites[j]) + "_" + str(cut_sites[j+1]) + "_chimera"
								outdict[outheader] = selected
							else:
								outheader = header + "|" + str(cut_sites[j]) + "_" + str(cut_sites[j+1])
								outdict[outheader] = selected
			else:
				if 5 <= len(seq) <= 55:
					if selected not in outdict.values():
						outheader = header + "|fulllength_chimera"
						outdict[outheader] = selected

	return outdict


def Trypsin_header(trp_outdict):

	counter = 1
	outdict = {}

	for k,v in trp_outdict.items():
		k2 = k.split("|")

		if not "WT" in k2[-2]:
			ID = k2[0]
			Gene = k2[1]
			cond = k2[3]
			cuts = k2[-1]

			cond2 = cond.split("_")
			cond3 = cond2[0] + "_" + cond2[1] + "_" +  cond2[4] + "_" + cuts
			outs = "_".join([ID,Gene,cond3, str(counter)]) 

			outdict[outs] = v
		else:
			ID = k2[0]
			Gene = k2[1]
			cond = k2[3]
			cuts = k2[-1]

			outs = "_".join([ID,Gene,cond,cuts, str(counter)])

			outdict[outs] = v

		counter = counter + 1

	return outdict

###############################

def main():

	logging.basicConfig(level = logging.INFO , format='{asctime} - {levelname} - {message}',filemode='w', style = "{", datefmt = "[%d-%m-%Y] [%H:%M] ")

	starttime = time.time()

	parser = argparse.ArgumentParser()
	allparsers.parse_all_findaberrant(parser)
	args = parser.parse_args()
	allparsers.validate_argparse_findaberrant_inputs(args)

	temp = tempfile.NamedTemporaryFile(delete=False)
	temp2 = tempfile.NamedTemporaryFile(delete=False)

##############################################

	if os.path.isfile(args.input): 	# input is a file [fasta , tsv, vcf, or mt.txt]

		if args.translation_error == 'mutation':  ### mutations in format gene:mt
			aberrantdict = call_transvar.cdna_mt(args.input, "")
			logging.info("Mutations mapped")
			basictools.write_fasta_to_temp(temp, aberrantdict, "aminoacid")

			args.aminoacid_outlevel = ""
			args.metadata = False
			args.trimm = False

		elif is_vcf_mode(args):
			if args.input.endswith(".vcf"): ### vcf file disable ...
				args.controltype = None				
				args.mRNA_codon_aminoacid = "mRNA"
				args.codon_aware = False
				aberrantdict = call_transvar.vcf_mt(args.input, "")
				logging.info("vcf mutations mapped")
				args.aminoacid_outlevel = ""
				args.metadata = False
				args.trimm = False

				basictools.write_fasta_to_temp(temp, aberrantdict, "aminoacid")

		elif args.position_string_type == "FILEpos":
		#	else: ### input is tsv file with gene\tpos

			if os.path.isfile(args.Reference): ### reference is given
				REFdict = basictools.get_fastadict(args.Reference)
				posdict, fastadict = validate_TSV_REF(REFdict,args.input, "default")

			else: ### reference is default
				REFlevel = {"hg38":"./validated_fasta/gencode_test_noncan.fasta"} ####
				REFfile = REFlevel[args.Reference]
				REFdict = basictools.get_fastadict(REFfile)

				if args.chr: ### chr based has always default reference ... fix in allparse ### chr disable ... 
					preposdict = call_transvar.chr_to_ENSTpos(args.input, "ganno")
					posdict, fastadict = call_transvar.validate_chrpos(preposdict, REFdict)
				else:
					posdict, fastadict = validate_TSV_REF(REFdict,args.input, "given")

			logging.info("Fasta dict created")

		else: # input is not a file pos, so positions should be found.
			fastadict = basictools.get_fastadict(args.input) # fastafile of all transcripts in dictionary { ENST : RNAseq }
			posdict = {}
			logging.info("Fasta dict created")
#####################
		if str(args.translation_error) != 'mutation' and not is_vcf_mode(args): 
			aberrant_translation_file(temp, temp2, args, fastadict, posdict)

	else: 				# input is a single sequence
		aberrant_translation_seq(args, temp)

########################################################################################

	aberranttime = time.time() #

	elapsed_time = aberranttime - starttime #
	logging.info("Aberrant sequences finished")

	temp.seek(0)
	temp2.seek(0)

	if args.metadata not in (False, None): ### if you want metadata
		runmeta(temp, temp2, args) # creates metadata file and writes it
		logging.info("Metadata file finished")

	if str(args.translation_error) != 'mutation' and not is_vcf_mode(args):
		if args.controltype != None: ### if you want controls
			controldict = runcontrols1(temp, args, fastadict)
			outdict = prepwriteout(args, temp, controldict) ### 
		else:
			outdict = prepwriteout(args, temp, {})
	else:
		outdict = prepwriteout(args, temp, {})

############################

	if "/" in args.outfilename:
		 outfilename = args.outfilename
	else:
		cwd = os.getcwd()
		outfilename = cwd + "/aberr_fasta/" + args.outfilename.split("/")[-1] #######

	infotext = "Results written to " + outfilename
	logging.info(infotext)

	basictools.write_fasta_to_file(outfilename, outdict, args.aminoacid_outlevel)

	########################################################################

	endtime = time.time() #

	endtime = endtime - starttime #
	endtime = round((endtime/60),2)
	xstring =  "Pipe finished in " + str(endtime) + " m"
	logging.info(xstring )

	temp.close()
	temp2.close()

main()
