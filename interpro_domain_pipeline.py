#!/usr/bin/env python3

import csv
import logging
import sqlite3
import time
import logging

from collections import defaultdict

from lib import allparsers


#############################################3

def mappos(abpos, start, end):

	abpos = int(abpos)
	start = int(start)
	end = int(end)

	if start < abpos < end:
		return "spanning"
	elif abpos > end:
		return "before"
	elif abpos <= start:
		return "after"
	elif abpos == end:
		return "truncated"

	return "NA"


def get_length_bin(length):
	length = int(length)

	if length < 100:
		return "0-100"
	elif length < 200:
		return "100-200"
	elif length < 300:
		return "200-300"
	elif length < 500:
		return "300-500"
	elif length < 750:
		return "500-750"
	elif length < 1000:
		return "750-1000"
	elif length < 2000:
		return "1000-2000"
	else:
		return "2000+"


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

		if "|" in descr:
			FSpos = descr.split("|")[0]
		else:
			FSpos = header.split("|")[-1].split("_")[-2]

		mca = header.split("|")[-1].split("_")[-1]
		if mca == "mRNA":
			FSpos = int(FSpos)
			FSpos = str(int(FSpos/ 3))

		return {"header":header, "ID":ID, "Gene":Gene, "FSdir": FSdir, "FSpos": FSpos, "mca": mca, "WT": False, "condition":condition}



def is_wt_id(protein_id):
	return protein_id.split("|")[-1] == "WT"


def min_evalue(values):
	valid = []

	for value in values:
		if value in ["NA", "-", "", None]:
			continue
		try:
			valid.append(float(value))
		except ValueError:
			continue

	if len(valid) == 0:
		return "NA"

	return min(valid)



# FASTA parsing
#########################################################

def parse_fasta(fasta):

	proteins = {}
	current_id = None
	current_seq = []

	def store_current():
		if current_id is None:
			return

		fields = current_id.split("|")
		enst = fields[0]
 
		if len(fields) > 1:
			gene = fields[1]
		else:
			gene = "NA" 

		is_wt = is_wt_id(current_id) ### condition

		if not is_wt:
			headerdict = header2dict(current_id, "")
			fspos = headerdict["FSpos"]

		length = len("".join(current_seq))
		proteins[current_id] = {"enst": enst, "gene": gene, "isWT": is_wt, "fspos": fspos, "length": length, "lengthbin": get_length_bin(length) if length > 0 else "NA"}

	with open(fasta, "r") as infile:
		for line in infile:
			line = line.strip()
			if line == "":
				continue

			if line.startswith(">"):
				store_current()
				current_id = line[1:]
				current_seq = []
			else:
				current_seq.append(line)

	store_current()

	return proteins



# InterProScan parsing
#########################################################

def parse_raw_interpro(raw_tsv, proteins=None, extracted_output=None):


	extracted_rows = []
	skipped = 0

	with open(raw_tsv, "r") as infile:
		for line in infile:
			line = line.rstrip("\n")
			if line == "":
				continue

			row = line.split("\t")

			# Skip already-normalized headers or malformed rows.
			if row[0] in ["ID", "ENSTID"]:
				 continue
			if len(row) < 9:
				skipped += 1
				continue

			protein_id = row[0]

			try:
				length = int(row[2])
				database = row[3]
				accession = row[4]
				featurename = row[5]
				start = int(row[6])
				end = int(row[7])
			except (ValueError, IndexError):
				skipped += 1
				continue

			evalue = row[8]
			if evalue == "-":
				evalue = "NA"

			fields = protein_id.split("|")
			enst = fields[0]
			gene = fields[1] if len(fields) > 1 else "NA"
			featureid = database + "|" + accession
			featurelength = end - start + 1
			lengthbin = get_length_bin(length)

			if is_wt_id(protein_id):
				region = "WT"
			else:
				fspos = None
				if proteins is not None and protein_id in proteins:
					fspos = proteins[protein_id]["fspos"]
				if fspos is None:
					fspos = parse_fspos_from_header(protein_id) ######parse_fspos_from_header header2dict(header, cond)
				region = mappos(fspos, start, end) if fspos is not None else "NA"

			feature = { "featureid": featureid "database": database, "accession": accession, "featurename": featurename, "start": start, "end": end, "featurelength": featurelength, "evalue": evalue, "region": region, "length": str(length), "lengthbin": lengthbin, "enst": enst, "gene": gene }
			features[protein_id].append(feature)

			if extracted_output is not None:
				extracted_rows.append([ enst, protein_id, gene, str(length), database, accession, str(start), str(end), featurename, str(featurelength), evalue, featureid, region, lengthbin])

	if extracted_output is not None:
		write_extracted_features(extracted_output, extracted_rows)

	if skipped > 0:
		logging.warning("Skipped %s malformed InterProScan rows", skipped)

	return features


def write_extracted_features(outfile, rows):
	header = ["ENSTID", "ID", "gene", "Length", "Database", "Accession", "Start", "End", "FeatureName", "FeatureLength", "Evalue", "FeatureID", "Region", "lengthbin"]

	with open(outfile, "w") as out:
		out.write("\t".join(header) + "\n")
		for row in sorted(rows, key=lambda x: (x[0], x[1], x[11], int(x[6]))):
			out.write("\t".join(row) + "\n")



# Domain count matrix
#########################################################

def build_enst_mapping(proteins):
	mapping = defaultdict(list)

	for protein_id, protein in proteins.items():
		mapping[protein["enst"]].append(protein_id)

	return mapping


def write_domain_count_matrix(proteins, features, outfile):

	enst_mapping = build_enst_mapping(proteins)
	header = [ "ENSTID", "WT_ID", "OOF_ID", "Gene","OOF_Length", "OOF_LengthBin", "FeatureID", "Database", "FeatureName", "Region", "WT_Count", "OOF_Count", "Difference", "AbsDifference", "WT_MinEvalue", "OOF_MinEvalue"]

	n_rows = 0

	with open(outfile, "w") as out:
		out.write("\t".join(header) + "\n")

		for enst in sorted(enst_mapping):
			ids = enst_mapping[enst]
			wt_ids = [x for x in ids if proteins[x]["isWT"]]

			if len(wt_ids) == 0:
				continue

			wt_id = wt_ids[0]
			oof_ids = [x for x in ids if not proteins[x]["isWT"]]
			wt_features = features.get(wt_id, [])

			for oof_id in oof_ids:
				fspos = proteins[oof_id]["fspos"]

				if fspos is None:
					continue

				gene = proteins[oof_id]["gene"]
				oof_features = features.get(oof_id, [])
				oof_length = proteins[oof_id]["length"]
				oof_lengthbin = proteins[oof_id]["lengthbin"]

				if len(oof_features) > 0:
					oof_length = oof_features[0]["length"]
					oof_lengthbin = oof_features[0]["lengthbin"]

				wt_counts = defaultdict(lambda: defaultdict(int))
				wt_evalues = defaultdict(list)
				feature_info = {}

				for feat in wt_features:
					region = mappos(fspos, feat["start"], feat["end"])
					featureid = feat["featureid"]

					wt_counts[featureid][region] += 1
					wt_evalues[featureid].append(feat["evalue"])
					feature_info[featureid] = (feat["database"], feat["featurename"])

				oof_counts = defaultdict(lambda: defaultdict(int))
				oof_evalues = defaultdict(list)

				for feat in oof_features:
					featureid = feat["featureid"]
					region = feat["region"]

					oof_counts[featureid][region] += 1
					oof_evalues[featureid].append(feat["evalue"])
					feature_info[featureid] = (feat["database"], feat["featurename"])

				all_featureids = set(wt_counts.keys())
				all_featureids.update(oof_counts.keys())

				for featureid in sorted(all_featureids):
					database, featurename = feature_info[featureid]

					regions = set(wt_counts[featureid].keys())
					regions.update(oof_counts[featureid].keys())

					for region in sorted(regions):
						wt_count = wt_counts[featureid].get(region, 0)
						oof_count = oof_counts[featureid].get(region, 0)
						difference = oof_count - wt_count
						absdifference = abs(difference)

						out.write("\t".join([ enst, wt_id, oof_id, gene, str(oof_length), str(oof_lengthbin), featureid, database, featurename, region, str(wt_count), str(oof_count), str(difference),str(absdifference), str(min_evalue(wt_evalues[featureid])), str(min_evalue(oof_evalues[featureid]))]) + "\n")

						n_rows += 1

	return n_rows



# Background statistics
#########################################################

def load_background(dbfile):
	conn = sqlite3.connect(dbfile)
	cur = conn.cursor()

	stats_exact = {}
	stats_global = {}

	cur.execute(""" SELECT lengthbin, featureid, mean_copies, sd_copies, total_copies, n_proteins FROM feature_copy_statistics """)

	global_tmp = defaultdict(lambda: {"copies": 0, "proteins": 0, "sd_sum": 0.0, "sd_n": 0})

	for lengthbin, featureid, mean_copies, sd_copies, total_copies, n_proteins in cur.fetchall():
		stats_exact[(featureid, lengthbin)] = (
			mean_copies,
			sd_copies,
			"LengthMatched",
			True
		)

		if total_copies is not None:
			global_tmp[featureid]["copies"] += total_copies

		if n_proteins is not None:
			global_tmp[featureid]["proteins"] += n_proteins

		if sd_copies is not None and n_proteins is not None:
			global_tmp[featureid]["sd_sum"] += float(sd_copies) * int(n_proteins)
			global_tmp[featureid]["sd_n"] += int(n_proteins)

	for featureid, d in global_tmp.items():
		if d["proteins"] > 0:
			mean = float(d["copies"]) / float(d["proteins"])
		else:
			mean = 0.0

		if d["sd_n"] > 0:
			sd = d["sd_sum"] / d["sd_n"]
		else:
			sd = 0.0

		stats_global[featureid] = ( mean, sd, "FeatureWide", True)

	dist_exact = defaultdict(list)
	dist_global = defaultdict(lambda: defaultdict(int))

	cur.execute(""" SELECT lengthbin, featureid, copy_number, n_proteins FROM feature_copy_distribution """)

	for lengthbin, featureid, copy_number, n_proteins in cur.fetchall():
		dist_exact[(featureid, lengthbin)].append((int(copy_number), int(n_proteins)))
		dist_global[featureid][int(copy_number)] += int(n_proteins)

	conn.close()

	return stats_exact, stats_global, dist_exact, dist_global


def get_statistics(featureid, lengthbin, stats_exact, stats_global):
	key = (featureid, lengthbin)

	if key in stats_exact:
		return stats_exact[key]

	if featureid in stats_global:
		return stats_global[featureid]

	return (0.0, 0.0, "Absent", False)


def get_empirical_pvalue(featureid, lengthbin, observed, dist_exact, dist_global):
	key = (featureid, lengthbin)

	if key in dist_exact:
		rows = dist_exact[key]
	elif featureid in dist_global:
		rows = list(dist_global[featureid].items())
	else:
		return 0.0

	total = 0
	ge = 0

	for copy_number, n_proteins in rows:
		total += n_proteins
		if copy_number >= observed:
			ge += n_proteins

	return (ge + 1) / (total + 1)


def zscore(observed, mean, sd):
	if sd is None or float(sd) == 0:
		return "NA"

	return (float(observed) - float(mean)) / float(sd)


def add_background_statistics(count_matrix, dbfile, outfile, absent_mode="na"):

	stats_exact, stats_global, dist_exact, dist_global = load_background(dbfile)
	totals = {}

	with open(count_matrix) as f:
		reader = csv.DictReader(f, delimiter="\t")

		for row in reader:
			key = (row["ENSTID"], row["WT_ID"], row["OOF_ID"], row["FeatureID"])

			if key not in totals:
				totals[key] = [0, 0]

			totals[key][0] += int(row["WT_Count"])
			totals[key][1] += int(row["OOF_Count"])

	n_rows = 0

	with open(count_matrix) as f, open(outfile, "w") as out:
		reader = csv.DictReader(f, delimiter="\t")

		header = reader.fieldnames + ["WT_TotalCopies", "OOF_TotalCopies", "BackgroundMean", "BackgroundSD", "BackgroundSource", "FeatureInBackground", "Zscore", "EmpiricalP"]

		out.write("\t".join(header) + "\n")

		for row in reader:
			key = (row["ENSTID"], row["WT_ID"], row["OOF_ID"], row["FeatureID"])

			wt_total, oof_total = totals[key]

			mean, sd, source, present = get_statistics(row["FeatureID"], row["OOF_LengthBin"], stats_exact, stats_global)

			p = get_empirical_pvalue( row["FeatureID"], row["OOF_LengthBin"], oof_total, dist_exact, dist_global)

			if source == "Absent":
				if absent_mode == "high" and oof_total > 0:
					p = 0
					z = "Inf"
				elif absent_mode == "zero":
					p = 0
					z = "NA"
				else:
					p = "NA"
					z = "NA"
			else:
				z = zscore(oof_total, mean, sd)

			vals = [row[c] for c in reader.fieldnames]
			vals.extend([ str(wt_total), str(oof_total),str(mean),str(sd), str(source), str(present), str(z), str(p)])

			out.write("\t".join(vals) + "\n")
			n_rows += 1

	return n_rows



#########################################################

def main():


	logging.basicConfig(level = logging.INFO , format='{asctime} - {levelname} - {message}',filemode='w', style = "{", datefmt = "[%d-%m-%Y] [%H:%M] ")
	starttime = time.time()
	logging.info("start run")
	logging.info("reading FASTA")

	parser=argparse.ArgumentParser()
	allparsers.parser_domain_pipeline(parser)
	args = parser.parse_args()
	proteins = parse_fasta(args.fasta)

	logging.info("loaded %s proteins from FASTA", len(proteins))
	logging.info("parsing InterProScan TSV")

	features = parse_raw_interpro(args.input, proteins=proteins, extracted_output=args.extracted_output) 

	logging.info("loaded features for %s proteins", len(features))
	logging.info("writing raw domain count matrix")

	count_rows = write_domain_count_matrix(proteins, features, args.count_matrix )
	logging.info("wrote %s rows to %s", count_rows, args.count_matrix)
	logging.info("adding background statistics")

	stat_rows = add_background_statistics(args.count_matrix, args.database, args.output )#, absent_mode=args.absent_background_mode)
	logging.info("wrote %s rows to %s", stat_rows, args.output)

	endtime = round((time.time() - starttime) / 60, 2)
	logging.info("run finished in %s m", endtime)


if __name__ == "__main__":
	main()
