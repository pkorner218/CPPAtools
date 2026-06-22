#!/usr/bin/env python3

import sys
import csv
import sqlite3
import statistics
from collections import defaultdict, Counter


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


def create_database(outfile):

	conn = sqlite3.connect(outfile)
	cur = conn.cursor()
	cur.execute(""" CREATE TABLE feature_annotation (featureid TEXT PRIMARY KEY, database_name TEXT, accession TEXT, feature_name TEXT ) """)
	cur.execute(""" CREATE TABLE feature_copy_distribution ( lengthbin TEXT, featureid TEXT, copy_number INTEGER, n_proteins INTEGER ) """)
	cur.execute(""" CREATE TABLE feature_copy_statistics ( lengthbin TEXT, featureid TEXT, mean_copies REAL, median_copies REAL, sd_copies REAL, total_copies INTEGER, n_proteins INTEGER) """)
	cur.execute(""" CREATE TABLE protein_feature_distribution ( lengthbin TEXT, feature_count INTEGER, n_proteins INTEGER ) """)
	cur.execute(""" CREATE TABLE length_statistics ( lengthbin TEXT PRIMARY KEY, n_sequences INTEGER, mean_features REAL, median_features REAL, sd_features REAL ) """)
	conn.commit()

	return conn, cur


def build_background(infile, cur):

	protein_lengthbin = {}
	protein_feature_count = defaultdict(int)
	protein_feature_copies = defaultdict(lambda: defaultdict(int))

	annotations = {}

	with open(infile, "r") as handle:

		reader = csv.reader(handle, delimiter="\t")

		for row in reader:
			if len(row) < 9:
				continue

			seqid = row[0]

			try:
				length = int(row[2])
			except:
				continue

			database = row[3]
			accession = row[4]
			feature_name = row[5]
			featureid = database + "|" + accession

			lengthbin = get_length_bin(length)
			protein_lengthbin[seqid] = lengthbin

			protein_feature_copies[seqid][featureid] += 1
			protein_feature_count[seqid] += 1

			annotations[featureid] = (database, accession, feature_name)


	for featureid, values in annotations.items():
		cur.execute(""" INSERT OR REPLACE INTO feature_annotation VALUES (?,?,?,?) """, (featureid, values[0], values[1], values[2]))

	proteins_per_bin = defaultdict(set)

	for seqid, lengthbin in protein_lengthbin.items():
		proteins_per_bin[lengthbin].add(seqid)

	counts_per_bin = defaultdict(list)

	for seqid, count in protein_feature_count.items():
		lengthbin = protein_lengthbin[seqid]
		counts_per_bin[lengthbin].append(count)

	for lengthbin, values in counts_per_bin.items():
		distribution = Counter(values)

		for feature_count, n_proteins in distribution.items():
			cur.execute(""" INSERT INTO protein_feature_distribution VALUES (?,?,?) """, (lengthbin, feature_count, n_proteins))

		mean_features = statistics.mean(values)
		median_features = statistics.median(values)

		if len(values) > 1:
			sd_features = statistics.stdev(values)
		else:
			sd_features = 0

		cur.execute(""" INSERT INTO length_statistics VALUES (?,?,?,?,?) """, (lengthbin, len(values), mean_features,median_features, sd_features))

	feature_counts_by_bin = defaultdict(lambda: defaultdict(list))

	for seqid, featuredict in protein_feature_copies.items():
		lengthbin = protein_lengthbin[seqid]

		for featureid, copies in featuredict.items():
			feature_counts_by_bin[lengthbin][featureid].append(copies)

	for lengthbin in feature_counts_by_bin:
		total_proteins = len(proteins_per_bin[lengthbin])

		for featureid in feature_counts_by_bin[lengthbin]:
			observed = feature_counts_by_bin[lengthbin][featureid]

			n_nonzero = len(observed)
			n_zero = total_proteins - n_nonzero
			distribution = Counter(observed)

			if n_zero > 0:
				distribution[0] += n_zero

			values = []

			for copy_number, n_proteins in distribution.items():

				cur.execute(""" INSERT INTO feature_copy_distribution VALUES (?,?,?,?) """, (lengthbin, featureid, copy_number, n_proteins))
				values.extend([copy_number] * n_proteins)

			mean_copies = statistics.mean(values)
			median_copies = statistics.median(values)

			if len(values) > 1:
				sd_copies = statistics.stdev(values)
			else:
				sd_copies = 0

			total_copies = sum(values)
			cur.execute(""" INSERT INTO feature_copy_statistics VALUES (?,?,?,?,?,?,?) """, (lengthbin, featureid, mean_copies, median_copies, sd_copies, total_copies, total_proteins ))


def add_indexes(cur):
	cur.execute(""" CREATE INDEX idx_feature_copy_distribution ON feature_copy_distribution(lengthbin, featureid) """)
	cur.execute(""" CREATE INDEX idx_feature_copy_statistics ON feature_copy_statistics(lengthbin, featureid) """)
	cur.execute(""" CREATE INDEX idx_feature_annotation ON feature_annotation(featureid) """)


def main():

	infile = sys.argv[1]
	outfile = sys.argv[2]

	conn, cur = create_database(outfile)

	build_background(infile, cur)

	add_indexes(cur)
	conn.commit()
	conn.close()


if __name__ == "__main__":
	main()
