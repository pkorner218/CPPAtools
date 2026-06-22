#!/usr/bin/env python3

import csv
import argparse
import logging
import time 
import numpy as np

from collections import defaultdict
from statsmodels.stats.multitest import multipletests

from lib import allparsers


def load_weights(weightfile):
	weights = {}
	infile = open(weightfile)
	reader = csv.DictReader(infile, delimiter="\t")

	for line in reader:
		database = line["Database"]
		feature_class = line["FeatureClass"]
		weight = float(line["Weight"])
		weights[(database, feature_class)] = weight

	infile.close()

	return weights


FEATURE_CLASS_MAP = {
	"CDD": "Domain",
	"Pfam": "Domain",
	"Gene3D": "Domain",
	"PANTHER": "Domain",
	"PIRSF": "Domain",
	"SMART": "Domain",
	"SUPERFAMILY": "Domain",
	"PRINTS": "Motif",
	"ProSiteProfiles": "Motif",
	"ProSitePatterns": "Motif",
	"TMHMM": "Transmembrane",
	"Coils": "CoiledCoil",
	"MobiDBLite": "Disorder"
}


ADDED_FIELDS = [
	"EmpiricalP_BH",
	"FeatureClass",
	"ToolWeight",
	"RankEligible",
	"AbsentBackgroundScored",
	"BackgroundScoringMode",
	"BaseRankScore",
	"alternativeWeightedRankScore",
	"ConsensusWeightedRankScore",
	"WeightedRankScore",
	"ConsensusGroup",
	"ConsensusScore",
	"ConsensusSupportCount",
	"ConsensusSupportTools",
	"ScoreMode",
	"BackgroundEffectSize",
	"BackgroundSignificant"
]


def is_missing(value):
	return value in ["NA", "", None, "None"]


def parse_float(value, default=0.0):
	if is_missing(value):
		return default
	try:
		return float(value)
	except ValueError:
		return default


def clean_tsv_field(value):

	if value is None:
		return "NA"

	value = str(value)

	if value == "":
		return "NA"

	value = value.replace("\t", " ")
	value = value.replace("\n", " ")
	value = value.replace("\r", " ")

	while "  " in value:
		value = value.replace("  ", " ")

	value = value.strip()

	if value == "":
		return "NA"

	return value


def clean_output_row(row, fieldnames):

	return {field: clean_tsv_field(row.get(field,"NA")) for field in fieldnames}


def get_direction(difference):

	if difference > 0:
		return "gain"

	if difference < 0:
		return "loss"

	return "unchanged"


def safe_group_value(row, column):

	value = row.get(column, "")

	if is_missing(value):
		return row.get("FeatureID", "NA")

	return value


def collect_pvalues(inputfile):

	pvalues = []
	valid_flags = []
	infile = open(inputfile,"r")
	reader = csv.DictReader(infile,delimiter="\t")

	for row in reader:
		p = row["EmpiricalP"]

		if p in ["NA", "", None]:
			valid_flags.append(False)
		else:
			pvalues.append(float(p))
			valid_flags.append(True)

	infile.close()

	if len(pvalues) > 0:
		adj_pvalues = list(multipletests(pvalues,method="fdr_bh")[1])
	else:
		adj_pvalues = []

	return valid_flags, adj_pvalues


def get_bh_for_row(valid_flags, adj_pvalues, row_index, p_index):

	if not valid_flags[row_index]:
		return "NA", p_index

	p = adj_pvalues[p_index]
	p_index += 1

	return p, p_index


def initialize_scoring_fields(row, weights, args):

	database = row["Database"]
	feature_class = FEATURE_CLASS_MAP.get(database,"Unknown")

	row["FeatureClass"] = feature_class

	tool_weight = weights.get((database,feature_class),0.0)

	row["ToolWeight"] = tool_weight
	row["RankEligible"] = False

	row["AbsentBackgroundScored"] = False
	row["BackgroundScoringMode"] = "observed_background"
	row["BaseRankScore"] = 0.0
	row["alternativeWeightedRankScore"] = 0.0
	row["ConsensusWeightedRankScore"] = 0.0
	row["WeightedRankScore"] = 0.0
	row["ConsensusScore"] = 0.0
	row["ConsensusSupportCount"] = 0
	row["ConsensusSupportTools"] = "NA"
	row["ConsensusGroup"] = safe_group_value(row,args.consensus_group_column)

	if args.use_alternative_tool_weight_score:
		row["ScoreMode"] = "alternative_tool_weight"
	else:
		row["ScoreMode"] = "consensus"

	return row, tool_weight


def calculate_row_background_scores(row, tool_weight, args):

	if row["BackgroundSource"] == "Absent":

		diff = parse_float(row["Difference"])

		row["BackgroundScoringMode"] = "absent_not_scored"

		if (args.score_absent_background and diff != 0 and tool_weight > 0 ):
			effect = (diff * args.absent_effect_per_copy)
			p = args.absent_pvalue
			row["BackgroundEffectSize"] = effect
			row["EmpiricalP_BH"] = p

			row["BackgroundSignificant"] = True

			base_score = (abs(effect) * (-np.log10(p + 1e-300)))

			row["BaseRankScore"] = base_score
			row["alternativeWeightedRankScore"] = (base_score * tool_weight)
			row["RankEligible"] = True
			row["AbsentBackgroundScored"] = True
			row["BackgroundScoringMode"] = "absent_background_pseudo_score"

		else:
			row["BackgroundEffectSize"] = "NA"
			row["BackgroundSignificant"] = False

		return row

	if row["BackgroundSD"] in ["NA", "", None]:
		row["BackgroundEffectSize"] = "NA"
		row["BackgroundSignificant"] = False

		return row

	sd = parse_float(row["BackgroundSD"])

	if sd <= 0:
		row["BackgroundEffectSize"] = "NA"
		row["BackgroundSignificant"] = False

		return row

	diff = parse_float(row["Difference"])
	mean = parse_float(row["BackgroundMean"])
	effect = (diff - mean) / sd

	row["BackgroundEffectSize"] = effect
	p = row["EmpiricalP_BH"]

	if p == "NA":
		row["BackgroundSignificant"] = False
		return row

	p = float(p)
	significant = (p < 0.05)

	row["BackgroundSignificant"] = significant

	base_score = (abs(effect) * (-np.log10(p + 1e-300)))

	row["BaseRankScore"] = base_score
	row["alternativeWeightedRankScore"] = (base_score * tool_weight)

	if tool_weight > 0:
		row["RankEligible"] = True

	return row


def build_consensus_support(inputfile, weights, valid_flags, adj_pvalues, args):

	support_by_group = defaultdict(dict)

	infile = open(inputfile,"r")

	reader = csv.DictReader(infile,delimiter="\t")

	p_index = 0

	for row_index, row in enumerate(reader):
		p, p_index = get_bh_for_row(valid_flags,adj_pvalues,row_index,p_index)

		row["EmpiricalP_BH"] = p

		row, tool_weight = initialize_scoring_fields(row,weights,args)
		row = calculate_row_background_scores(row,tool_weight,args)

		if not row["RankEligible"]:
			continue

		diff = parse_float(row["Difference"])

		if diff == 0:
			continue

		direction = get_direction(diff)

		key = (row["OOF_ID"],row["Region"],direction,row["ConsensusGroup"])
		tool_key = (row["Database"],row["FeatureClass"])

		if tool_weight <= 0:
			continue

		if tool_key not in support_by_group[key]:
			support_by_group[key][tool_key] = tool_weight

		elif tool_weight > support_by_group[key][tool_key]:
			support_by_group[key][tool_key] = tool_weight

	infile.close()

	consensus_lookup = {}

	for key, support in support_by_group.items():

		product_not_correct = 1.0
		support_labels = []

		for tool_key in sorted(support):
			weight = support[tool_key]

			product_not_correct *= (1.0 - weight)
			support_labels.append(tool_key[0] + "|" + tool_key[1] + "=" + str(weight))

		consensus_score = (1.0 - product_not_correct )

		consensus_lookup[key] = {"ConsensusScore": consensus_score, "ConsensusSupportCount": len(support), "ConsensusSupportTools": ";".join(support_labels)}

	return consensus_lookup


def update_row_with_consensus(row, consensus_lookup, args):

	diff = parse_float(row["Difference"])

	if diff != 0 and row["RankEligible"]:

		direction = get_direction(diff)

		key = (row["OOF_ID"], row["Region"], direction, row["ConsensusGroup"])

		if key in consensus_lookup:
			info = consensus_lookup[key]
			row["ConsensusScore"] = info["ConsensusScore"]
			row["ConsensusSupportCount"] = info["ConsensusSupportCount"]
			row["ConsensusSupportTools"] = info["ConsensusSupportTools"]

	consensus_weighted_score = (parse_float(row["BaseRankScore"]) * parse_float(row["ConsensusScore"]))
	row["ConsensusWeightedRankScore"] = consensus_weighted_score

	if args.use_alternative_tool_weight_score:
		row["WeightedRankScore"] = row["alternativeWeightedRankScore"]
	else:
		row["WeightedRankScore"] = consensus_weighted_score

	return row



class OOFConsensusCountAggregator:

	def __init__(self):
		self.groups = {}

	def _empty(self, row, key):
		return {"ENSTID": row.get("ENSTID", "NA"), "Gene": row.get("Gene", "NA"), "WT_ID": row.get("WT_ID", "NA"), "OOF_ID": row.get("OOF_ID", "NA"), "Region": key[4], "WT_TotalDomainCopies_raw": 0.0, "OOF_TotalDomainCopies_raw": 0.0, "WT_FeatureIDs": set(), "OOF_FeatureIDs": set(), "WT_ConsensusGroups": set(), "OOF_ConsensusGroups": set()}

	def add(self, row):
		region = row.get("Region","NA")
		keys = [(row.get("ENSTID", "NA"), row.get("Gene", "NA"), row.get("WT_ID", "NA"), row.get("OOF_ID", "NA"), region ),(row.get("ENSTID", "NA"), row.get("Gene", "NA"), row.get("WT_ID", "NA"), row.get("OOF_ID", "NA"), "Total" )]

		wt_count = parse_float(row.get("WT_Count", 0))
		oof_count = parse_float(row.get("OOF_Count", 0))

		feature_id = row.get("FeatureID","NA")
		consensus_group = row.get("ConsensusGroup",feature_id)

		if is_missing(consensus_group):
			consensus_group = feature_id

		for key in keys:

			if key not in self.groups:
				self.groups[key] = self._empty(row,key)

			group = self.groups[key]
			group["WT_TotalDomainCopies_raw"] += wt_count
			group["OOF_TotalDomainCopies_raw"] += oof_count

			if wt_count > 0:
				group["WT_FeatureIDs"].add(feature_id)
				group["WT_ConsensusGroups"].add(consensus_group)

			if oof_count > 0:
				group["OOF_FeatureIDs"].add(feature_id)
				group["OOF_ConsensusGroups"].add(consensus_group)

	def results(self):

		output = []
		region_order = {"before": 0,"spanning": 1,"truncated": 2,"after": 3,"Total": 4}

		for key in sorted(self.groups,key=lambda x: (x[0],x[3],region_order.get(x[4], 99))):
			group = self.groups[key]

			wt_copies = group["WT_TotalDomainCopies_raw"]
			oof_copies = group["OOF_TotalDomainCopies_raw"]
			wt_features = len(group["WT_FeatureIDs"])

			oof_features = len(group["OOF_FeatureIDs"])
			wt_consensus = len(group["WT_ConsensusGroups"])
			oof_consensus = len(group["OOF_ConsensusGroups"])

			output.append({"ENSTID": group["ENSTID"],"Gene": group["Gene"],"WT_ID": group["WT_ID"],"OOF_ID": group["OOF_ID"],"Region": group["Region"],"WT_TotalDomainCopies_raw": wt_copies,"OOF_TotalDomainCopies_raw": oof_copies,"DomainCopyDifference_raw": oof_copies - wt_copies,"WT_UniqueFeatures_raw": wt_features,"OOF_UniqueFeatures_raw": oof_features,"UniqueFeatureDifference_raw": oof_features - wt_features,"WT_UniqueConsensusDomains_raw": wt_consensus,"OOF_UniqueConsensusDomains_raw": oof_consensus,"ConsensusDomainDifference_raw": oof_consensus - wt_consensus})

		return output


def write_per_oof_consensus_counts(count_rows, args):

	outfilename = (args.outfilename + "_per_oof_consensus_counts.tsv")
	fieldnames = ["ENSTID","Gene","WT_ID","OOF_ID","Region","WT_TotalDomainCopies_raw","OOF_TotalDomainCopies_raw","DomainCopyDifference_raw","WT_UniqueFeatures_raw","OOF_UniqueFeatures_raw","UniqueFeatureDifference_raw","WT_UniqueConsensusDomains_raw","OOF_UniqueConsensusDomains_raw","ConsensusDomainDifference_raw"]

	outfile = open(outfilename,"w")

	writer = csv.DictWriter(outfile,fieldnames=fieldnames,delimiter="\t",lineterminator="\n")
	writer.writeheader()

	for row in count_rows:
		writer.writerow(clean_output_row(row,fieldnames))
	outfile.close()


class FeatureAggregator:

	def __init__(self, args):
		self.args = args
		self.groups = {}

	def _empty(self, row, key):
		return {"key": key, "FeatureID": key[0], "Region": key[1], "N": 0, "SumDifference": 0.0, "Diffs": [], "NumGain": 0, "NumLoss": 0, "NumUnchanged": 0, "Effects": [], "NumSignificant": 0, "NumRankEligible": 0, "SumBaseRankScore": 0.0, "SumalternativeWeightedRankScore": 0.0, "SumConsensusWeightedRankScore": 0.0, "SumWeightedRankScore": 0.0, "SumConsensusScore": 0.0, "RankN": 0, "FeatureIDs": set(), "Databases": set(), "FeatureNames": set(), "FeatureClasses": set(), "ConsensusGroups": set(), "ConsensusSupportTools": set(), "ToolWeight": row.get("ToolWeight", 0.0), "ScoreMode": row.get("ScoreMode", "consensus")}

	def add(self, row):
		if self.args.feature_summary_level == "consensus_group":
			key = (row["ConsensusGroup"],row["Region"])

		else:
			key = (row["FeatureID"],row["Region"])

		if key not in self.groups:
			self.groups[key] = self._empty(row, key)

		group = self.groups[key]
		diff = parse_float(row["Difference"])
		group["N"] += 1
		group["SumDifference"] += diff
		group["Diffs"].append(diff)

		if diff > 0:
			group["NumGain"] += 1
		elif diff < 0:
			group["NumLoss"] += 1
		else:
			group["NumUnchanged"] += 1

		if row["BackgroundEffectSize"] != "NA":
			group["Effects"].append(parse_float(row["BackgroundEffectSize"]))

		if str(row["BackgroundSignificant"]) in ["True", "TRUE", "true", "1"]:
			group["NumSignificant"] += 1

		if row["RankEligible"]:
			group["NumRankEligible"] += 1
			group["RankN"] += 1
			group["SumBaseRankScore"] += parse_float(row["BaseRankScore"])
			group["SumalternativeWeightedRankScore"] += parse_float(row["alternativeWeightedRankScore"])
			group["SumConsensusWeightedRankScore"] += parse_float(row["ConsensusWeightedRankScore"])
			group["SumWeightedRankScore"] += parse_float(row["WeightedRankScore"])
			group["SumConsensusScore"] += parse_float(row["ConsensusScore"])

		group["FeatureIDs"].add(row["FeatureID"])
		group["Databases"].add(row["Database"])
		group["FeatureNames"].add(row["FeatureName"])
		group["FeatureClasses"].add(row["FeatureClass"])
		group["ConsensusGroups"].add(row["ConsensusGroup"])

		if row["ConsensusSupportTools"] != "NA":
			group["ConsensusSupportTools"].add(row["ConsensusSupportTools"])

	def results(self):
		output = []

		for key in self.groups:
			group = self.groups[key]

			n = group["N"]
			rank_n = group["RankN"]
			diffs = group["Diffs"]

			if len(group["Effects"]) > 0:
				mean_effect = float(np.mean(group["Effects"]))
			else:
				mean_effect = np.nan

			if rank_n > 0:
				mean_base = group["SumBaseRankScore"] / rank_n
				mean_alternative = group["SumalternativeWeightedRankScore"] / rank_n
				mean_consensus_weighted = group["SumConsensusWeightedRankScore"] / rank_n
				mean_weighted = group["SumWeightedRankScore"] / rank_n
				mean_consensus = group["SumConsensusScore"] / rank_n
			else:
				mean_base = 0.0
				mean_alternative = 0.0
				mean_consensus_weighted = 0.0
				mean_weighted = 0.0
				mean_consensus = 0.0

			output.append({"FeatureID": group["FeatureID"],"Database": ";".join(sorted(group["Databases"])),"FeatureName": ";".join(sorted(group["FeatureNames"])),"FeatureClass": ";".join(sorted(group["FeatureClasses"])),"ConsensusGroup": ";".join(sorted(group["ConsensusGroups"])),"ConsensusSupportTools": ";".join(sorted(group["ConsensusSupportTools"])) if len(group["ConsensusSupportTools"]) > 0 else "NA","Region": group["Region"],"N": n,"MeanDifference": group["SumDifference"] / n,"MedianDifference": float(np.median(diffs)),"MeanAbsDifference": float(np.mean([abs(x) for x in diffs])),"NumGain": group["NumGain"],"NumLoss": group["NumLoss"],"NumUnchanged": group["NumUnchanged"],"MeanBackgroundEffect": mean_effect,"NumSignificant": group["NumSignificant"],"FractionSignificant": group["NumSignificant"] / n,"NumRankEligible": group["NumRankEligible"],"FractionRankEligible": group["NumRankEligible"] / n,"ToolWeight": group["ToolWeight"],"MeanConsensusScore": mean_consensus,"MeanBaseRankScore": mean_base,"MeanalternativeWeightedRankScore": mean_alternative,"MeanConsensusWeightedRankScore": mean_consensus_weighted,"MeanWeightedRankScore": mean_weighted,"FeatureIDsInGroup": ";".join(sorted(group["FeatureIDs"])),"ScoreMode": group["ScoreMode"]})

		return output


def write_ranked_output(feature_stats, args):

	ranked = sorted(
		feature_stats,
		key=lambda x:(-x["MeanWeightedRankScore"],-x["FractionSignificant"],-x["MeanAbsDifference"]))

	outfilename = (args.outfilename + "_feature_ranked.tsv" )

	fieldnames = ["Rank","FeatureID","Database","FeatureName","FeatureClass","ConsensusGroup","ConsensusSupportTools","Region","N","MeanDifference","MedianDifference","MeanAbsDifference","NumGain","NumLoss","NumUnchanged","MeanBackgroundEffect","NumSignificant","FractionSignificant","NumRankEligible","FractionRankEligible","ToolWeight","MeanConsensusScore","MeanBaseRankScore","MeanalternativeWeightedRankScore","MeanConsensusWeightedRankScore","MeanWeightedRankScore","FeatureIDsInGroup","ScoreMode"]

	outfile = open(outfilename,"w")

	writer = csv.DictWriter(outfile,fieldnames=fieldnames,delimiter="\t",lineterminator="\n")
	writer.writeheader()

	rank = 1

	for row in ranked:
		row["Rank"] = rank
		writer.writerow(clean_output_row(row,fieldnames))
		rank += 1

	outfile.close()


def write_per_oof_and_collect_stats(inputfile, weights, valid_flags, adj_pvalues, consensus_lookup, args):

	outfilename = (args.outfilename + "_per_oof_background_weighted.tsv" )

	infile = open(inputfile,"r")
	reader = csv.DictReader(infile,delimiter="\t")
	input_fields = list(reader.fieldnames)

	added_fields = [field for field in ADDED_FIELDS if field not in input_fields]

	outfile = open(outfilename, "w")

	writer = csv.DictWriter(outfile,fieldnames=input_fields + added_fields,delimiter="\t",lineterminator="\n")
	writer.writeheader()

	aggregator = FeatureAggregator(args)
	count_aggregator = OOFConsensusCountAggregator()
	p_index = 0

	for row_index, row in enumerate(reader):

		p, p_index = get_bh_for_row(valid_flags,adj_pvalues,row_index,p_index)
		row["EmpiricalP_BH"] = p

		row, tool_weight = initialize_scoring_fields(row,weights,args)
		row = calculate_row_background_scores(row,tool_weight,args)
		row = update_row_with_consensus(row,consensus_lookup,args)

		writer.writerow(clean_output_row(row,input_fields + added_fields))

		aggregator.add(row)
		count_aggregator.add(row)

	infile.close()
	outfile.close()

	return aggregator.results(), count_aggregator.results()


def main():


	parser=argparse.ArgumentParser()
	allparsers.parser_consensus_weighting(parser)
	args = parser.parse_args()

	logging.basicConfig(level = logging.INFO , format='{asctime} - {levelname} - {message}',filemode='w', style = "{", datefmt = "[%d-%m-%Y] [%H:%M] ")
	starttime = time.time()
	logging.info("start run")
	logging.info("reading input")


	logging.info("Loading weights")

	weights = load_weights(args.weights)

	logging.info("Pass 1: collecting p-values for BH correction")

	valid_flags, adj_pvalues = collect_pvalues(args.inputfile)

	logging.info("Pass 2: collecting consensus support")

	consensus_lookup = build_consensus_support(args.inputfile,weights,valid_flags,adj_pvalues,args)

	logging.info("Pass 3: writing per-OOF output and feature ranking statistics")

	feature_stats, consensus_count_rows = write_per_oof_and_collect_stats(args.inputfile,weights,valid_flags,adj_pvalues,consensus_lookup,args)

	logging.info("Writing per-OOF consensus count output")

	write_per_oof_consensus_counts(consensus_count_rows,args)

	logging.info("Writing ranked output")

	write_ranked_output(feature_stats,args)
	endtime = round((time.time() - starttime) / 60, 2)
	logging.info("run finished in %s m", endtime)


if __name__ == "__main__":
	main()