#!/usr/bin/env python3

import csv
import time
import argparse
import logging
from collections import defaultdict
from lib import allparsers
import numpy as np
from scipy.stats import wilcoxon, t, percentileofscore
from statsmodels.stats.multitest import multipletests

def get_config(inputtype):

	if inputtype == "chemproperties":

		comparisons = [("full", "fullOOF", "fullWT"),
                        ("after", "afterOOF", "afterWT"),
                        ("perc10", "perc10OOF", "perc10WT")]

		property_column = "chemproperty"

	elif inputtype == "secondary_structure":

		comparisons = [("full", "Coil", "WT_Coil"),
                        ("full", "Helix", "WT_Helix"),
                        ("full", "E_turn", "WT_E_turn"),
                        ("after", "after_Coil", "after_WT_Coil"),
                        ("after", "after_Helix", "after_WT_Helix"),
                        ("after", "after_E_turn", "after_WT_E_turn")]

		property_column = None

	elif inputtype == "secondary_structure_diff":

		comparisons = [("full", "diff_Coil"),
                        ("full", "diff_Helix"),
                        ("full", "diff_E_turn"),
                        ("after", "after_diff_Coil"),
                        ("after", "after_diff_Helix"),
                        ("after", "after_diff_E_turn")]

		property_column = None

	elif inputtype == "pdb_tmalign":
		comparisons = [("metric", None, None)]
		property_column = None

	elif inputtype == "domains":
		comparisons = []
		property_column = None 

	else:
		raise ValueError(f"Unknown inputtype: {inputtype}" )

	return comparisons, property_column

def get_tmalign_metrics():
	return ["TM_after_WTnorm", "TM_full_WTnorm", "TM_after_aligned_fraction", "after_length_ratio", "OOF_after_conf_mean", "delta_conf_after"]


def calculate_domain_stats(values):

	values = np.array(values, dtype=float)

	n = len(values)

	mean_diff = np.mean(values)
	median_diff = np.median(values)

	if n > 1:
		sd_diff = np.std(values, ddof=1)
	else:
		sd_diff = 0.0

	if n > 1 and sd_diff > 0:
		sem = sd_diff / np.sqrt(n)
		ci_low, ci_high = t.interval(confidence=0.95, df=n-1, loc=mean_diff, scale=sem)

	else:
		ci_low = mean_diff
		ci_high = mean_diff

	return {"N": n, "MeanDifference": mean_diff, "MedianDifference": median_diff, "SDDifference": sd_diff, "CI95_Lower": ci_low, "CI95_Upper": ci_high}

def run_domain_statistics(args):

	rows = []

	region_diffs = defaultdict(list)
	total_diffs = defaultdict(list)
	total_seen = set()

	feature_names = {}
	feature_databases = {}

	region_wt_counts = defaultdict(list)
	region_oof_counts = defaultdict(list)

	total_wt_counts = defaultdict(list)
	total_oof_counts = defaultdict(list)

	infile = open(args.input)

	reader = csv.DictReader(infile, delimiter="\t")

	for line in reader:
		rows.append(line)

		feature = line["FeatureID"]
		region = line["Region"]
		diff = float(line["Difference"])

		region_diffs[(feature, region)].append(diff)
		region_wt_counts[(feature, region)].append(float(line["WT_Count"]))
		region_oof_counts[(feature, region)].append(float(line["OOF_Count"]))

		feature_names[feature] = line["FeatureName"]
		feature_databases[feature] = line["Database"]

		total_key = (line["OOF_ID"],feature)

		if total_key not in total_seen:

			total_seen.add(total_key)
			total_diff = (float(line["OOF_TotalCopies"]) - float(line["WT_TotalCopies"]))
			total_diffs[feature].append(total_diff)
			total_wt_counts[feature].append( float(line["WT_TotalCopies"]))
			total_oof_counts[feature].append(float(line["OOF_TotalCopies"]))

	infile.close()

	region_stats = {}

	for key in region_diffs:
		region_stats[key] = calculate_domain_stats(region_diffs[key])

	total_stats = {}

	for feature in total_diffs:
		total_stats[feature] = calculate_domain_stats(total_diffs[feature])

	outfile = open(args.outfilename + "_domains_feature_region_statistics.tsv","w")
	outfile.write("FeatureID\tDatabase\tFeatureName\tRegion\tN\tMeanDifference\tMedianDifference\tSDDifference\tCI95_Lower\tCI95_Upper\tNumGain\tNumLoss\tNumUnchanged\tMeanWTCount\tMeanOOFCount\n")

	for feature, region in sorted(region_stats):
		stats = region_stats[(feature, region)]
		values = region_diffs[(feature, region)]

		num_gain = sum(1 for x in values if x > 0)
		num_loss = sum(1 for x in values if x < 0)
		num_unchanged = sum(1 for x in values if x == 0)

		mean_wt = np.mean(region_wt_counts[(feature, region)])
		mean_oof = np.mean(region_oof_counts[(feature, region)])

		outfile.write("\t".join([feature,feature_databases[feature],feature_names[feature],region,str(stats["N"]), str(stats["MeanDifference"]),str(stats["MedianDifference"]),str(stats["SDDifference"]),str(stats["CI95_Lower"]),str(stats["CI95_Upper"]),str(num_gain),str(num_loss), str(num_unchanged),str(mean_wt), str(mean_oof)]) + "\n")
	outfile.close()

	outfile = open(args.outfilename + "_domains_feature_total_statistics.tsv","w")
	outfile.write("FeatureID\tDatabase\tFeatureName\tN\tMeanDifference\tMedianDifference\tSDDifference\tCI95_Lower\tCI95_Upper\tNumGain\tNumLoss\tNumUnchanged\tMeanWTTotalCopies\tMeanOOFTotalCopies\n")

	for feature in sorted(total_stats):
		stats = total_stats[feature]
		values = total_diffs[feature]
		num_gain = sum(1 for x in values if x > 0)
		num_loss = sum(1 for x in values if x < 0)
		num_unchanged = sum(1 for x in values if x == 0)

		mean_wt = np.mean(total_wt_counts[feature])
		mean_oof = np.mean(total_oof_counts[feature])

		outfile.write("\t".join([feature, feature_databases[feature], feature_names[feature],str(stats["N"]),str(stats["MeanDifference"]),str(stats["MedianDifference"]),str(stats["SDDifference"]),str(stats["CI95_Lower"]),str(stats["CI95_Upper"]), str(num_gain),str(num_loss),str(num_unchanged),str(mean_wt),str(mean_oof)]) + "\n")
	outfile.close()

	outfilename = (args.outfilename + "_domains_per_oof_statistics.tsv")
	outfile = open(outfilename, "w")

	fieldnames = list(rows[0].keys())
	fieldnames.extend(["Direction","Region_Percentile", "Region_EffectSize","Region_N", "Region_MeanDifference","Region_MedianDifference", "Region_SDDifference", "Region_CI95_Lower","Region_CI95_Upper", "Total_Percentile","Total_EffectSize","Total_N","Total_MeanDifference","Total_MedianDifference","Total_SDDifference", "Total_CI95_Lower", "Total_CI95_Upper"])

	writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter="\t")
	writer.writeheader()

	for row in rows:
		feature = row["FeatureID"]
		region = row["Region"]
		diff = float(row["Difference"])
		total_diff = (float(row["OOF_TotalCopies"]) - float(row["WT_TotalCopies"]))

		if diff > 0:
			direction = "gain"
		elif diff < 0:
			direction = "loss"
		else:
			direction = "unchanged"

		rstats = region_stats[(feature, region)]
		region_percentile = percentileofscore(np.abs(region_diffs[(feature, region)]), abs(diff))

		if rstats["SDDifference"] > 0:
			region_effect = (diff - rstats["MeanDifference"]) / rstats["SDDifference"]
		else:
			if diff == rstats["MeanDifference"]:
				region_effect = 0
			else:
				region_effect = np.nan

		tstats = total_stats[feature]

		total_percentile = percentileofscore(np.abs(total_diffs[feature]),abs(total_diff))

		if tstats["SDDifference"] > 0:
			total_effect = (total_diff - tstats["MeanDifference"]) / tstats["SDDifference"]
		else:
			if total_diff == tstats["MeanDifference"]:
				total_effect = 0
			else:
				total_effect = np.nan

		row["Direction"] = direction
		row["Region_Percentile"] = region_percentile
		row["Region_EffectSize"] = region_effect
		row["Region_N"] = rstats["N"]
		row["Region_MeanDifference"] = rstats["MeanDifference"]
		row["Region_MedianDifference"] = rstats["MedianDifference"]
		row["Region_SDDifference"] = rstats["SDDifference"]
		row["Region_CI95_Lower"] = rstats["CI95_Lower"]
		row["Region_CI95_Upper"] = rstats["CI95_Upper"]
		row["Total_Percentile"] = total_percentile
		row["Total_EffectSize"] = total_effect
		row["Total_N"] = tstats["N"]
		row["Total_MeanDifference"] = tstats["MeanDifference"]
		row["Total_MedianDifference"] = tstats["MedianDifference"]
		row["Total_SDDifference"] = tstats["SDDifference"]
		row["Total_CI95_Lower"] = tstats["CI95_Lower"]
		row["Total_CI95_Upper"] = tstats["CI95_Upper"]

		writer.writerow(row)

	outfile.close()
	global_domain_burden_statistics( rows, args )


def get_datastructure(args):

	comparisons, property_column = get_config(args.inputtype)

	infile = open(args.input, "r")
	readerdict = csv.DictReader(infile, delimiter="\t")

	alllines = []
	global_diffs = defaultdict(list)
	condition_prop_diffs = defaultdict(lambda: defaultdict(list))

	for line in readerdict:
		condition = line["condition"]

		if args.inputtype == "pdb_tmalign":
			for metric in get_tmalign_metrics():
				value = line.get(metric)

				if value in ["", "NA", "None", None]:
					continue

				diff = float(value)

				longline = dict(line)
				longline["property"] = metric
				longline["comparison"] = "metric"
				longline["diff"] = diff

				global_key = ("metric", condition, metric)
				global_diffs[global_key].append(diff)

				local_key = (condition, metric)
				condition_prop_diffs[local_key]["metric"].append(diff)
				alllines.append(longline)

			continue


		if args.inputtype == "chemproperties":
			properties = [(line[property_column], comparisons)]

		elif args.inputtype == "secondary_structure":
			properties = [("Coil", [comparisons[0], comparisons[3]]),
					("Helix", [comparisons[1], comparisons[4]]),
					("E_turn", [comparisons[2], comparisons[5]])]

		elif args.inputtype == "secondary_structure_diff":
			properties = [("Coil", [comparisons[0], comparisons[3]]),
					("Helix", [comparisons[1], comparisons[4]]),
					("E_turn", [comparisons[2], comparisons[5]])]

		for prop, prop_comparisons in properties:

			if args.inputtype == "secondary_structure_diff":

				for comparison, diff_column in prop_comparisons:
					value = line[diff_column]

					if value in ["", "NA", "None"]:
						continue

					diff = float(value)
					longline = dict(line)
					longline["diff"] = diff

					global_key = (comparison, condition, prop)
					global_diffs[global_key].append(diff)

					local_key = (condition, prop)
					condition_prop_diffs[local_key][comparison].append(diff)

					longline["property"] = prop
					longline["comparison"] = comparison

					alllines.append(longline)

			else:

				for comparison, oof_column, wt_column in prop_comparisons:
					diff = float(line[oof_column]) - float(line[wt_column])
					longline = dict(line)
					longline["OOF_value"] = line[oof_column]
					longline["WT_value"] = line[wt_column]
					longline["diff"] = diff

					global_key = (comparison, condition, prop)
					global_diffs[global_key].append(diff)

					local_key = (condition, prop)
					condition_prop_diffs[local_key][comparison].append(diff)

					longline["property"] = prop
					longline["comparison"] = comparison
					alllines.append(longline)

	return alllines, global_diffs, condition_prop_diffs

def global_statistics(global_diffs, args):

	results = []
	all_pvalues = []

	outfilename = args.outfilename + "_" + args.inputtype + "_global_stats.txt"
	outfile = open(outfilename, "w")

	if args.inputtype == "pdb_tmalign":
		outfile.write("comparison\tcondition\tproperty\tN\tmean_value\tmedian_value\tsd\tci_low\tci_high\tcohens_d\tpvalue\tadj_pvalue\n")
	else:
		outfile.write("comparison\tcondition\tproperty\tN\tmean_diff\tmedian_diff\tsd\tci_low\tci_high\tcohens_d\tpvalue\tadj_pvalue\n")                      

	for (comparison, condition, prop), diffs in sorted(global_diffs.items()):

		diffs = np.array(diffs, dtype=float)
		n = len(diffs)

		mean_diff = np.mean(diffs)
		median_diff = np.median(diffs)
		sd_diff = np.std(diffs, ddof=1)
		cohens_d = np.nan

		if sd_diff > 0:
			cohens_d = mean_diff / sd_diff

		if n > 1:
			sem = sd_diff / np.sqrt(n)
			ci_low, ci_high = t.interval(confidence=0.95, df=n - 1, loc=mean_diff, scale=sem)
		else:
			ci_low = np.nan
			ci_high = np.nan

		pvalue = np.nan

		if args.inputtype == "pdb_tmalign":
			pvalue = np.nan
		elif n > 1 and np.any(diffs != 0):
			pvalue = wilcoxon(diffs)[1]

		all_pvalues.append(pvalue)

		results.append({"comparison": comparison,"condition": condition, "property": prop, "N": n, "mean_diff": mean_diff, "median_diff": median_diff, "sd": sd_diff, "ci_low": ci_low, "ci_high": ci_high, "cohens_d": cohens_d, "pvalue": pvalue})

	valid_p = [x for x in all_pvalues if not np.isnan(x)]

	if len(valid_p) > 0:
		adj_p = multipletests( valid_p, method="fdr_bh")[1]
	else:
		adj_p = []

	idx = 0

	for res in results:
		if np.isnan(res["pvalue"]):
			res["adj_pvalue"] = np.nan
		else:
			res["adj_pvalue"] = adj_p[idx]
			idx += 1

		outfile.write("\t".join([str(res["comparison"]),str(res["condition"]),str(res["property"]), str(res["N"]),str(res["mean_diff"]),str(res["median_diff"]),str(res["sd"]),str(res["ci_low"]), str(res["ci_high"]), str(res["cohens_d"]), str(res["pvalue"]),str(res["adj_pvalue"])]) + "\n")

	outfile.close()


def global_domain_burden_statistics(rows, args):

	fieldnames = ["Region", "Metric","N","MeanDifference", "MedianDifference",  "SDDifference", "CI95_Lower","CI95_Upper","EffectSize","WilcoxonP", "BH_P" ]

	per_oof = defaultdict(lambda: defaultdict(list))

	total_counts = {}

	for row in rows:
		oof_id = row["OOF_ID"]
		region = row["Region"]
		diff = float(row["Difference"])

		per_oof[region][oof_id].append(diff)

		if oof_id not in total_counts:

			total_counts[oof_id] = (float(row["OOF_TotalCopies"]) - float(row["WT_TotalCopies"]))

	region_results = []
	total_diffs = list(total_counts.values())
	region_results.append(("Total","DomainCopies",total_diffs))

	for region in ["before", "after", "spanning", "truncated"]:
		if region not in per_oof:
			continue

		net_diffs = []
		gain_counts = []
		loss_counts = []

		for oof_id in per_oof[region]:

			diffs = per_oof[region][oof_id]
			net_diffs.append(sum(diffs))

			gain_counts.append(sum(1 for x in diffs if x > 0 ))
			loss_counts.append(sum(1 for x in diffs if x < 0 ))

		region_results.append((region,"NetDifference",net_diffs))
		region_results.append((region,"GainCount",gain_counts))
		region_results.append((region,"LossCount",loss_counts))

	all_pvalues = []
	temp_results = []

	for region, metric, values in region_results:

		values = np.array(values, dtype = float)
		n = len(values)
		mean_diff = np.mean(values)
		median_diff = np.median(values)

		if n > 1:
			sd_diff = np.std(values, ddof=1)
		else:
			sd_diff = 0.0

		if n > 1 and sd_diff > 0:
			sem = (sd_diff / np.sqrt(n))
			ci_low, ci_high = t.interval(confidence = 0.95, df = n - 1, loc = mean_diff, scale = sem)

		else:
			ci_low = mean_diff
			ci_high = mean_diff

		if n > 1 and np.any(values != 0):
			pvalue = wilcoxon(values).pvalue
		else:
			pvalue = np.nan


		all_pvalues.append(pvalue)

		if sd_diff > 0:
			effect_size = (mean_diff / sd_diff)
		else: 
			effect_size = np.nan 

		temp_results.append({"Region":region, "Metric":metric, "N":n, "MeanDifference":mean_diff, "MedianDifference":median_diff, "SDDifference":sd_diff, "CI95_Lower":ci_low, "CI95_Upper":ci_high, "EffectSize":effect_size, "WilcoxonP":pvalue})

	valid_p = [     x for x in all_pvalues if not np.isnan(x)]

	if len(valid_p) > 0:
		adj_p = multipletests(       valid_p,     method = "fdr_bh"   )[1]
	else:
		adj_p = []

	idx = 0

	for result in temp_results:

		if np.isnan(  result[ "WilcoxonP" ]):
			result[ "BH_P"] = np.nan
		else:
			result[ "BH_P" ] = adj_p[idx]
			idx += 1

	outfile = open(args.outfilename + "_domain_burden_statistics.tsv","w")

	writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter="\t")
	writer.writeheader()

	for result in temp_results:
		writer.writerow(result)
	outfile.close()






def per_peptide_statistics(condition_prop_diffs, alllines, args):

	means = {}
	sds = {}
	percentiles = {}

	outfilename = args.outfilename + "_" + args.inputtype + "_perpeptide_stats.txt"
	outfile = open(outfilename, "w")

	for local_key in condition_prop_diffs:

		means[local_key] = {}
		sds[local_key] = {}
		percentiles[local_key] = {}

		for comparison in condition_prop_diffs[local_key]:
			vals = np.array(condition_prop_diffs[local_key][comparison],dtype=float)

			means[local_key][comparison] = np.mean(vals)
			sds[local_key][comparison] = np.std(vals, ddof=1)
			background = [round(abs(x), 8) for x in vals]
			percentiles[local_key][comparison] = {}
			sorted_bg = np.sort(background)
			unique_vals, counts = np.unique(sorted_bg ,return_counts=True)
			cumulative = np.cumsum(counts)

			for value, rank in zip(unique_vals, cumulative):
				percentiles[local_key][comparison][value] = (rank / len(sorted_bg)) * 100.0


	for line in alllines:

		local_key = (line["condition"], line["property"])
		comparison = line["comparison"]
		diff = float(line["diff"])
		sd = sds[local_key][comparison]
		mean = means[local_key][comparison]

		if sd > 0:
			zscore = (diff - mean) / sd
		else:
			zscore = 0.0

		line["zscore"] = zscore
		line["abs_zscore"] = abs(zscore)

		key = round(abs(diff), 8)
		line["percentile"] = percentiles[local_key][comparison].get(key, np.nan)

	fieldnames = list(alllines[0].keys())

	writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter="\t")
	writer.writeheader()

	for line in alllines:
		writer.writerow(line)

	outfile.close()


def main():

	logging.basicConfig(level = logging.INFO , format='{asctime} - {levelname} - {message}',filemode='w', style = "{", datefmt = "[%d-%m-%Y] [%H:%M] ")
	logging.info("start run")

	starttime = time.time()

	parser = argparse.ArgumentParser()
	allparsers.parser_getstats(parser)
	args = parser.parse_args()

	if args.inputtype == "domains":
		run_domain_statistics(args)
	else:
		alllines, global_diffs, condition_prop_diffs = get_datastructure(args)

		global_statistics(global_diffs, args)
		per_peptide_statistics(condition_prop_diffs, alllines, args)

	endtime = round((time.time() - starttime) / 60, 2)

	infotext = "run finished in " + str(endtime) + " m"
	logging.info(infotext)

main()
