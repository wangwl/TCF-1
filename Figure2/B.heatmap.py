import os
import pandas as pd
from collections import defaultdict,Counter
import matplotlib
matplotlib.use('agg')
matplotlib.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import argparse


def make_heatmap(infile, outdir):
	bin_df = pd.read_csv(infile, index_col = 0)
	CTCF_df = bin_df.loc[(bin_df["DP_CTCF"] > 0) & (bin_df['Thy_TCF1'] == 0)].mean()
	TCF1_df = bin_df.loc[(bin_df['Thy_TCF1'] > 0) & (bin_df['DP_CTCF'] == 0)].mean()
	CTCF1_df = bin_df.loc[(bin_df["DP_CTCF"] > 0) & (bin_df['Thy_TCF1'] > 0)].mean()
	noCTCF1_df = bin_df.loc[(bin_df['Thy_TCF1'] == 0) & (bin_df['DP_CTCF'] == 0)].mean()

	is_dict = defaultdict(dict)
	Cells = ["CLP", "ETP", "DN2", "DN3", "DN4", "DP"]

	for i, cell in enumerate(Cells):
		is_dict['CTCF+TCF-1-'][cell] = CTCF_df[cell]
		is_dict['CTCF-TCF-1+'][cell] = TCF1_df[cell]
		is_dict['CTCF+TCF-1+'][cell] = CTCF1_df[cell]
		is_dict['CTCF-TCF-1-'][cell] = noCTCF1_df[cell]

	is_df = pd.DataFrame.from_dict(is_dict, orient='index')
	is_df = is_df.sort_index()

	plt.figure(figsize=[6,6])
	# is_df.dropna(inplace=True, how='all', axis=0)
	# is_df.index.name = "TFs"
	# is_df.to_csv("TFs_insulation.csv")

	data_df = is_df[["CLP", "ETP", "DN2", "DN3", "DN4", "DP"]]
	# zscore = stats.zscore(data_df.values, axis=1, ddof=1)
	# zscore_df = pd.DataFrame(data=zscore, index=data_df.index, columns=data_df.columns)
	# zscore_df.dropna(inplace=True)
	# zscore_df.to_csv("TFs_insulation.zscore.csv")
	# snsplot = sns.heatmap(zscore_df[["CLP","ETP", "DN2", "DN3", "DN4", "DP"]], cmap='viridis')
	# plt.savefig("TF_cells.zscore.pdf")
	data_df = data_df.T[["CTCF+TCF-1-", "CTCF+TCF-1+", "CTCF-TCF-1+", "CTCF-TCF-1-"]]
	plt.figure(figsize=[6,6])
	snsplot = sns.heatmap(data_df, cmap='PiYG', center=0)
	plt.savefig(f"{outdir}/Figure2B.IS.pdf")


def regression(infile):
	bin_df = pd.read_csv(infile, index_col=0)
	bin_df.dropna(inplace=True)



def main():
	parser = argparse.ArgumentParser(description="plot heatmap to show IS")
	parser.add_argument("infile")
	parser.add_argument("--outdir", default=".")
	args = parser.parse_args()

	make_heatmap(args.infile, args.outdir)


if __name__ == '__main__':
	main()
