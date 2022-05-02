import os
import pandas as pd
import numpy as np
from scipy import stats
from collections import defaultdict

import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt
import seaborn as sns
import argparse


def ks_test(data, outdir):
	data_df = pd.read_csv(data)
	group_df = data_df.groupby('bound')
	Group_H = group_df.get_group('H')

	KS_Pvalue = defaultdict(dict)
	for group in list('ABCDEFGJ'):
		data = group_df.get_group(group)

		for feature in data.columns:
			if feature == 'bound':
				continue

			control_data = Group_H[feature].to_list()
			group_data = data[feature].to_list()

			stat, pvalue = stats.kstest(group_data, control_data, alternative='less')
			KS_Pvalue[group][feature] = -1 * np.log10(pvalue)

	KS_df = pd.DataFrame.from_dict(KS_Pvalue)

	fig = plt.figure(figsize=(6, 12))
	sns.heatmap(KS_df, vmin=0, vmax=10, center=2, cmap='bwr')
	plt.savefig(f"{outdir}/FigureS3I.pdf")



def main():
	parser = argparse.ArgumentParser(description="Perform KS test")
	parser.add_argument("data")
	parser.add_argument("-d", default=".")
	args = parser.parse_args()

	ks_test(args.data, args.d)


if __name__ == '__main__':
	main()
