import os
import pandas as pd
import numpy as np
from collections import defaultdict, Counter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import argparse


def main():
	parser = argparse.ArgumentParser(description="Plot heatmap from bin_pairs.bed")
	parser.add_argument("bin_pairs")
	args = parser.parse_args()

	bin_pairs_dict = defaultdict(dict)
	A_bin_pairs_count = defaultdict(Counter)
	B_bin_pairs_count = defaultdict(Counter)
	for line in open(args.bin_pairs, 'r'):
		row = line.rstrip().split()
		group1, count1 = row[0].split(".")
		group2, count2 = row[1].split(".")

		row[2] = int(row[2])
		row[3] = int(row[3])

		count1 = int(count1)
		count2 = int(count2)

		if count1 > 13:
			count1 = 14

		if count2 > 13:
			count2 = 14

		# if count1 <= 13 and count2 <= 13:
		# 	# bin_pairs_dict[row[0]][row[1]] = float(row[4])
		# else:
		index1 = "{}.{}".format(group1, count1)
		index2 = "{}.{}".format(group2, count2)
		A_bin_pairs_count[index1][index2] += row[2]
		B_bin_pairs_count[index1][index2] += row[3]

	for index1 in A_bin_pairs_count:
		for index2 in A_bin_pairs_count[index1]:
			bin_pairs_dict[index1][index2] = np.log2((A_bin_pairs_count[index1][index2] + A_bin_pairs_count[index2][index1]) / (B_bin_pairs_count[index1][index2] + B_bin_pairs_count[index2][index1]))

	bin_pairs_df = pd.DataFrame.from_dict(bin_pairs_dict)

	order = ["No." + str(x) for x in range(15)] + ["TCF1." + str(x) for x in reversed(range(1,15))]

	bin_pairs_df = bin_pairs_df.loc[order]
	bin_pairs_df = bin_pairs_df[order]
	# bin_pairs_df = bin_pairs_df.sort_index(kind="stable")
	# bin_pairs_df = bin_pairs_df.sort_index(axis=1, kind="stable")

	sns_plot = sns.heatmap(bin_pairs_df, cmap="YlGnBu", vmax=2)
	plt.savefig("Figure3H.pdf")


if __name__ == '__main__':
	main()
