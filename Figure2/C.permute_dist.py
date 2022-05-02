import os
import sys
import re
import pandas as pd
import numpy as np
from subprocess import check_output
from collections import defaultdict
from pybedtools import BedTool
import argparse


def permute_distance(bedfile, peaks, size=50, n=1000):
	name = re.sub(".bed", "", os.path.basename(bedfile))

	genesBed = BedTool(bedfile)
	peaksBed = BedTool(peaks)
	for i in range(n):
		sampleBed = genesBed.sample(n=size)
		sampleBed = sampleBed.sort()
		closest = sampleBed.closest(peaksBed, d=True, t='first')
		distance = list()
		for row in closest:
			dist = int(row[6])
			if dist == -1:
				continue
			distance.append(dist)

		# print(distance)
		average_dist = np.mean(distance)
		print(f'{name}\t{average_dist}')


def main():
	parser = argparse.ArgumentParser(description="permutate Genes to calculate distance to peaks")
	parser.add_argument("peakfile", help="The set of peaks to calculate distance to genes")
	parser.add_argument("--genes", nargs="+", help="bed files of gene TSS")
	parser.add_argument("--outdir", default=".", help="The output directory")
	args = parser.parse_args()

	for bedfile in args.genes:
		permute_distance(bedfile, args.peakfile)


if __name__ == '__main__':
	main()
