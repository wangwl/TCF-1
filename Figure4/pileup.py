from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
from collections import defaultdict
import multiprocess as mp
import numpy as np
import pandas as pd
from scipy import stats
import argparse
import bioframe
import sparse
import cooltools
import cooler
import os
import re

mpl.style.use('seaborn-white')


def pileup(clrs, conditions, windows, supports, center, flank, gs, outdir, name, row_i):
	from cooltools import snipping
	from cooltools.expected import diagsum
	from enrichment import get_inter

	piles = {}
	sites_enrich = defaultdict(dict)
	inter_matrix = defaultdict(dict)
	vmin = 0
	vmax = 0
	with mp.Pool(10) as pool:
		for cond in conditions:
			cis_exp = diagsum(clrs[cond], supports, transforms={'balanced': lambda p: p['count'] * p['weight1'] * p['weight2']}, map=pool.map)
			# cis_exp = pd.concat([tables[support] for support in supports], keys=[support[0] for support in supports], names=['chrom'])
			cis_exp['balanced.avg'] = cis_exp['balanced.sum'] / cis_exp['n_valid']
			# cis_exp.reset_index(level=['chrom', 'diag'], inplace=True)
			# cis_exp = cooltools.expected.cis_expected(clrs[cond], supports)

			snipper = snipping.ObsExpSnipper(clrs[cond], cis_exp, regions=supports)
			stack = snipping.pileup(windows, snipper.select, snipper.snip, map=pool.map)
			pile_matrix = np.log2(np.nanmean(stack, axis=2))
			piles[cond] = pile_matrix
			xmin = np.nanmin(pile_matrix)
			xmax = np.nanmax(pile_matrix)
			vmin = xmin if abs(vmin) < abs(xmin) else vmin
			vmax = xmax if abs(vmax) < abs(xmax) else vmax

			x, y, z = stack.shape
			for z_i in range(z):
				if center > 0:
					center_mean = np.nanmean(stack[x // 2 - center + 1: x // 2 + center, y // 2 - center + 1: y // 2 + center, z_i])
					rest_matrix = stack[:, :, z_i]
					rest_matrix[x // 2 - center + 1: x // 2 + center, y // 2 - center + 1: y // 2 + center] = np.nan
					rest_mean = np.nanmean(rest_matrix)
					# sites_enrich[z_i][cond] = center_mean / rest_mean
					sites_enrich[z_i][cond] = center_mean
				elif center == 0:
					sites_enrich[z_i][cond] = get_inter(stack[:, :, z_i])
					inter_matrix[z_i][cond] = np.nan_to_num(np.array(stack[x // 2 + 1: x,  0: x // 2, z_i]).flatten())

	for z_i in sorted(inter_matrix.keys()):
		for i, cond_i in enumerate(conditions):
			for j in range(i+1, len(conditions)):
				cond_j = conditions[j]
				t, p = stats.ttest_rel(inter_matrix[z_i][cond_i], inter_matrix[z_i][cond_j])
				log2FC = np.log2(sites_enrich[z_i][cond_j] / sites_enrich[z_i][cond_i])
				sites_enrich[z_i][f'{cond_j}_{cond_i}_P'] = p
				sites_enrich[z_i][f'{cond_j}_{cond_i}_log2FC'] = log2FC

	if abs(vmin) > abs(vmax):
		vmax = -1 * vmin
	else:
		vmin = -1 * vmax
	gsl = gs.subgridspec(nrows=1, ncols=len(conditions) + 1, width_ratios=[20] * len(conditions) + [1])
	opts = dict(vmin=vmin, vmax=vmax, extent=[-flank // 1000, flank // 1000, -flank // 1000, flank // 1000], cmap='coolwarm')

	for i, cond in enumerate(conditions):
		matrix = piles[cond]
		ax = plt.subplot(gsl[i])
		img = ax.matshow(matrix, **opts)
		if center > 1:
			center_mean = np.nanmean(matrix[x // 2 - center + 1: x // 2 + center, y // 2 - center + 1: y // 2 + center])
			matrix[x // 2 - center + 1: x // 2 + center, y // 2 - center + 1: y // 2 + center] = np.nan
			rest_mean = np.nanmean(matrix)
			# enrichment = round(center_mean / rest_mean, 3)
			enrichment = round(center_mean, 3)
			ax.text(x // 10, y - y // 10, str(enrichment))
		ax.xaxis.tick_bottom()
		if i > 0:
			ax.yaxis.set_visible(False)
		if row_i == 0:
			ax.set_title(cond)
		ax.set_ylabel(name)
		# ax.set_xlabel(f'enrichment:{enrichment}')

	ax = plt.subplot(gsl[len(conditions)])
	plt.colorbar(img, cax=ax)
	return sites_enrich


def bedpileup(bedfile, bedname, clrs, conditions, gs, outdir, flank=250000, res=5000, row_i=0):
	from cooltools import snipping

	# mm10 = bioframe.fetch_chromsizes('mm10')
	chromsizes = bioframe.fetch_chromsizes('mm10', as_bed=True)
	supports = chromsizes.set_index("chrom").reset_index()
	supports = bioframe.parse_regions(supports) 

	sites = pd.read_csv(bedfile, sep="\t", names=['chrom', 'start', 'end'])
	# sites = bioframe.read_table(bedfile, schema='bed')
	windows = snipping.make_bin_aligned_windows(res, sites['chrom'], (sites['start'] + sites['end']) // 2, flank_bp=flank)
	windows = snipping.assign_regions(windows, supports)
	windows = windows.dropna()

	sites_enrich = pileup(clrs, conditions, windows, supports, 0, flank, gs, outdir, bedname, row_i)

	sites_enrich_df = pd.DataFrame.from_dict(sites_enrich, orient="index")
	sites_enrich_df = pd.concat([sites, sites_enrich_df], axis=1)
	sites_enrich_df.to_csv(f'{outdir}/{bedname}.enrich.csv', sep="\t", index=False)



def bedpepileup(bedpefile, bpname, clrs, conditions, gs, outdir, center=3, flank=100000, res=5000, row_i=0):
	from cooltools import snipping
	chromsizes = bioframe.fetch_chromsizes('mm10', as_bed=True)
	supports = chromsizes.set_index("chrom").reset_index()
	supports = bioframe.parse_regions(supports) 

	sites = pd.read_csv(bedpefile, sep="\t", names=['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2'])
	snippet_flank = flank

	windows1 = snipping.make_bin_aligned_windows(res, sites['chrom1'], (sites['start1'] + sites['end1']) // 2, flank_bp=snippet_flank)
	windows2 = snipping.make_bin_aligned_windows(res, sites['chrom2'], (sites['start2'] + sites['end2']) // 2, flank_bp=snippet_flank)

	windows = pd.merge(windows1, windows2, left_index=True, right_index=True, suffixes=('1', '2'))
	windows = snipping.assign_regions(windows, supports)
	windows = windows.dropna()

	sites_enrich = pileup(clrs, conditions, windows, supports, center, flank, gs, outdir, bpname, row_i)
	sites_enrich_df = pd.DataFrame.from_dict(sites_enrich, orient="index")
	sites_enrich_df = pd.concat([sites, sites_enrich_df], axis=1)
	sites_enrich_df.to_csv(f'{outdir}/{bpname}.enrich.csv', sep="\t", index=False)


def bedpair(target, bedfile, bedname, clrs, conditions, gs, outdir, mindist=100000, maxdist=1000000, center=3, flank=100000, res=5000, row_i=0):
	from cooltools import snipping
	import itertools
	chromsizes = bioframe.fetch_chromsizes('mm10', as_bed=True)
	supports = chromsizes.set_index("chrom").reset_index()
	supports = bioframe.parse_regions(supports) 
	snippet_flank = flank

	sites1 = pd.read_csv(target, sep="\t", names=['chrom', 'start', 'end'])
	sites2 = pd.read_csv(bedfile, sep="\t", names=['chrom', 'start', 'end'])
	pair_bed = list(itertools.product(sites1.values.tolist(), sites2.values.tolist()))

	filter_pair = [x[0] + x[1] for x in pair_bed if x[0][0] == x[1][0] and abs(x[0][1] - x[1][1]) > mindist and abs(x[0][1] - x[1][1]) < maxdist]
	sites = pd.DataFrame.from_records(filter_pair, columns=['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2'])
	windows1 = snipping.make_bin_aligned_windows(res, sites['chrom1'], (sites['start1'] + sites['end1']) // 2, flank_bp=snippet_flank)
	windows2 = snipping.make_bin_aligned_windows(res, sites['chrom2'], (sites['start2'] + sites['end2']) // 2, flank_bp=snippet_flank)

	windows = pd.merge(windows1, windows2, left_index=True, right_index=True, suffixes=('1', '2'))
	windows = snipping.assign_regions(windows, supports)
	windows = windows.dropna()

	sites_enrich = pileup(clrs, conditions, windows, supports, center, flank, gs, outdir, bedname, row_i)
	sites_enrich_df = pd.DataFrame.from_dict(sites_enrich, orient="index")
	sites_enrich_df = pd.concat([sites, sites_enrich_df], axis=1)
	sites_enrich_df.to_csv(f'{outdir}/{bedname}.enrich.csv', sep="\t", index=False)


def main():
	parser = argparse.ArgumentParser(description="plot pileups and export enrichment of each region")
	parser.add_argument("--bed", nargs="+")
	parser.add_argument("--bname", nargs="+")
	parser.add_argument("--bedpe", nargs="+")
	parser.add_argument("--bpname", nargs="+")
	parser.add_argument("--bed2", help="Another bed file form paired peaks")
	parser.add_argument("--cool", required=True, nargs="+")
	parser.add_argument("--cname", nargs="+")
	parser.add_argument("--flank", type=int, default=250000, help="The flank region size")
	parser.add_argument("--res", type=int, default=5000, help="Resolution of the cool file")
	parser.add_argument("--center", type=int, default=1, help="Size of the center to calculate enrichment")
	parser.add_argument("--outdir", default=".")
	parser.add_argument("--prefix", default="pileup", help="prefix of the figure file")
	args = parser.parse_args()

	cname = list()
	clrs = dict()
	for i, coolfile in enumerate(args.cool):
		if args.cname:
			basename = args.cname[i]
			cname = args.cname 
		else:
			basename = re.sub("\.\S+", "", os.path.basename(coolfile))
			cname.append(basename)
		clrs[basename] = cooler.Cooler(coolfile)

	cols = len(cname)
	rows = 0
	if args.bed:
		rows += len(args.bed)
	if args.bedpe:
		rows += len(args.bedpe)
	fig = plt.figure(figsize=[cols*2, rows*2])
	gs = fig.add_gridspec(rows, 1, figure=fig, hspace=0.1)

	row_i = 0
	if args.bed and args.bed2:
		bnames = list()
		if args.bname:
			bnames = args.bname
		for i, bedfile in enumerate(args.bed):
			bedname = bnames[i] if len(bnames) > 0 else re.sub("\.bed", "", os.path.basename(bedfile))
			bedpair(target=args.bed2, bedfile=bedfile, bedname=bedname, clrs=clrs, conditions=cname, gs=gs[row_i], outdir=args.outdir, row_i=row_i)
			row_i += 1
	elif args.bed:
		bnames = list()
		if args.bname:
			bnames = args.bname
		for i, bedfile in enumerate(args.bed):
			bedname = bnames[i] if len(bnames) > 0 else re.sub("\.bed", "", os.path.basename(bedfile))
			bedpileup(bedfile, bedname, clrs, cname, gs[row_i], args.outdir, flank=args.flank, row_i=row_i)
			row_i += 1

	if args.bedpe:
		bpnames = list()
		if args.bpname:
			bpnames = args.bpname
		for bedpefile in args.bedpe:
			bpname = bpnames[i] if len(bpnames) > 0 else re.sub("\.bedpe", "", os.path.basename(bedpefile))
			bedpepileup(bedpefile, bpname, clrs, cname, gs[row_i], args.outdir, row_i=row_i)
			row_i += 1
	plt.savefig(f'{args.outdir}/{args.prefix}.pdf')


if __name__ == '__main__':
	main()
