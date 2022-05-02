import os
import re
import sys
import pandas as pd
import numpy as np
from sklearn.neighbors import KernelDensity
from sklearn.preprocessing import normalize
from collections import defaultdict
import pybedtools
import bbi

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle

import argparse


def data_density(data, bins, xmin, xmax):
	data = np.array(data)
	X = data[:, np.newaxis]
	X_plot = np.linspace(xmin, xmax, bins)[:,np.newaxis]
	kde = KernelDensity(kernel='gaussian', bandwidth=0.75).fit(X)
	log_dens = kde.score_samples(X_plot)
	return np.exp(log_dens)


def group_boundaries(boundary):
	data_df = pd.read_csv(boundary)
	data_df[["DN2", "DN3", "ETP", "CLP", "DN4", "DP"]] = data_df[["DN2", "DN3", "ETP", "CLP", "DN4", "DP"]].astype(int)
	data_df[["DPPC1", "CLPPC1", "ETPPC1", "DN2PC1", "DN3PC1", "DN4PC1", "strength"]] = \
		data_df[["DPPC1", "CLPPC1", "ETPPC1", "DN2PC1", "DN3PC1", "DN4PC1", "strength"]].astype(float)

	decode = {
	'000001': "A", '000010': 'A', '000100': 'A', '100000': 'C', '010000': 'C', 
'001000': 'C',\
	'011000': 'D', '110000': 'D', '101000': 'D', '000011': 'B', '000101': 'B', '000110': 'B', '000111': 'B',\
	'111000': 'E', '111100': 'F', '111110': 'G', '111111': 'H'
	}

	groupdict = defaultdict(list)
	group_df = data_df.groupby(["CLP", "ETP", "DN2", "DN3", "DN4", "DP"])
	for code, intersect_df in group_df:
		code_str = "".join([str(x) for x in code])
		if code_str == '000000':
			continue
		if code_str not in decode:
			group = 'J'
		else:
			group = decode[code_str]
		groupdict[group].append(intersect_df)

	for group in groupdict:
		groupdict[group] = pd.concat(groupdict[group])
	return groupdict


def barcode_plot(boundary, outdir, gs, bins=1000, xmin=-1.5, xmax=1.5, columns=["CLPPC1", "ETPPC1", "DN2PC1", "DN3PC1", "DN4PC1", "DPPC1"]):
	'''
	barcode plot of some quantitative features at target regions
	'''
	gsl = gs.subgridspec(nrows=len(columns) + 1, height_ratios=[1, 2, 2, 2, 2, 2, 2], ncols=9, hspace=0, wspace=0)
	colors = {
		'A': 'blue', 'C': 'green', 'D': '#C7EA46', 'B': '#63C5DA', 'E': 'yellow', 'F': 'orange', 'G': 'red', 'H': 'purple', 'J': 'black'
	}
	group_bounds = group_boundaries(boundary)
	pc_dict = defaultdict(dict)
	all_density = np.empty(0)
	for i, group in enumerate(list("ABCDEFGHJ")):
		group_df = group_bounds[group]
		for j in range(len(columns)):
			j = int(j)
			density_np = data_density(group_df[columns[j]].to_list(), bins, xmin, xmax)
			pc_dict[i][j] = density_np
			all_density = np.append(all_density, density_np)
	axes = list()
	vmin, vmax = np.quantile(all_density, [0.01, 0.99])
	for i, group in enumerate(list("ABCDEFGHJ")):
		ax = ax1 = plt.subplot(gsl[0, i])
		rect = ax.patch
		rect.set_facecolor(colors[group])
		ax.xaxis.set_visible(False)
		ax.yaxis.set_visible(False)
		if i == 8:
			axes.append(ax1)
		for j in range(len(columns)):
			density_np = pc_dict[i][j]
			density_np = (density_np - vmin) / (vmax - vmin)
			barprops = dict(aspect='auto', cmap='coolwarm', interpolation='nearest', vmin=0, vmax=1)
			ax = plt.subplot(gsl[j + 1, i], sharex=ax1)
			im = ax.imshow(density_np.reshape((1, -1)), **barprops)
			ax.axvline(bins // 2, c='black', lw=0.5)
			ax.axvline(bins // 4, c='grey', lw=0.5)
			ax.axvline(bins * 3 // 4, c='grey', lw=0.5)
			# ax.yaxis.set_visible(False)
			ax.set_yticks([])
			ax.xaxis.set_visible(False)
			ax.grid(False)
			if i == 0:
				tag = re.sub("PC1", "", columns[j])
				ax.set_ylabel(tag, fontsize=4, rotation=0, labelpad=15)
			if i == 8:
				axes.append(ax)
	cax, kw = matplotlib.colorbar.make_axes(axes)
	plt.colorbar(im, cax=cax, **kw)


def quanti_center(stacks, vmin, vrange, flank):
	regions, bins = stacks.shape
	avg_bw = np.nanmean(stacks, axis=0)
	center_bins = bins * 100000 // (2 * flank)

	enrichment = list()
	mid = bins // 2
	for record in stacks:
		enrich = 0
		# record = (record - vmin) / vrange
		# record = [x if x>0 else 0 for x in record]
		# record = [x if x<1 else 1 for x in record]
		center_bw = record[mid - center_bins // 2: mid + center_bins // 2]
		enrich = np.mean(center_bw)
		# for i in range(mid):
			# distance = i + 1
			# enrich += (record[mid - 1 - i] + record[mid + i]) / distance
		enrichment.append(enrich)

	avg_enrich = 0
	# avg_bw = (avg_bw - vmin) / vrange
	# avg_bw = [x if x>0 else 0 for x in avg_bw]
	# avg_bw = [x if x<1 else 1 for x in avg_bw]
	avg_center_bw = avg_bw[mid - center_bins // 2: mid + center_bins // 2]
	avg_enrich = np.mean(avg_center_bw)
	# for i in range(mid):
	# 	distance = i + 1
	# 	avg_enrich += (avg_bw[mid - 1 - i] + avg_bw[mid + i]) / distance

	return enrichment, avg_enrich


def bw_plot(group_bounds, bwfiles, tags, gs, nbins=200, flank=250000, cmap="YlGnBu"):
	height_ratios = [2] * len(tags)
	gsl = gs.subgridspec(nrows=len(tags), height_ratios=height_ratios, ncols=9, hspace=0, wspace=0)

	bw_enrich = defaultdict(dict)
	bw_dict = defaultdict(dict)
	ranges = dict()
	for i, bwfile in enumerate(bwfiles):
		tag = tags[i]
		all_avg = np.empty(0)
		for j, group in enumerate(list("ABCDEFGHJ")):
			group_df = group_bounds[group]
			mids = (group_df['start'] + group_df['end']) // 2
			stacks = bbi.stackup(bwfile, group_df['chrom'], mids - flank, mids + flank, bins=nbins)
			avg_bw = np.nanmean(stacks, axis=0)
			avg_bw = np.nan_to_num(avg_bw)
			all_avg = np.append(all_avg, avg_bw)
			bw_dict[tag][group] = stacks
		xmin, xmax = np.quantile(all_avg, [0.02, 0.98])
		ranges[tag] = [xmin, xmax - xmin]

	axes = list()
	vmin = 0
	vmax = 1
	for i, tag in enumerate(tags):
		for j, group in enumerate(list("ABCDEFGHJ")):
			group_df = group_bounds[group]
			stacks = bw_dict[tag][group]
			avg_bw = np.nanmean(stacks, axis=0)
			avg_bw = np.nan_to_num(avg_bw)
			avg_bw = (avg_bw - ranges[tag][0]) / ranges[tag][1]
			avg_bw = [x if x>0 else 0 for x in avg_bw]
			avg_bw = [x if x<1 else 1 for x in avg_bw]
			avg_bw = np.array(avg_bw)
			
			enrichment, avg_enrich = quanti_center(stacks, ranges[tag][0], ranges[tag][1], flank)
			bw_enrich[group][tag] = avg_enrich
			try:
				group_df[tag] = enrichment
			except:
				print(group_df.shape, len(enrichment))
			group_df['bound'] = group

			ax = plt.subplot(gsl[i, j])
			barprops = dict(aspect='auto', cmap=cmap, interpolation='nearest', vmin=vmin, vmax=vmax)
			im = ax.imshow(avg_bw.reshape((1, -1)), **barprops)
			ax.set_yticks([])
			ax.xaxis.set_visible(False)
			# ax.yaxis.set_visible(False)
			ax.grid(False)
			if j == 0:
				ax.set_ylabel(tag, fontsize=4, rotation=0, labelpad=15)
			if j == 8:
				axes.append(ax)
	cax, kw = matplotlib.colorbar.make_axes(axes)
	plt.colorbar(im, cax=cax, **kw)

	return bw_enrich


def scale_TAD_plot(all_df, bwfiles, tags, outdir, chrs_info, nbins=200):
	gs = GridSpec(nrows=3, ncols=len(tags) + 1, height_ratios=[15, 2, 0.5], hspace=0)
	plt.figure(figsize=(3 * len(tags) + 1, 10))

	stacks = defaultdict(list)
	for i, boundary in all_df.iterrows():
		this_bound = (boundary.start + boundary.end) // 2

		if i == 0 or all_df.iloc[i - 1].chrom != boundary.chrom:
			up_bound = chrs_info[boundary.chrom][0]
			down_bound = (all_df.iloc[i + 1].start + all_df.iloc[i + 1].end) // 2
		elif i == all_df.shape[0] - 1 or all_df.iloc[i + 1].chrom != boundary.chrom:
			down_bound = chrs_info[boundary.chrom][1]
			up_bound = (all_df.iloc[i - 1].start + all_df.iloc[i - 1].end) // 2
		else:
			up_bound = (all_df.iloc[i - 1].start + all_df.iloc[i - 1].end) // 2
			down_bound = (all_df.iloc[i + 1].start + all_df.iloc[i + 1].end) // 2

		for i, bwfile in enumerate(bwfiles):
			tag = tags[i]
			left = bbi.fetch(bwfile, boundary.chrom, up_bound, this_bound, bins=nbins)
			right = bbi.fetch(bwfile, boundary.chrom, this_bound, down_bound, bins=nbins)
			tad_cov = np.concatenate([left, right])
			if np.nanmean(left) < np.nanmean(right):
				tad_cov = np.flip(tad_cov)
			stacks[tag].append(tad_cov)

	X = np.array(stacks[tags[0]])
	idx = np.argsort(X.sum(axis=1))
	x = np.linspace(-1, 1, 2 * nbins)
	cmap = plt.cm.get_cmap('coolwarm')
	cmap.set_bad('#777777')
	im_opts = dict(cmap=cmap)

	bounds = all_df['bound'].to_numpy()
	bound_fh = open(f'{outdir}/bound.order', 'w')
	for bound in bounds:
		print(bound, file=bound_fh)
	for i, name in enumerate(stacks):
		# heatmap
		ax = ax1 = plt.subplot(gs[0, i])
		X = np.array(stacks[name])
		img = ax.matshow(X[idx, :], **im_opts, rasterized=True)
		ax.axvline(0, c='grey', lw=0.5)
		ax.grid('off')
		ax.set_aspect('auto')
		ax.set_title(name)
		if i > 0:
			ax.yaxis.set_visible(False)

		# summary
		ax = plt.subplot(gs[1, i], sharex=ax1)
		ax.axhline(0, c='#777777', lw=1, ls='--')
		ax.plot(x, np.nanmean(stacks[name], axis=0), c='k', lw=2)
		ax.set_xlim(-1, 1)
		ax.xaxis.set_visible(False)
		ax.set_ylim(0, 1)
		if i > 0:
			ax.yaxis.set_visible(False)

		# color bar
		cax = plt.subplot(gs[2, i])
		cb = plt.colorbar(img, cax=cax, orientation='horizontal')
		cb.locator = matplotlib.ticker.MaxNLocator(nbins=3)
		cb.update_ticks()
		name_df = pd.DataFrame(X, columns=range(2 * nbins))
		name_df.to_csv(f'{outdir}/{name}.cov.csv')

	plt.savefig(f'{outdir}/boundary.TFs.pdf')


def main():
	parser = argparse.ArgumentParser(description="Plot 3D and 1D features at different boundaries")
	parser.add_argument("boundary", help="The boundary annotation file")
	parser.add_argument("--bed", nargs="+", help="Bed files of features")
	parser.add_argument("--bnames", nargs="+", help="")
	parser.add_argument("--pc1", default=False)
	parser.add_argument("--TFs", nargs="+", help="The TFs bw files")
	parser.add_argument("--TFtags", nargs="+", help="Bigwig file names")
	parser.add_argument("--ATAC", nargs="+", help="ATAC bw files")
	parser.add_argument("--Atags", nargs="+", help="ATAC names")
	parser.add_argument("--histones", nargs="+", help="Histone bw files")
	parser.add_argument("--Htags", nargs="+", help="Histone tags")
	parser.add_argument("--RNAseq", nargs="+", help='RNAseq bw files')
	parser.add_argument("--Rtags", nargs="+", help='RNAseq tags')
	parser.add_argument("--rescale", nargs="+", help="bw files to plot rescaled TADs")
	parser.add_argument("--ReTags", nargs="+", help="tags for rescaled plots")
	parser.add_argument("--chrs", default="/mnt/data0/wenliang/database/mm10/mm10.chrStart.bed", help="Effective chr start and end")
	parser.add_argument("--outdir", default=".")
	args = parser.parse_args()

	group_bounds = group_boundaries(args.boundary)
	enrich_dict = defaultdict(dict)

	rows = 1
	height_ratios = [6]
	if args.TFs:
		TFlen = len(args.TFs)
		height_ratios.append(TFlen)
		rows += 1
	if args.histones:
		Histlen = len(args.histones)
		height_ratios.append(Histlen)
		rows += 1
	if args.ATAC:
		ATAClen = len(args.ATAC)
		height_ratios.append(ATAClen)
		rows += 1
	if args.RNAseq:
		RNAlen = len(args.RNAseq)
		height_ratios.append(RNAlen)
		rows += 1

	fig = plt.figure(figsize=[6, rows])
	gs = fig.add_gridspec(rows, 1, figure=fig, height_ratios=height_ratios, hspace=0.1)

	if args.pc1:
		barcode_plot(args.boundary, args.outdir, gs[0])

	gs_i = 0	
	if args.TFs:
		gs_i += 1
		TFtags = []
		if args.TFtags:
			TFtags = args.TFtags
		else:
			for i, bwfile in enumerate(args.TFs):
				# tag = re.sub("\.bw", "", os.path.basename(bwfile))
				tag = ".".join(os.path.basename(bwfile).split(".")[0:2])
				TFtags.append(tag)
		bw_dict = bw_plot(group_bounds, args.TFs, TFtags, gs[gs_i], cmap='Blues')
		for x in bw_dict:
			enrich_dict[x].update(bw_dict[x])

	if args.histones:
		gs_i += 1
		Htags = []
		if args.Htags:
			Htags = args.Htags
		else:
			for i, bwfile in enumerate(args.histones):
				# tag = re.sub("\.bw", "", os.path.basename(bwfile))
				tag = ".".join(os.path.basename(bwfile).split(".")[0:2])
				Htags.append(tag)
		bw_dict = bw_plot(group_bounds, args.histones, Htags, gs[gs_i], cmap='Purples')
		for x in bw_dict:
			enrich_dict[x].update(bw_dict[x])

	if args.ATAC:
		gs_i += 1
		Atags = []
		if args.Atags:
			Atags = args.Atags
		else:
			for bwfile in args.ATAC:
				# tag = re.sub("\.bw", "", os.path.basename(bwfile))
				tag = ".".join(os.path.basename(bwfile).split(".")[0:2])
				Atags.append(tag)
		bw_dict = bw_plot(group_bounds, args.ATAC, Atags, gs[gs_i], cmap='Greens')
		for x in bw_dict:
			enrich_dict[x].update(bw_dict[x])

	if args.RNAseq:
		gs_i += 1
		Rtags = []
		if args.Rtags:
			Rtags = args.Rtags
		else:
			for bwfile in args.RNAseq:
				# tag = re.sub("\.bw", "", os.path.basename(bwfile))
				tag = ".".join(os.path.basename(bwfile).split(".")[0:2])
				Atags.append(tag)
		bw_dict = bw_plot(group_bounds, args.RNAseq, Rtags, gs[gs_i], cmap='Oranges')
		for x in bw_dict:
			enrich_dict[x].update(bw_dict[x])

	enrich_df = pd.DataFrame.from_dict(enrich_dict, orient='index')
	enrich_df.index.name = 'bound'
	all_boundaries = list()
	for group in group_bounds:
		group_df = group_bounds[group]
		all_boundaries.append(group_df)
	all_df = pd.concat(all_boundaries)
	all_df = all_df.sort_values(by=['chrom', 'start'])
	all_df = all_df.reset_index(drop=True)
	if gs_i > 0:
		plt.savefig(f'{args.outdir}/integrate_plot.pdf')
		enrich_df.to_csv(f'{args.outdir}/bound_enrich.csv')
		all_df.to_csv(f'{args.outdir}/boundary_enrich.csv', index=False)

	if args.rescale:
		ReTags = []
		if args.ReTags:
			ReTags = args.ReTags
		else:
			for i, bwfile in enumerate(args.rescale):
				tag = ".".join(os.path.basename(bwfile).split(".")[0:2])
				ReTags.append(tag)
		chrs_df = pd.read_csv(args.chrs, sep="\t", names=['chrom', 'start', 'end'])
		chrs_info = dict(zip(chrs_df.chrom, zip(chrs_df.start, chrs_df.end)))
		scale_TAD_plot(all_df, args.rescale, ReTags, args.outdir, chrs_info)


if __name__ == '__main__':
	main()
