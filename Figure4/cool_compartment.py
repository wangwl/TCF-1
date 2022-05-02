from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
mpl.style.use('seaborn-white')

import multiprocess as mp
import numpy as np
import pandas as pd
import bioframe
import cooltools
import cooler
import argparse
import os
import re


def scatterplot(eigs, tags, outdir):
	chromsizes = bioframe.fetch_chromsizes('mm10')
	chromosomes = list(chromsizes.index)

	conditions = tags
	long_names = dict(zip(tags, tags))
	pal = sns.color_palette('colorblind')
	colors = {tags[i] : pal[i] for i in range(len(tags))}


	from scipy.stats import rankdata

	ncols = len(tags) - 1
	nrows = 4
	gs = plt.GridSpec(nrows=nrows, ncols=ncols)
	plt.figure(figsize=(ncols*4, nrows*3))

	condx = tags[0]

	for i, condy in enumerate(tags[1:]):
		plt.subplot(gs[0,i])
		lo, hi = -2 , 2
		plt.hexbin(
		    eigs[condx]['E1'],
		    eigs[condy]['E1'],
		    vmax=50,
		)
		plt.xlabel('E1 ' + long_names[condx])
		plt.ylabel('E1 ' + long_names[condy])
		plt.gca().set_aspect(1)
		plt.xlim(lo, hi)
		plt.ylim(lo, hi)
		plt.axvline(0, c='b', lw=0.5, ls='--')
		plt.axhline(0, c='b', lw=0.5, ls='--')
		plt.plot([lo, hi], [lo, hi], c='b', lw=0.5, ls='--')
		plt.colorbar(shrink=0.6)


		plt.subplot(gs[1,i])
		mask = eigs[condx]['E1'].notnull() & eigs[condy]['E1'].notnull() 
		vx = eigs[condx]['E1'].loc[mask].values
		vy = eigs[condy]['E1'].loc[mask].values
		lo, hi = 0 , len(vx)

		plt.hexbin(
		    rankdata(vx),
		    rankdata(vy),
		    vmax=20,
		)
		plt.xlabel('E1 rank ' + long_names[condx])
		plt.ylabel('E1 rank ' + long_names[condy])
		plt.gca().set_aspect(1)
		plt.xlim(lo, hi)
		plt.ylim(lo, hi)
		plt.plot([lo, hi], [lo, hi], c='b', lw=0.5, ls='--')
		plt.colorbar(shrink=0.6)

	condx = tags[-1]
	for i, condy in enumerate(tags[:-1]):
		plt.subplot(gs[2,i])
		lo, hi = -2 , 2
		plt.hexbin(
		    eigs[condx]['E1'],
		    eigs[condy]['E1'],
		    vmax=50,
		)
		plt.xlabel('E1 ' + long_names[condx])
		plt.ylabel('E1 ' + long_names[condy])
		plt.gca().set_aspect(1)
		plt.xlim(lo, hi)
		plt.ylim(lo, hi)
		plt.axvline(0, c='b', lw=0.5, ls='--')
		plt.axhline(0, c='b', lw=0.5, ls='--')
		plt.plot([lo, hi], [lo, hi], c='b', lw=0.5, ls='--')
		plt.colorbar(shrink=0.6)


		plt.subplot(gs[3,i])
		mask = eigs[condx]['E1'].notnull() & eigs[condy]['E1'].notnull() 
		vx = eigs[condx]['E1'].loc[mask].values
		vy = eigs[condy]['E1'].loc[mask].values
		lo, hi = 0 , len(vx)

		plt.hexbin(
		    rankdata(vx),
		    rankdata(vy),
		    vmax=20,
		)
		plt.xlabel('E1 rank ' + long_names[condx])
		plt.ylabel('E1 rank ' + long_names[condy])
		plt.gca().set_aspect(1)
		plt.xlim(lo, hi)
		plt.ylim(lo, hi)
		plt.plot([lo, hi], [lo, hi], c='b', lw=0.5, ls='--')
		plt.colorbar(shrink=0.6)

	plt.savefig(f"{outdir}/PC1.scatter.pdf")


def plot_chr(eigs, tags, gg, outdir):
	PC1_df = eigs[tags[0]][['chrom', 'start']]
	for tag in tags:
		PC1_df[tag] = eigs[tag]['E1']
	PC1_df = PC1_df.groupby('chrom')

	chromsizes = bioframe.fetch_chromsizes(gg)
	chromsizes = chromsizes.drop("chrM")
	chromosomes = list(chromsizes.index)

	pal = sns.color_palette('colorblind', n_colors=len(tags))
	colors = {tags[i]: pal[i] for i in range(len(tags))}

	fig_cols = 2
	fig_rows = len(chromosomes) // fig_cols + 1
	gs = plt.GridSpec(nrows=fig_rows, ncols=fig_cols)
	fig = plt.figure(figsize=(10 * fig_cols, 3 * fig_rows))
	i = 0

	for chrom, chrom_pc1 in PC1_df:
		if chrom not in chromosomes:
			continue
		coordinates = [x/1000000 for x in chrom_pc1['start']]
		ax = plt.subplot(gs[i])
		for tag in tags:
			PC1 = chrom_pc1[tag]
			ax.plot(coordinates, PC1, '-', color=colors[tag], label=tag, lw=0.1)
			# ax.set_xlabel(f"Position on {chrom}")
			ax.set_ylabel("PC1 value")
			ax.set_title(f"{chrom}")
		ax.legend()
		i += 1
	plt.savefig(f'{outdir}/chromPC1.pdf')



def computePC1(cools, tags, binsize, fa, gg, outdir):
	coolfiles = [f'{x}::/resolutions/{binsize}' if x.endswith("mcool") else x for x in cools]

	cooler_paths = dict(zip(tags, coolfiles))
	clrs = {
		cond: cooler.Cooler(cooler_paths[cond]) for cond in tags
	}

	chromsizes = bioframe.fetch_chromsizes(gg)
	chromsizes = chromsizes.drop("chrM")
	chromosomes = list(chromsizes.index)
	bins = cooler.binnify(chromsizes, binsize)
	fasta_records = bioframe.load_fasta(fa)
	bins = bioframe.frac_gc(bins, fasta_records)

	from cooltools.eigdecomp import cooler_cis_eig
	lam = {}
	eigs = {}

	for cond in tags:
		if os.path.exists(f'{outdir}/{cond}.{binsize//1000}kb.eigs.cis.vecs.txt'):
			eigs[cond] = pd.read_csv(f'{outdir}/{cond}.{binsize//1000}kb.eigs.cis.vecs.txt', sep="\t")
			continue
		lam[cond], eigs[cond] = cooler_cis_eig(clrs[cond], bins, n_eigs=3, phasing_track_col='GC', sort_metric='var_explained')
		# Save text files
		lam[cond].to_csv(f'{outdir}/{cond}.{binsize//1000}kb.eigs.cis.lam.txt', sep='\t')
		eigs[cond].to_csv(f'{outdir}/{cond}.{binsize//1000}kb.eigs.cis.vecs.txt', sep='\t', index=False)

		# Save bigwig track
		bioframe.to_bigwig(eigs[cond], chromsizes, f'{outdir}/{cond}.{binsize//1000}kb.eigs.cis.vecs.E1.bw', 'E1')

	# scatterplot(eigs, tags, outdir)
	plot_chr(eigs, tags, gg, outdir)


def interaction_ratio(S, C, segments, square):
	"""
	Parameters
	----------
	S, C : 2D arrays, square, same shape
		Saddle sums and counts, respectively
	segments: 1D array segments of the matrix to calculate interaction
	Returns
	-------
	1D array
	Ratios of cumulative corner interaction scores, where the saddle data is
	grouped over the AA+BB corners and AB+BA corners with increasing extent.
	"""
	m, n = S.shape
	if m != n:
		raise ValueError("`saddledata` should be square.")

	if sum(segments) != 1:
		raise ValueError("segments should summed up to 1.")

	ratios = np.zeros(n)
	for k in range(1, n):
		intra_sum = np.nansum(S[k - 1: k + square - 1, k - 1: k + square - 1])
		intra_count = np.nansum(C[k - 1: k + square - 1, k - 1: k + square - 1])
		intra = intra_sum / intra_count

		interB_sum = np.nansum(S[0: square, k - 1: k + square - 1])
		interB_count = np.nansum(C[0: square, k - 1: k + square - 1])
		interB = interB_sum / interB_count

		interA_sum = np.nansum(S[k - 1: k + square - 1, n - square: n])
		interA_count = np.nansum(C[k - 1: k + square - 1, n - square: n])
		interA = interA_sum / interA_count

		ratios[k] = intra / max(interA, interB)

	i_start = 0
	interactions = dict()
	for i in range(len(segments)):
		i_size = int(segments[i] * n)
		i_intra_sum = np.nansum(S[i_start: i_start + i_size, i_start: i_start + i_size])
		i_intra_count = np.nansum(C[i_start: i_start + i_size, i_start: i_start + i_size])
		i_intra = i_intra_sum / i_intra_count
		interactions[i] = i_intra

		j_start = i_start + i_size
		for j in range(i + 1, len(segments)):
			j_size = int(segments[j] * n)
			inter_sum = np.nansum(S[i_start: i_start + i_size, j_start: j_start + j_size])
			inter_count = np.nansum(C[i_start: i_start + i_size, j_start: j_start + j_size])
			inter = inter_sum / inter_count
			interactions[f'{i}-{j}'] = inter
			j_start += j_size
		i_start += i_size

	seg_len = len(segments) - 1
	interactions['strength'] = (interactions[0] + interactions[seg_len]) / (interactions[f'0-{seg_len}'] * 2)

	return ratios, interactions


def saddleplot(cools, tags, outdir, binsize):
	from cooltools import saddle
	from cooltools.expected import diagsum

	QUANTILE_BINNING = True

	binedges = {}
	digitized = {}
	hist = {}
	sums = {}
	counts = {}
	saddledata = {}

	coolfiles = [f'{x}::/resolutions/{binsize}' if x.endswith("mcool") else x for x in cools]
	cooler_paths = dict(zip(tags, coolfiles))
	clrs = {
		cond: cooler.Cooler(cooler_paths[cond]) for cond in tags
	}
	pal = sns.color_palette('colorblind', n_colors=len(tags))
	colors = {tags[i]: pal[i] for i in range(len(tags))}

	gs = plt.GridSpec(nrows=1, ncols=len(tags))
	fig = plt.figure(figsize=(6 * len(tags), 6))
	histbins = 50

	chromsizes = bioframe.fetch_chromsizes('mm10', as_bed=True)
	supports = chromsizes.set_index("chrom").reset_index()
	supports = bioframe.parse_regions(supports) 

	for i, cond in enumerate(tags):
		if os.path.exists(f'{outdir}/{cond}.{binsize//1000}kb.expected.cis.tsv'):
			exp = pd.read_csv(f'{outdir}/{cond}.{binsize//1000}kb.expected.cis.tsv', sep="\t")
		else:
			exp = diagsum(clrs[cond], supports, transforms={'balanced': lambda p: p['count'] * p['weight1'] * p['weight2']})
			# exp = pd.concat([tables[support] for support in supports], keys=[support[0] for support in supports], names=['chrom'])
			exp['balanced.avg'] = exp['balanced.sum'] / exp['n_valid']
			exp.to_csv(f'{outdir}/{cond}.{binsize//1000}kb.expected.cis.tsv', sep="\t")

		eig = pd.read_csv(f'{outdir}/{cond}.{binsize//1000}kb.eigs.cis.vecs.txt', sep="\t")
		# Determine how to bin the range of the E1 signal
		if QUANTILE_BINNING:
			q_binedges = np.linspace(0, 1, histbins)
			binedges[cond] = saddle.quantile(eig['E1'], q_binedges)
		else:
			qlo, qhi = saddle.quantile(eig['E1'], [0.02, 0.98])  # trim outliers
			binedges[cond] = np.linspace(qlo, qhi, histbins)

		# Digitize the signal into integers
		digitized[cond], hist[cond] = saddle.digitize_track(
			binedges[cond],
			track=(eig, 'E1'))
		# Construct a function that fetches and calculates observed/expected
		getmatrix = saddle.make_cis_obsexp_fetcher(clrs[cond], (exp, 'balanced.avg'))

		# Build the saddle histogram
		sums[cond], counts[cond] = saddle.make_saddle(
			getmatrix,
			binedges[cond],
			(digitized[cond], 'E1.d'),
			contact_type='cis')
		saddledata[cond] = sums[cond] / counts[cond]

		# Make the saddle plot
		g = saddle.saddleplot(
			q_binedges if QUANTILE_BINNING else binedges[cond],
			hist[cond],
			saddledata[cond],
			# np.log10(saddledata[cond]),
			color=colors[cond],
			# heatmap_kws={'vmin': -0.5, 'vmax': 0.5},
			fig=fig, subplot_spec=gs[i])
	plt.savefig(f'{outdir}/All.{binsize//1000}kb.saddle.pdf')

	strength = {
		cond: saddle.saddle_strength(sums[cond], counts[cond]) for cond in tags
	}

	gs = plt.GridSpec(nrows=1, ncols=3)
	plt.figure(figsize=(21, 6))

	plt.subplot(gs[0])
	x = np.arange(histbins + 2)
	for cond in tags:
		plt.step(x[:-1], strength[cond], where='pre', color=colors[cond], label=cond)

	plt.legend()
	plt.xlabel('extent')
	plt.ylabel('(AA + BB) / (AB + BA)')
	plt.title('saddle strength profile')
	plt.axhline(0, c='grey', ls='--', lw=1)
	plt.xlim(0, len(x) // 2)

	plt.subplot(gs[1])
	plt.step(x[:-1], strength[tags[-1]] / strength[tags[0]], where='pre', c='k')
	plt.axhline(1, c='grey', ls='--', lw=1)
	plt.xlim(0, len(x) // 2)
	plt.xlabel('extent')
	plt.ylabel('enrichment')
	plt.title(f'{tags[-1]} / {tags[0]}')

	enrich1 = dict()
	enrich2 = dict()
	segments = [0.3, 0.4, 0.3]
	for cond in tags:
		ratios, interactions = interaction_ratio(sums[cond], counts[cond], segments, 4)
		enrich1[cond] = ratios
		enrich2[cond] = interactions

	plt.subplot(gs[2])
	x = np.arange(histbins + 2)
	for cond in tags:
		plt.step(x[:-1], enrich1[cond], where='pre', color=colors[cond], label=cond)

	plt.legend()
	plt.xlabel('Quantile')
	plt.ylabel('Local enrichment')
	plt.title('Local quantile enrichment profile')
	plt.axhline(0, c='grey', ls='--', lw=1)
	plt.xlim(0, len(x))
	plt.savefig(f'{outdir}/All.{binsize//1000}kb.saddle.strength.pdf')

	fig = plt.figure(figsize=(6, 20))
	enrich2_df = pd.DataFrame.from_dict(enrich2, orient='index')
	enrich2_df.to_csv(f'{outdir}/Segments.interaction.csv')
	axes = enrich2_df.plot.bar(rot=0, subplots=True)
	plt.savefig(f'{outdir}/Segments.interaction.pdf')


def main():
	parser = argparse.ArgumentParser(description="Do compartment analysis for group of cool files")
	parser.add_argument("--cools", nargs="+", required=True, help="cool/mcool files")
	parser.add_argument("--tags", nargs="+", help="Names for the cool files")
	parser.add_argument("--bin", type=int, default=20000)
	parser.add_argument("-g", default="mm10", help="The genome version (default mm10)")
	parser.add_argument("--ref", default="/mnt/data0/wenliang/database/iGenomes/Mus_musculus/Ensembl/GRCm38/Sequence/WholeGenomeFasta/genome.chr.fa", help="Fasta file of the reference genome")
	parser.add_argument("--outdir", default=".", help="Directory to put the results in")
	args = parser.parse_args()

	if args.tags:
		tags = args.tags
	else:
		tags = []
		for coolfile in args.cools:
			filename = os.path.basename(coolfile)
			tags.append(re.sub("\.\S+", "", filename))

	outdir = args.outdir
	if not os.path.exists(f'{outdir}/PC1.scatter.pdf'):
		computePC1(args.cools, tags, args.bin, args.ref, args.g, outdir)

	if not os.path.exists(f'{outdir}/All.{args.bin//1000}kb.saddle.pdf'):
		saddleplot(args.cools, tags, outdir, args.bin)


if __name__ == '__main__':
	main()
