"""Calculate HITindex metrics, classifications, and PSI tables."""

# This project uses the HITindex library (https://github.com/thepailab/HITindex)
# Please refer to the repository for more details.

import math
import random
import subprocess
import collections
import gzip
import os
import numpy as np
import scipy
from scipy import optimize
import pandas as pd
from util.common.write_logs import log_message

def _run_intersectbed_to_gzip(out_path, command):
	with gzip.open(out_path, "wb") as out_fh:
		proc = subprocess.Popen(command, stdout=subprocess.PIPE)
		try:
			for chunk in iter(lambda: proc.stdout.read(1024 * 1024), b""):
				if not chunk:
					break
				out_fh.write(chunk)
		finally:
			if proc.stdout is not None:
				proc.stdout.close()
		return_code = proc.wait()
		if return_code != 0:
			raise subprocess.CalledProcessError(return_code, command)

def _run_intersectbed_to_plain(out_path, command):
	with open(out_path, "wb") as out_fh:
		subprocess.run(command, check=True, stdout=out_fh)

def _concat_plain_beds_to_gzip(out_path, input_paths):
	with gzip.open(out_path, "wb") as out_fh:
		for path in input_paths:
			with open(path, "rb") as in_fh:
				while True:
					chunk = in_fh.read(1024 * 1024)
					if not chunk:
						break
					out_fh.write(chunk)

def getExons(bedname):
	"""Load metaexon BED records into the nested HITindex exon dictionary."""
	beddict = collections.defaultdict(lambda:collections.defaultdict(dict))
	exoncount = 0
	if bedname.endswith(".gz"):
		bed_context = gzip.open(bedname, "rt")
	else:
		bed_context = open(bedname, "r", encoding="utf-8")
	with bed_context as bedfile:
		for rawline in bedfile:
			if not rawline.strip():
				continue
			bedline = rawline.strip().split('\t')
			if len(bedline) < 6:
				continue
			name_parts = bedline[3].split(';')
			if len(name_parts) < 2:
				continue
			gene = name_parts[1]
			exon = name_parts[0]
			exoncount += 1
			info = ';'.join(name_parts[2:7])
			start = int(bedline[1])
			end = int(bedline[2])
			beddict[gene][exon]['info'] = info
			beddict[gene][exon]['start'] = start
			beddict[gene][exon]['end'] = end
			beddict[gene][exon]['strand'] = bedline[5]
	log_message("[INFO]", f"{exoncount} exons processed in {len(beddict)} genes.", color="info")
	return(beddict)

def getJuncBed(bed, juncbam, readtype, readstrand):
	"""Intersect junction-read BAMs with exon BED intervals for HITindex metrics."""
	out_gz = juncbam + "_exonjuncs.bed.gz"
	if readstrand == 'none':
		# single bam file, no strand
		log_message("[INFO]", "Read layout: unstranded reads.", color="info")
		command = ["intersectBed", "-abam", juncbam + "_read.bam", "-b", bed, "-bed", "-wo", "-split"]
		_run_intersectbed_to_gzip(out_gz, command)
	elif readtype == 'single' and readstrand == 'r':
		# single bam file, read1=-S
		log_message("[INFO]", "Read layout: fr-firststrand single-end reads.", color="info")
		command = ["intersectBed", "-abam", juncbam + "_read.bam", "-b", bed, "-bed", "-wo", "-S", "-split"]
		_run_intersectbed_to_gzip(out_gz, command)
	elif readtype == 'single' and readstrand == 'f':
		# single bam file, read1=-s
		log_message("[INFO]", "Read layout: fr-secondstrand single-end reads.", color="info")
		command = ["intersectBed", "-abam", juncbam + "_read.bam", "-b", bed, "-bed", "-wo", "-s", "-split"]
		_run_intersectbed_to_gzip(out_gz, command)
	elif readtype == 'paired' and readstrand == 'rf':
		# two bams, read1=-S, read2=-s-
		log_message("[INFO]", "Read layout: fr-firststrand paired-end reads.", color="info")
		read1_path = juncbam + "_read1.bed"
		read2_path = juncbam + "_read2.bed"
		try:
			_run_intersectbed_to_plain(read1_path, ["intersectBed", "-abam", juncbam + "_read1.bam", "-b", bed, "-bed", "-wo", "-S", "-split"])
			_run_intersectbed_to_plain(read2_path, ["intersectBed", "-abam", juncbam + "_read2.bam", "-b", bed, "-bed", "-wo", "-s", "-split"])
			_concat_plain_beds_to_gzip(out_gz, [read1_path, read2_path])
		finally:
			for path in (read1_path, read2_path):
				if os.path.exists(path):
					os.remove(path)
	elif readtype == 'paired' and readstrand == 'fr':
		# two bams, read1=-s, read2=-S
		log_message("[INFO]", "Read layout: fr-secondstrand paired-end reads.", color="info")
		read1_path = juncbam + "_read1.bed"
		read2_path = juncbam + "_read2.bed"
		try:
			_run_intersectbed_to_plain(read1_path, ["intersectBed", "-abam", juncbam + "_read1.bam", "-b", bed, "-bed", "-wo", "-s", "-split"])
			_run_intersectbed_to_plain(read2_path, ["intersectBed", "-abam", juncbam + "_read2.bam", "-b", bed, "-bed", "-wo", "-S", "-split"])
			_concat_plain_beds_to_gzip(out_gz, [read1_path, read2_path])
		finally:
			for path in (read1_path, read2_path):
				if os.path.exists(path):
					os.remove(path)
	else:
		raise RuntimeError(f"Unsupported read layout for junction extraction: readtype={readtype}, strand={readstrand}")

def getExonReads(juncbam, overlap):
	"""Extract junction read start/end support per exon from intersectBed output."""
	readstartdict = collections.defaultdict(dict)
	readenddict = collections.defaultdict(dict)
	with gzip.open(juncbam + "_exonjuncs.bed.gz", 'rt') as bedfile:
		for rawline in bedfile:
			bedline = rawline.strip().split('\t')
			if not bedline or bedline[0] == '':
				continue
			# check that read overlaps both ends of split region by at least overlap
			segments = bedline[10].split(',')
			if int(segments[0]) < int(overlap) or int(segments[1]) < int(overlap): continue
			# check that read overlaps this exon by at least overlap
			if int(bedline[18]) > (int(overlap) - 1):
				exon = bedline[15].split(';')[0] +';'+ bedline[15].split(';')[1]
				readstart = bedline[1]
				readend = bedline[2]
				if bedline[17] == '-':
					readstart = bedline[2]
					readend = bedline[1]
				# beddict[exon_name][read_name][read_position] = value
				readstartdict[exon][bedline[3]] = readstart
				readenddict[exon][bedline[3]] = readend
	return(readstartdict, readenddict)

def calculateMetric(beddict, startdict, enddict, readnum):
	"""Calculate per-exon HITindex values from upstream/downstream junction reads."""
	count = 0
	log_message("[INFO]", "Calculating HITindex metrics for exons.", color="info")
	for gene in beddict:
		for exon in beddict[gene].keys():
			count += 1
			if count % 1000 == 0:
				log_message("[INFO]", f"HITindex metrics processed {count} exons.", color="info")
			exonname = exon +';'+ gene
			exstart = beddict[gene][exon]['start']
			exend = beddict[gene][exon]['end']
			startreads = {k for k,v in startdict[exonname].items() if exstart <= int(v) < exend}
			endreads = {k for k,v in enddict[exonname].items() if exstart <= int(v) < exend}
			bothstartend = startreads & endreads
			startreads_only = startreads - bothstartend
			endreads_only = endreads - bothstartend
			nright = len(startreads_only)
			nleft = len(endreads_only)
			HITindex = -2.0
			if nright + nleft > readnum:
				HITindex =  round(float(nleft - nright)/float(nright + nleft), 3)
			beddict[gene][exon]['name'] = exonname
			beddict[gene][exon]['nleft'] = nleft
			beddict[gene][exon]['nright'] = nright
			beddict[gene][exon]['HIT'] = HITindex
	log_message("[INFO]", f"HITindex metrics processed {count} exons total.", color="info")
	return(beddict)

def piecewise_linear(x, x0, y0, k1, k2):
    """Piecewise linear model used to estimate edge-flagging distance cutoff."""
    return np.piecewise(x, [x < x0, x >= x0], [lambda x:k1*x + y0-k1*x0, lambda x:k2*x + y0-k2*x0])

def edge_flagging(HITdict):
	"""Estimate edge-distance cutoff and annotate HITindex exon dictionaries."""
	flaggingdata = []
	for gene in HITdict:
		exonReads = [k for k in HITdict[gene].keys() if (HITdict[gene][k]['HIT'] > -2)]
		if len(exonReads) == 0: continue
		exonnames = [HITdict[gene][k]['name'] for k in exonReads] 
		exonleft = [HITdict[gene][k]['nleft'] for k in exonReads]
		exonright = [HITdict[gene][k]['nright'] for k in exonReads]
		exonratio = [exonright[i]/(exonright[i]+exonleft[i]) for i in range(len(exonReads))]
		strand = HITdict[gene][exonReads[0]]['strand']
		exstarts = [int(HITdict[gene][k]['start']) for k in exonReads]
		exends = [int(HITdict[gene][k]['end']) for k in exonReads]
		if strand == "+":
			startsite = min(exstarts)
			start_distances = [i - startsite for i in exstarts]
			endsite = max(exends)
			end_distances = [i - endsite for i in exends]
		if strand == "-":
			startsite = max(exends)
			start_distances = [startsite - i for i in exends]
			endsite = min(exstarts)
			end_distances = [endsite - i for i in exstarts]
		if strand not in {"+", "-"}:
			continue
		data = {'name': exonnames, 'ratio': exonratio, 'startdistance': start_distances, 'enddistance': end_distances}
		df = pd.DataFrame(data=data)
		flaggingdata.append(df)
		for i in range(len(exonnames)):
			HITdict[gene][exonReads[i]]['startdistance'] = start_distances[i]
			HITdict[gene][exonReads[i]]['enddistance'] = end_distances[i]

	if len(flaggingdata) == 0:
		log_message("[WARNING]", "No data was added to HITindex edge flagging data.", color="warning")
		return(HITdict, 0)

	flaggingdata = pd.concat(flaggingdata, ignore_index=True)
	ratio = np.array(flaggingdata['ratio'].values, dtype=float)
	distance = np.array(flaggingdata['startdistance'].values, dtype=float)
	try:
		optimresult = optimize.curve_fit(piecewise_linear, distance, ratio)
		param = optimresult[0][0]
	except (RuntimeError, TypeError, ValueError, optimize.OptimizeWarning):
		param = 0
	return(HITdict, param)

def classify_observations(N, D, k, omega):
    """Compute posterior class probabilities for observed HITindex read counts."""
    omega = omega[:-1]
    mode = np.array([1, 0.5, 0, 0])
    k = np.array(list(k) + [-1.0])
    k+=3
    omega = np.array(list(omega) + [1 - np.sum(omega)])

    alpha = mode * (k-2)+1
    beta = (1-mode)*(k-2)+1
    probs = []
    for i in range(len(k)):
        probs.append(omega[i] * scipy.stats.betabinom.pmf(np.array(D), np.array(N), alpha[i], beta[i]))
    
    probs = np.array(probs)
    probs = probs/probs.sum(0)
    return(probs)

def q_posterior(D, N, k_I, mesh=1000):
    """Calculate posterior density over hybrid downstream-read fraction q."""
    q = np.linspace(0, 1, mesh)[1:-1]
    alpha = q * k_I
    beta = (1-q) * k_I
    lk = scipy.stats.betabinom.logpmf(D, N, alpha, beta)
    z = scipy.special.logsumexp(lk)
    return np.exp(lk-z)

def estimate_q(probs, D, N, k_I, mesh = 1000):
    """Summarize hybrid posterior probabilities and downstream-read means."""
    summary_dict = {'post_mean': [],
                    'prob_FI': [],
                    'prob_IL': [],
                    'CI_hyb_lo': [],
                    'CI_hyb_hi': []}
    NH_means = np.array([1, 0.5, 0])
    center = int(mesh/2)-1
    for i in range(len(D)):
        post = q_posterior(D[i], N[i], k_I, mesh=1000)
        prob_FI = post[center:].sum() * probs[-1, i]
        prob_IL = post[:center].sum() * probs[-1, i]
        hybrid_mean = np.dot(post, np.linspace(0, 1, 1000)[1:-1])
        post_mean = (probs[:3,i] * NH_means).sum() + hybrid_mean*probs[-1,i]
        summary_dict['post_mean'].append(post_mean)
        summary_dict['prob_FI'].append(prob_FI)
        summary_dict['prob_IL'].append(prob_IL)
    for key in summary_dict.keys():
        summary_dict[key] = np.array(summary_dict[key])
    return summary_dict

def running_genmodel(exon_outname, fit_iter=100000, seed=None):
	"""Fit the HITindex generative model and append posterior columns to exon output."""
	df_exon = pd.read_csv(exon_outname, sep='\t')
	replicate = os.path.basename(exon_outname)
	log_message(
		"[INFO]",
		f"HITindex model start: {replicate}, exons={len(df_exon)}, ADVI iterations={fit_iter}",
		color="info",
	)
	if df_exon.empty:
		for col in ['PofF', 'PofI', 'PofL', 'PofH', 'PofFI', 'PofIL', 'downstream_fraction', 'HIT_postmean']:
			df_exon[col] = []
		df_exon.to_csv(exon_outname, sep='\t', index=False)
		log_message("[INFO]", f"HITindex model skipped: {replicate}, no exons", color="info")
		return(df_exon)
	if seed is not None:
		random.seed(int(seed))
		np.random.seed(int(seed))
	try:
		import theano
		theano.config.compute_test_value = 'off'
	except ImportError:
		pass
	import pymc3 as pm
	D = np.array(df_exon['nDOWN'].values)
	N = np.array((df_exon['nUP'] + df_exon['nDOWN']).values)

	# ADVI with minibatch fitting.
	min_batch_size = min(1000, len(D))
	Dmini = pm.Minibatch(D, batch_size = min_batch_size)
	Nmini = pm.Minibatch(N, batch_size = min_batch_size)

	with pm.Model() as model_mini:
		# modes for each clas of metaexons held as constants
		mu = np.array([0, 0.5, 1])
		# prior reflecting how concentrated the beta distributions are
		# k at least 2
		k = pm.HalfNormal('k', 1000, shape=3)+3
		# convert mode and concentration parameters to alpha and beta parameters
		alpha = pm.Deterministic('alpha', mu*(k-2)+1)
		beta = pm.Deterministic('beta', (1-mu)*(k-2)+1)
		# omega is frequency of each class, modeled by uniform, Dirichlet prior
		omega = pm.Dirichlet('omega', np.array([1, 1, 1, 1]))
		# explictly define each distribution in the mixture
		betaF = pm.BetaBinomial.dist(alpha[2], beta[2], Nmini)
		betaI = pm.BetaBinomial.dist(alpha[1], beta[1], Nmini)
		betaL = pm.BetaBinomial.dist(alpha[0], beta[0], Nmini)
		betaH = pm.BetaBinomial.dist(1, 1, Nmini)
		# likelihood: observed number of downstream reads arises from mixture of beta-binom dists
		obs = pm.Mixture('D', w = omega, comp_dists = [betaF, betaI, betaL, betaH], observed = Dmini)

	with model_mini:
		fit_kwargs = {"progressbar": False}
		if seed is not None:
			fit_kwargs["random_seed"] = int(seed)
		mean_field = pm.fit(fit_iter, **fit_kwargs)

	mf_trace = mean_field.sample()
	# write back out lines with probabilities
	probs_sum = None
	trace_n = mf_trace['k'].shape[0]
	for i in range(trace_n):
		p_i = classify_observations(N, D, mf_trace['k'][i,:], mf_trace['omega'][i,:])
		if probs_sum is None:
			probs_sum = p_i
		else:
			probs_sum += p_i
	probs = probs_sum / float(trace_n)
	df_exon['PofF'] = np.round(probs[0], 5)
	df_exon['PofI'] = np.round(probs[1], 5)
	df_exon['PofL'] = np.round(probs[2], 5)
	df_exon['PofH'] = np.round(probs[3], 5)

	mean_kI = mf_trace['k'][0].mean()
	summary = estimate_q(probs, D, N, mean_kI)
	df_exon['PofFI'] = np.round(summary['prob_FI'], 5)
	df_exon['PofIL'] = np.round(summary['prob_IL'], 5)
	df_exon['downstream_fraction'] = np.round(summary['post_mean'], 5)
	df_exon['HIT_postmean'] = np.round(-1*(np.array(summary['post_mean'])*2-1), 5)

	df_exon.to_csv(exon_outname, sep='\t', index=False)
	log_message("[INFO]", f"HITindex model done: {replicate}", color="info")

	return(df_exon)

def writeMetrics(HITdict, param, outname):
	"""Write per-exon HITindex metrics to the .exon table."""
	with open(outname, 'w', encoding="utf-8") as outexon:
		outexon.write('exon\tgene\tstrand\tnTXPT\tnFE\tnINTERNAL\tnLE\tnSINGLE\tnUP\tnDOWN\tHITindex\tdist_to_TSS\tdist_to_TES\tedge\n')
		for gene in HITdict:
			for exon in HITdict[gene].keys():
				if HITdict[gene][exon]['HIT'] == -2:
					continue
				infolist = [x.split(':')[1] for x in HITdict[gene][exon]['info'].split(';')]
				edge = 'no'
				if HITdict[gene][exon]['startdistance'] <= param and HITdict[gene][exon]['HIT'] not in {-1, 1}:
					edge = 'yes'
				row = [
					exon,
					gene,
					HITdict[gene][exon]['strand'],
					*infolist,
					str(HITdict[gene][exon]['nleft']),
					str(HITdict[gene][exon]['nright']),
					str(HITdict[gene][exon]['HIT']),
					str(HITdict[gene][exon]['startdistance']),
					str(HITdict[gene][exon]['enddistance']),
					edge,
				]
				outexon.write('\t'.join(row) + '\n')

def readMetrics(metrics):
	"""Read a HITindex .exon metrics table."""
	df = pd.read_csv(metrics, sep='\t')
	return(df)

# Junction-centric AFE/ALE PSI.
def calculatePSI(HITidentify, outname):
	"""Calculate AFE/ALE PSI tables from HITindex-classified exon metrics."""
	afepsidata = []
	alepsidata = []
	afe_cols = ['gene', 'exon', 'strand', 'nTXPT','nFE', 'nUP', 'nDOWN', 'HITindex', 'sumR-L', 'AFEPSI']
	ale_cols = ['gene', 'exon', 'strand', 'nTXPT','nLE', 'nUP', 'nDOWN', 'HITindex', 'sumL-R', 'ALEPSI']

	genelist = HITidentify.gene.unique()
	for genehere in genelist:
		dfhere = HITidentify.loc[HITidentify.gene == genehere]

		afehere = dfhere.loc[(dfhere.ID == 'first') | (dfhere.ID == 'FirstInternal_medium') | (dfhere.ID == 'FirstInternal_high')]
		alehere = dfhere.loc[(dfhere.ID == 'last') | (dfhere.ID == 'InternalLast_medium') | (dfhere.ID == 'InternalLast_high')]

		# AFE PSI
		RLlist = afehere['nDOWN'] - afehere['nUP']
		total_junc = RLlist.sum()
		if len(afehere) > 0 and total_junc > 0:
			if len(afehere) == 1:
				AFEPSI = [1.0]
			else:
				AFEPSI = [float(x)/float(total_junc) for x in RLlist]
			afehere = afehere.copy()
			afehere['sumR-L'] = total_junc
			afehere['AFEPSI'] = AFEPSI
			afehereparse = afehere[afe_cols]
			afepsidata.append(afehereparse)

		# ALE PSI
		LRlist = alehere['nUP'] - alehere['nDOWN']
		total_junc = LRlist.sum()
		if len(alehere) > 0 and total_junc > 0:
			if len(alehere) == 1:
				ALEPSI = [1.0]
			else:
				ALEPSI = [float(x)/float(total_junc) for x in LRlist]
			alehere = alehere.copy()
			alehere['sumL-R'] = total_junc
			alehere['ALEPSI'] = ALEPSI
			alehereparse = alehere[ale_cols]
			alepsidata.append(alehereparse)

	if len(afepsidata) > 0:
		pd.concat(afepsidata, ignore_index=True).to_csv(outname +'.AFEPSI', sep='\t', index=False)
	else:
		pd.DataFrame(columns=afe_cols).to_csv(outname +'.AFEPSI', sep='\t', index=False)
	if len(alepsidata) > 0:
		pd.concat(alepsidata, ignore_index=True).to_csv(outname +'.ALEPSI', sep='\t', index=False)
	else:
		pd.DataFrame(columns=ale_cols).to_csv(outname +'.ALEPSI', sep='\t', index=False)

	return

def probabilistic_round(x):
    """Stochastically round a numeric value using its fractional part."""
    return int(math.floor(x + random.random()))

def sampleReads(allreads, nBOTH):
    """Bootstrap sample synthetic upstream/downstream reads and return HITindex."""
    nchoice = np.array(np.random.choice(allreads, nBOTH, replace=True))
    nchoice_up = np.sum(nchoice < 0)
    nchoice_down = np.sum(nchoice > 0)
    HIT = round(float(nchoice_up - nchoice_down) / float(nBOTH), 3)
    return(HIT)

def bootstrap(exon, n_boot, cutoff):
    """Bootstrap confidence intervals and p-values for one HITindex exon row."""
    nBOTH = int(exon.nUP) + int(exon.nDOWN)
    # Represent upstream/downstream pseudo-reads as signed integers for bootstrap sampling.
    nUP_internal = probabilistic_round(nBOTH/2)
    up_internal = list(range((nUP_internal+1)*-1,0))
    nDOWN_internal = probabilistic_round(nBOTH/2)
    down_internal = list(range(1,(nDOWN_internal+2)))
    up_CI = list(range((int(exon.nUP)+1)*-1, 0))
    down_CI = list(range(1,(int(exon.nDOWN) + 2)))
    allreads_internal = up_internal + down_internal
    allreads_CI = up_CI + down_CI
    bootHITs_internal = []
    bootHITs_CI = []
    for boot in range(0, n_boot): 
        bootHITs_internal.append(sampleReads(allreads_internal, nBOTH))
        bootHITs_CI.append(sampleReads(allreads_CI, nBOTH))
    bootHITs_internal = np.array(bootHITs_internal)
    bootHITs_CI = np.array(bootHITs_CI)
    if float(exon.HITindex) == 0.0:
        pval_internal = np.sum(np.abs(bootHITs_internal) > 0.0) / float(n_boot)
    if float(exon.HITindex) < 0.0:
        pval_internal = np.sum(bootHITs_internal <= exon.HITindex) / float(n_boot)
    if float(exon.HITindex) > 0.0:
        pval_internal = np.sum(bootHITs_internal >= exon.HITindex) / float(n_boot)
    CI75 = str(round(np.percentile(bootHITs_CI, 12.5), 3)) +','+ str(round(np.percentile(bootHITs_CI, 87.5), 3))
    CI90 = str(round(np.percentile(bootHITs_CI, 5), 3)) +','+ str(round(np.percentile(bootHITs_CI, 95), 3))
    CI95 = str(round(np.percentile(bootHITs_CI, 2.5), 3)) +','+ str(round(np.percentile(bootHITs_CI, 97.5), 3))
    pval_CI = np.sum(np.abs(bootHITs_CI) < abs(cutoff)) / float(n_boot)
    return CI75, CI90, CI95, pval_CI, pval_internal

def significance(HITcombo, n_boot, cutoff, outname, seed=None):
    """Run HITindex bootstrap significance and write confidence columns."""
    replicate = os.path.basename(outname)
    log_message(
        "[INFO]",
        f"HITindex bootstrap start: {replicate}, exons={len(HITcombo)}, bootstrap_n={n_boot}",
        color="info",
    )
    if seed is not None:
        random.seed(int(seed))
        np.random.seed(int(seed))
    if HITcombo.empty:
        for col in ['CI_75', 'CI_90', 'CI_95', 'pval_CI', 'pval_internal']:
            HITcombo[col] = []
        HITcombo.to_csv(outname, sep='\t', index=False)
        log_message("[INFO]", f"HITindex bootstrap skipped: {replicate}, no exons", color="info")
        return HITcombo
    HITcombo[['CI_75', 'CI_90', 'CI_95', 'pval_CI', 'pval_internal']] = HITcombo[['nUP', 'nDOWN', 'HITindex']].apply(bootstrap, axis=1, result_type="expand", n_boot=n_boot, cutoff=cutoff)
    HITcombo.to_csv(outname, sep='\t', index=False)	
    log_message("[INFO]", f"HITindex bootstrap done: {replicate}", color="info")
    return HITcombo

def _resolve_hit_ci_column(hit_ci):
    value = str(hit_ci or "none").strip()
    if value.lower() == "none":
        return None
    normalized = value.lower().replace("ci_", "").replace("ci", "").replace("%", "")
    aliases = {
        "0.75": "CI_75",
        ".75": "CI_75",
        "75": "CI_75",
        "0.9": "CI_90",
        "0.90": "CI_90",
        ".9": "CI_90",
        ".90": "CI_90",
        "90": "CI_90",
        "0.95": "CI_95",
        ".95": "CI_95",
        "95": "CI_95",
    }
    if normalized not in aliases:
        raise RuntimeError(f"Unsupported HIT_CI value: {hit_ci}. Expected none, 0.75, 0.90, or 0.95.")
    return aliases[normalized]

def _ci_spans_zero(ci_value):
    try:
        low_s, high_s = str(ci_value).split(",", 1)
        return float(low_s) < 0 < float(high_s)
    except (TypeError, ValueError):
        return False

def call_terminal(HITcombo, paramdict, outname):
    """Assign terminal/internal HITindex labels and write the classified table."""
    # set all to internal
    HITcombo['ID'] = "internal"
    # FI_med
    HITcombo.loc[(HITcombo.HITindex <= float(paramdict['HIThybrid'])*-1) & (HITcombo.PofFI >= float(paramdict['prob_med'])), 'ID'] = 'FirstInternal_medium'
    # FI_high
    HITcombo.loc[(HITcombo.ID == 'FirstInternal_medium') & (HITcombo.PofFI >= float(paramdict['prob_high'])), 'ID'] = 'FirstInternal_high'
    # IL_med
    HITcombo.loc[(HITcombo.HITindex >= float(paramdict['HIThybrid'])) & (HITcombo.PofIL >= float(paramdict['prob_med'])), 'ID'] = 'InternalLast_medium'
    # IL_high
    HITcombo.loc[(HITcombo.ID == 'InternalLast_medium') & (HITcombo.PofIL >= float(paramdict['prob_high'])), 'ID'] = 'InternalLast_high'
    # first
    HITcombo.loc[(HITcombo.HITindex <= float(paramdict['HITterminal'])*-1) & (HITcombo.pval_internal <= float(paramdict['HITpval'])), 'ID'] = 'first'
    # last
    HITcombo.loc[(HITcombo.HITindex >= float(paramdict['HITterminal'])) & (HITcombo.pval_internal <= float(paramdict['HITpval'])), 'ID'] = 'last'
    # confidence intervals
    ci_column = _resolve_hit_ci_column(paramdict.get('HIT_CI', 'none'))
    if ci_column:
        if ci_column not in HITcombo.columns:
            raise RuntimeError(f"HIT_CI requested {ci_column}, but column is missing from HITindex table.")
        HITcombo.loc[HITcombo[ci_column].apply(_ci_spans_zero), 'ID'] = 'internal'
    # add new column with assignments based on position in gene
    HITcombo['ID_position'] = HITcombo['ID'] 
    # first if upstream most exon
    HITcombo.loc[(HITcombo.dist_to_TSS == 0), 'ID_position'] = 'first'
    # last if downstream most exon
    HITcombo.loc[(HITcombo.dist_to_TES == 0), 'ID_position'] = 'last'
    # write file to csv
    HITcombo.to_csv(outname, sep='\t', index=False)	
    return
