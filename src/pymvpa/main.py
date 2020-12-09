#!/usr/bin/python2

# author: Isaac David <isacdaavid@at@isacdaavid@dot@info>
# license: GPLv3 or later

import socket
print(socket.gethostname())

import sys
sys.path.append('/home/inb/lconcha/fmrilab_software/miniconda2/lib/python2.7/site-packages')
import matplotlib
matplotlib.use('Agg') # force matplotlib not to use any Xwindows backend
import matplotlib.pyplot as plt
import re
from mvpa2.suite import *

# local imports
from libmvpa import *

# arguments passed to script
ATTR_FNAME = sys.argv[1]
BOLD_FNAME = sys.argv[2]
MASK_FNAME = sys.argv[3]
OUTDIR = sys.argv[4]
REDUCED_BOLD_FNAME = sys.argv[5]
CLASSES = eval("['" + re.sub(",", "','", sys.argv[6]) + "']")
PERMUTATIONS = int(sys.argv[7])

# other global constants
# TODO: don't hard-code. infer from file
STEP = 2000 # time step between different Hemodynamic Response delays (ms)
TIME_START = 0 # first HR delay to test for (ms)
TIME_LIMIT = 10001 # maximum HR delay to test for (ms)
SLICE_TIMING_REFERENCE = +1000 # ms
DELAYS = range(TIME_START, TIME_LIMIT, STEP)
PREDEFINED_DELAY = 4000 # ms

ANOVA_SELECTION = 1 # proportion of features to work with

# WARNING: assigning an existing pyMVPA Dataset object (or one of its attributes)
#          to a new variable/attribute is a call by reference. for actual copies
#          you must use dataset.copy(deep=True, ...)

################################################################################
# main
################################################################################

# load psychopy events/design matrix (aka target attributes)
attr = SampleAttributes(ATTR_FNAME,
                        header = ['onset_time', 'label', 'block'])

result_dist = []
fo = open(OUTDIR + "/result-time-series.csv", "w+")
fo.writelines("sample_size mean_accuracy voxel_prop ms fold_accuracies\n")

for delay in DELAYS:
        orig = fmri_dataset(BOLD_FNAME, mask = MASK_FNAME)
        ds = orig
        if REDUCED_BOLD_FNAME != BOLD_FNAME :
                # load low-dimensional version of dataset
                ds = Dataset(np.genfromtxt(REDUCED_BOLD_FNAME, delimiter=' '))
                ds.sa = orig.sa
        ds2 = label(ds, attr, SLICE_TIMING_REFERENCE, delay)
        ds3 = subsample(ds2, CLASSES)
        model,validator = train(ANOVA_SELECTION)
        results = validator(ds3)
        result_dist.append(np.mean(results))
        masks = sensitivity_maps(sensitivity_maps_aux(model, ds), CLASSES)
        name = sanitize_mask_name(str(sorted(CLASSES)))
        accuracies = ','.join([str(x) for x in np.ravel(results.samples)])
        line = str(ds3.nsamples / len(CLASSES)) + " " + str(np.mean(results)) \
               + " " + str(non_empty_weights_proportion(masks[name])) + " " \
               + str(delay) + " " + accuracies + "\n"
        fo.writelines(line)
        print(str(delay)+"/"+str(DELAYS[-1])+" ms : acc "+str(np.mean(results)))

fo.close()

plt.plot(result_dist)
plt.xlabel('labeling delay')
plt.ylabel('mean cross-validated classification accuracy')
plt.savefig(OUTDIR + '/result-time-series.svg')
plt.close()

plt.hist(result_dist, bins = 1000)
plt.xlabel('mean cross-validated classification accuracy')
plt.ylabel('frequency')
plt.savefig(OUTDIR + '/result-dist.svg')
plt.close()

# best model ###################################################################

optimal_delay = (result_dist.index(max(result_dist)) * STEP) + TIME_START
orig = fmri_dataset(BOLD_FNAME, mask = MASK_FNAME)
ds = orig
if REDUCED_BOLD_FNAME != BOLD_FNAME :
        # load low-dimensional version of dataset
        ds = Dataset(np.genfromtxt(REDUCED_BOLD_FNAME, delimiter=' '))
        ds.sa = orig.sa
ds2 = label(ds, attr, SLICE_TIMING_REFERENCE, PREDEFINED_DELAY)
ds3 = subsample(ds2, CLASSES)

# representation similarity analysis
ds4 = vstack([ds3[{'label': [l]}] for l in CLASSES])
# FIXME: Einsten sum inside RSA_matrix() consumes too much memory apparently
# plot_RSA_matrix(pattern = RSA_matrix(ds4.samples, distance = "euclidean"),
#                 labels = ds4.sa.label,
#                 path = OUTDIR + '/RSA_euclidean')
plot_RSA_matrix(RSA_matrix(ds4.samples - np.mean(ds4.samples, axis = 0),
                           distance = "pearson"),
                ds4.sa.label,
                OUTDIR + '/RSA_pearson')

# tune regularizing hyperparameter 'C'
hyper = []
# cees = [1, .1, .01, .001, .0001, -1]
cees = [-1]
for c in cees:
        oldmodel,oldvalidator = train(ANOVA_SELECTION, C = c)
        oldresults = oldvalidator(ds3)
        hyper.append(np.mean(oldresults))
        print("C={0} : acc {1}".format(c, np.mean(oldresults)))
optimal_hyper = cees[hyper.index(max(hyper))]
oldmodel,oldvalidator = train(ANOVA_SELECTION, C = optimal_hyper)
oldresults = oldvalidator(ds3)

np.savetxt(OUTDIR + "/results-folds-best.csv",
           oldresults.samples,
           fmt='%10.5f')

masks = sensitivity_maps(sensitivity_maps_aux(oldmodel, ds3), CLASSES)

# null accuracy estimation using Monte-Carlo method
model,validator,pvals,svm_max = null_cv(PERMUTATIONS,
                                        CLASSES,
                                        masks,
                                        ANOVA_SELECTION,
                                        C = optimal_hyper)
results = validator(ds3)
pvals = {k : pvals[k]/(PERMUTATIONS * max(ds3.sa.block)) for k in pvals.keys()}

fo = open(OUTDIR + "/conf-matrix-stats.txt", "w+")
fo.writelines(validator.ca.stats.as_string(description = True))
fo.close()

np.savetxt(OUTDIR + "/conf-matrix.csv",
           validator.ca.stats._ConfusionMatrix__matrix,
           header = ' '.join(validator.ca.stats._ConfusionMatrix__labels),
           fmt='%i')

validator.ca.stats.plot()
plt.savefig(OUTDIR + '/conf-matrix.svg')
plt.close()

fo = open(OUTDIR + '/null-dist.csv', "w+")
fo.writelines("\n".join(str(i) \
        for i in validator.null_dist.ca.dist_samples.samples.tolist()[0][0]))
fo.close()

plt.hist(np.ravel(validator.null_dist.ca.dist_samples),
         bins = 100, normed = True, alpha = 0.8)
plt.axvline(np.mean(results), color='red')
# a priori chance-level
plt.axvline(1.0 / len(CLASSES), color='black', ls='--')
# scale x-axis to full range of possible error values
plt.xlim(0,1)
plt.xlabel('mean cross-validated classification accuracy')
plt.ylabel('frequency')
plt.savefig(OUTDIR + '/null-dist.svg')
plt.close()

# sensitivity maps #############################################################

masks = sensitivity_maps(sensitivity_maps_aux(model, ds), CLASSES)
all_weights = masks[sanitize_mask_name(str(sorted(CLASSES)))]

# distribution of non-zero weights, normalized to the maximum one
plt.hist(all_weights[all_weights != 0] / max(all_weights), bins = 50)
plt.savefig(OUTDIR + '/weights-dist.svg')
plt.close()

# export sensitivity maps and p-value maps
for name in masks:
        if REDUCED_BOLD_FNAME != BOLD_FNAME :
                # load copies of parcellation to pour our values in
                parcells_s = fmri_dataset(MASK_FNAME)
                parcells_p = fmri_dataset(MASK_FNAME)
                parcells_s.samples = parcells_s.samples.astype(float)
                parcells_p.samples = parcells_p.samples.astype(float)
                for i in range(0, len(masks[name])):
                        parcells_s.samples[parcells_s.samples == i + 1] = \
                                masks[name][i]
                        parcells_p.samples[parcells_p.samples == i + 1] = \
                                pvals[name][i]
                nimg_s = map2nifti(parcells_s)
                nimg_p = map2nifti(parcells_p)
        else:
                nimg_s = map2nifti(ds, masks[name])
                nimg_p = map2nifti(ds, pvals[name])
        nimg_s.to_filename(OUTDIR + '/' + name + '-weights.nii.gz')
        nimg_p.to_filename(OUTDIR + '/' + name + '-uncorrected_pvals.nii.gz')
        fo = open(OUTDIR + '/' + name + "-max_weight_null_dist.csv", "w+")
        fo.writelines("\n".join(str(x) for x in svm_max[name]))
        fo.close()
