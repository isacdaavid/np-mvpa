#!/usr/bin/python2

# author: Isaac David <isacdaavid@at@isacdaavid@dot@info>
# license: GPLv3 or later

import sys
sys.path.append('/home/inb/lconcha/fmrilab_software/miniconda2/lib/python2.7/site-packages')
import matplotlib
matplotlib.use('Agg') # force matplotlib not to use any Xwindows backend
import matplotlib.pyplot as plt
import re
from mvpa2.suite import *

# local imports
from rsa import *

ATTR_FNAME = "../../data/psychopy/01.csv"
BOLD_FNAME = "../../data/pymvpa/01/concat-brain-norm.nii.gz"
MASK_FNAME = "../../data/feat/01/volbrain-mask.tmp.nii.gz"
OUTDIR = "../../out/pymvpa/whole/01/faces/scrambled,neutral"
REDUCED_BOLD_FNAME = "../../data/pymvpa/01/concat-brain-norm.nii.gz"
CLASSES = ["scrambled", "neutral"]
PERMUTATIONS = 5

# arguments passed to script
ATTR_FNAME = sys.argv[1]
BOLD_FNAME = sys.argv[2]
MASK_FNAME = sys.argv[3]
OUTDIR = sys.argv[4]
REDUCED_BOLD_FNAME = sys.argv[5]
CLASSES = eval("['" + re.sub(",", "','", sys.argv[6]) + "']")
PERMUTATIONS = int(sys.argv[7])

# other global constants
STEP = 2000 # time step between different Hemodynamic Response delays (ms)
TIME_START = 0 # first HR delay to test for (ms)
TIME_LIMIT = 10001 # maximum HR delay to test for (ms)
SLICE_TIMING_REFERENCE = +1000 # ms
DELAYS = range(TIME_START, TIME_LIMIT, STEP)

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
fo = open(OUTDIR + "/result-time-series.txt", "w+")
fo.writelines("sample_size mean_accuracy voxel_prop ms\n")

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
        masks = sensitivity_maps(model, ds3, CLASSES)
        name = sanitize_mask_name(str(sorted(CLASSES)))
        line = str(ds3.nsamples / len(CLASSES)) + " " + str(np.mean(results)) \
               + " " + str(non_empty_weights_proportion(masks[name])) + " " \
               + str(delay) + "\n"
        fo.writelines(line)
        print(str(delay)+"/"+str(DELAYS[-1])+" ms : acc "+str(np.mean(results)))

fo.close()

plt.plot(result_dist)
plt.savefig(OUTDIR + '/result-time-series.svg')
plt.close()
plt.hist(result_dist, bins = 1000)
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
ds2 = label(ds, attr, SLICE_TIMING_REFERENCE, optimal_delay)
ds3 = subsample(ds2, CLASSES)

# representation similarity analysis
ds4 = vstack([ds3[{'label': [l]}] for l in CLASSES])
plot_RSA_matrix(pattern = RSA_matrix(ds4.samples, distance = "euclidean"),
                labels = ds4.sa.label,
                path = OUTDIR + '/rsa_euclidean')
plot_RSA_matrix(RSA_matrix(ds4.samples - np.mean(ds4.samples, axis = 0),
                           distance = "pearson"),
                ds4.sa.label,
                OUTDIR + '/rsa_pearson')

# null accuracy estimation using Monte-Carlo method
model,validator = null_cv(PERMUTATIONS, ANOVA_SELECTION)
results = validator(ds3)

fo = open(OUTDIR + "/conf-matrix.txt", "w+")
fo.writelines(validator.ca.stats.as_string(description = True))
fo.close()

validator.ca.stats.plot()
plt.savefig(OUTDIR + '/conf-matrix.svg')
plt.close()

fo = open(OUTDIR + '/null-dist.txt', "w+")
fo.writelines("\n".join(str(i) \
        for i in validator.null_dist.ca.dist_samples.samples.tolist()[0][0]))
fo.close()

make_null_dist_plot(np.ravel(validator.null_dist.ca.dist_samples),
                    np.mean(results),
                    len(CLASSES))
plt.savefig(OUTDIR + '/null-dist.svg')
plt.close()

# sensitivity maps #############################################################

masks = sensitivity_maps(model, ds3, CLASSES)
all_weights = masks[sanitize_mask_name(str(sorted(CLASSES)))]

# distribution of non-zero weights, normalized to the maximum one
plt.hist(all_weights[all_weights != 0] / max(all_weights), bins = 50)
plt.savefig(OUTDIR + '/weights-dist.svg')
plt.close()

# export sensitivity maps
for name in masks:
        if REDUCED_BOLD_FNAME != BOLD_FNAME :
                parcellation = fmri_dataset(MASK_FNAME)
                parcellation.samples = parcellation.samples.astype(float)
                for i in range(0, len(masks[name])):
                        parcellation.samples[parcellation.samples == i + 1] = \
                                masks[name][i]
                nimg = map2nifti(parcellation)
        else:
                nimg = map2nifti(ds, masks[name])
        nimg.to_filename(OUTDIR + '/' + name + '.nii.gz')
