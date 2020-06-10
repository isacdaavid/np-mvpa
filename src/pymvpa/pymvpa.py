#!/usr/bin/python2

# author: Isaac David <isacdaavid@at@isacdaavid@dot@info>
# license: GPLv3 or later

import sys
sys.path.append('/home/inb/lconcha/fmrilab_software/miniconda2/lib/python2.7/site-packages')
import matplotlib
matplotlib.use('Agg') # force matplotlib not to use any Xwindows backend
import matplotlib.pyplot as plt
from mvpa2.suite import *
from itertools import chain, combinations
import re

# arguments passed to script
ATTR_FNAME = sys.argv[1]
BOLD_FNAME = sys.argv[2]
MASK_FNAME = sys.argv[3]
OUTDIR = sys.argv[4]
REDUCED_BOLD_FNAME = sys.argv[5]
CLASSES = eval("['" + re.sub(",", "','", sys.argv[6]) + "']")

STEP = 2000 # time step between different Hemodynamic Response delays (ms)
TIME_START = 0 # first HR delay to test for (ms)
TIME_LIMIT = 10001 # maximum HR delay to test for (ms)
SLICE_TIMING_REFERENCE = +1000 # ms
DELAYS = range(TIME_START, TIME_LIMIT, STEP)

PERMUTATIONS = 5000 # label permutations used to estimate null accuracy distrib
ANOVA_SELECTION = 1 # proportion of features to work with

# WARNING: assigning an existing pyMVPA Dataset object (or one of its attributes)
#          to a new variable/attribute is a call by reference. for actual copies
#          you must use dataset.copy(deep=True, ...)

################################################################################
# volume labeling
################################################################################

def label(ds, attr, slice_timing_reference, hrf_delay):
        # correct for preambulus time-shift, start onset_times count at 0
        onset_times = np.array(attr.onset_time) - attr.onset_time[0]
        attr.pop('onset_time')
        attr['onset_time'] = np.ndarray.tolist(onset_times)

        # convert volume times to ms, as eprime, then time-slice-shift
        # (tr=2s, interleaved)
        ds.sa.time_coords = (ds.sa.time_coords * 1000) + slice_timing_reference

        indices = list()
        # exclude volumes prior to 1 hrf_delay, since they are not supposed
        # to reflect a maximum of task-related HRF activity
        for vol_time in ds[ds.sa.time_coords >= hrf_delay].sa.time_coords:
                # select latest event to occur before (vol_time - hrf_delay)
                match_floor = \
                        max(onset_times[onset_times <= (vol_time - hrf_delay)])
                indices.append(attr.onset_time.index(match_floor))

        # although we won't return them, excluded premature volumes still are
        # assigned the attributes of the first non-excluded one
        useless_vols = 0
        while(len(indices) < len(ds)):
                indices.insert(0, indices[0])
                useless_vols = useless_vols + 1

        # label volumes according to attribute sublist
        for at in attr.keys():
                ds.sa[at] = [attr[at][index] for index in indices]

        return ds[useless_vols:]

################################################################################
# classification
################################################################################

def subsample(ds0):
        ds1 = ds0[{'emotion': CLASSES}]
        ds1.sa['targets'] = ds1.sa.emotion # this is the label/target attribute
        return ds1

def train():
        clf = LinearCSVMC()
        fsel = SensitivityBasedFeatureSelection(
                        OneWayAnova(),
                        FractionTailSelector(ANOVA_SELECTION,
                                             mode = 'select',
                                             tail = 'upper'))
        fclf = FeatureSelectionClassifier(clf, fsel)
        cv = CrossValidation(fclf, NFoldPartitioner(attr = 'block'),
                             errorfx = lambda p, t: np.mean(p == t),
                             enable_ca=['stats'])
        return fclf,cv

################################################################################
# Monte-Carlo null hypothesis estimation
################################################################################

def null_cv(permutations = PERMUTATIONS):
        repeater = Repeater(count = permutations)
        permutator = AttributePermutator('targets',
                                         limit={'partitions': 1}, count = 1)
        clf = LinearCSVMC()
        fsel = SensitivityBasedFeatureSelection(
                        OneWayAnova(),
                        FractionTailSelector(ANOVA_SELECTION,
                                             mode = 'select',
                                             tail = 'upper'))
        fclf = FeatureSelectionClassifier(clf, fsel)
        partitioner = NFoldPartitioner(attr = 'block')
        cv = CrossValidation(fclf,
                             ChainNode([partitioner, permutator],
                                       space = partitioner.get_space()),
                             errorfx = lambda p, t: np.mean(p == t),
                             postproc = mean_sample())
        distr_est = MCNullDist(repeater, tail = 'right',
                               measure = cv,
                               enable_ca = ['dist_samples'])
        cv_mc = CrossValidation(fclf,
                                partitioner,
                                errorfx = lambda p, t: np.mean(p == t),
                                postproc = mean_sample(),
                                null_dist = distr_est,
                                enable_ca = ['stats'])
        return fclf,cv_mc

def make_null_dist_plot(dist_samples, empirical, nclasses):
     pl.hist(dist_samples, bins = 100, normed = True, alpha = 0.8)
     pl.axvline(empirical, color='red')
     # a priori chance-level
     pl.axvline(1.0 / nclasses, color='black', ls='--')
     # scale x-axis to full range of possible error values
     pl.xlim(0,1)
     pl.xlabel('Average cross-validated classification error')

################################################################################
# sensitivity analysis
################################################################################

# - remove sign (there's no interpretation to feature-importance direction in
#   the orthogonal vector to SVM hyperplane, other than encoding class)
# - L2-normalize to make sure vector sum is meaningful
# - sum
# - rescale to maximum weight (summed masks will be comparable operators)
# - optionally, return n most significant weights
def normalize_weights(weight_lists, significance = 1):
        for i in range(0, len(weight_lists)):
                weight_lists[i] = l2_normed(abs(weight_lists[i]))
        if len(weight_lists) > 1:
                total = np.sum(weight_lists, axis = 0)
        else:
                total = weight_lists[0]
        total /= max(total)
        ntile = np.sort(total)[-int(round(len(total) * significance))]
        return np.array([(0 if (x < ntile) else x) for x in total])

def sensitivity_maps_aux(model, ds):
        analyzer = model.get_sensitivity_analyzer()
        return analyzer(ds)

def powerset(iterable):
        s = list(iterable)
        return chain.from_iterable(combinations(s, r) for r in range(2, len(s) + 1))

def sanitize_mask_name(string):
        return re.sub("[' \[\]]", '', string)

# outputs the computed "activation" maps (rather, sensitivity masks)
def sensitivity_maps(model, ds):
        sens = sensitivity_maps_aux(model, ds)
        masks = dict()
        for comb in powerset(CLASSES):
                subset_pairs = []
                for i in range(0, len(sens.targets)):
                        if set(sens.targets[i]).issubset(comb):
                                subset_pairs.append(sens[i].samples[0])
                name = sanitize_mask_name(str(sorted(comb)))
                masks[name] = normalize_weights(np.array(subset_pairs))
        return masks

# percentage of voxels with non-zero weights
def non_empty_weights_proportion(weights):
        return len(weights[weights != 0.0]) / float(len(weights))

################################################################################
# main
################################################################################

# load psychopy events/design matrix (aka target attributes)
attr = SampleAttributes(ATTR_FNAME,
                        header = ['onset_time', 'emotion', 'block'])

result_dist = []
fo = open(OUTDIR + "/result-time-series.txt", "w+")

for delay in DELAYS:
        orig = fmri_dataset(BOLD_FNAME, mask = MASK_FNAME)
        ds = orig
        if REDUCED_BOLD_FNAME != BOLD_FNAME :
                # load low-dimensional version of dataset
                ds = Dataset(np.genfromtxt(REDUCED_BOLD_FNAME, delimiter=' '))
                ds.sa = orig.sa
        ds2 = label(ds, attr, SLICE_TIMING_REFERENCE, delay)
        ds3 = subsample(ds2)
        model,validator = train()
        results = validator(ds3)
        result_dist.append(np.mean(results))
        masks = sensitivity_maps(model, ds3)
        name = sanitize_mask_name(str(sorted(CLASSES)))
        print(np.mean(results))
        fo.writelines(str(ds3.nsamples / len(CLASSES)) + " " + str(np.mean(results)) + " " +
                      str(non_empty_weights_proportion(masks[name])) + "\n")

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
ds3 = subsample(ds2)
# null accuracy estimation using Monte-Carlo method
model,validator = null_cv(PERMUTATIONS)
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

masks = sensitivity_maps(model, ds3)
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
