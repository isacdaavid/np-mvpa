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
from itertools import chain, combinations

def label(ds, attr, slice_timing_reference, hrf_delay):
    """Append labels to a dataset sample attributes at ds.sa['label']

    Args:
      ds (pyMVPA Dataset): dataset which labels will be added to
      attr (pyMVPA SampleAttributes): dict with labels, onset_time and block
      slice_timing_reference (int): time of reference volume in miliseconds
      hrf_delay: withold label by this amount of miliseconds after onset

    Returns: pyMVPA Dataset with .sa['label'] attribute and adjusted
      .sa['onset_times'] and sa['time_coords']
    """

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
        match_floor = max(onset_times[onset_times <= (vol_time - hrf_delay)])
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

def subsample(ds0, classes):
        ds1 = ds0[{'label': classes}]
        ds1.sa['targets'] = ds1.sa.label # this is the label/target attribute
        return ds1

def train(anova_selection):
        clf = LinearCSVMC()
        fsel = SensitivityBasedFeatureSelection(
                        OneWayAnova(),
                        FractionTailSelector(anova_selection,
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

def null_cv(permutations, anova_selection = 1.0):
        repeater = Repeater(count = permutations)
        permutator = AttributePermutator('targets',
                                         limit={'partitions': 1}, count = 1)
        clf = LinearCSVMC()
        fsel = SensitivityBasedFeatureSelection(
                        OneWayAnova(),
                        FractionTailSelector(anova_selection,
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
def sensitivity_maps(model, ds, classes):
        sens = sensitivity_maps_aux(model, ds)
        masks = dict()
        for comb in powerset(classes):
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
# representation similarity analysis
################################################################################

def RSA_matrix(lds, distance = "euclidean"):
    """Computes some distance between all pattern/sample pairs
    and stores them in matrix form.

    Args:
      lds (ndarray): brain samples as row-vectors
      distance (string): "euclidean" or "pearson"

    Returns:
      ndarray: representation-similarity matrix with distances between patterns
    """
    if distance == "pearson":
        return np.corrcoef(lds)
    lds2 = lds.reshape(lds.shape[0], 1, lds.shape[1])
    return np.sqrt(np.einsum('ijk, ijk->ij', lds - lds2, lds - lds2))

def plot_RSA_matrix(pattern, labels, path):
    """Plot representation similarity matrix and also save as CSV

    Args:
      pattern (ndarray): pattern distance matrix as computed by RSA_matrix()
      labels (ndarray): labels vector in the same order as in `pattern`
      path (string): file name where to save path.svg and path.csv matrices
    """
    maxlen = len(np.array2string(labels))
    header = filter(lambda ch: ch not in "[]",
                    np.array2string(labels, max_line_width = maxlen))
    np.savetxt(path + ".csv", pattern, header = header)
    counts = np.unique(labels, return_counts = True)
    classes = counts[0]
    classes_counts = counts[1]
    plt.matshow(pattern)
    plt.xticks(np.cumsum(classes_counts), classes)
    plt.yticks(np.cumsum(classes_counts), classes)
    plt.colorbar()
    plt.savefig(path + ".svg")
    plt.close()
