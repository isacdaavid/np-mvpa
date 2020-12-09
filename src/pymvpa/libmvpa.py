#!/usr/bin/python2

# author: Isaac David <isacdaavid@at@isacdaavid@dot@info>
# license: GPLv3 or later

import sys
# sys.path.append('/home/inb/lconcha/fmrilab_software/miniconda2/lib/python2.7/site-packages')

import matplotlib
matplotlib.use('Agg') # force matplotlib not to use any Xwindows backend
import matplotlib.pyplot as plt
import re
from mvpa2.suite import *
from itertools import chain, combinations

def label(ds, attr, slice_timing_reference, hrf_delay):
    """Append labels to a dataset sample attributes at ds.sa['label']

    Args:
      ds (pyMVPA Dataset): dataset to which labels will be added
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
        """Subselect sample vectors based on class label

        Args:
          ds0 (pyMVPA Dataset): dataset to be subsampled
          classes (list): list of class/label strings

        Returns: (pyMVPA Dataset)
        """
        ds1 = ds0[{'label': classes}]
        ds1.sa['targets'] = ds1.sa.label # this is the label/target attribute
        return ds1

def train(anova_selection, C = -1.0):
        """Train and cross-validate SVM with optional voxel selection

        Args:
          anova_selection (float): select top F-scoring voxels by this amount (0, 1]
          C (float): SVM regularization hyperparameter (default -1. See LinearCSVMC())

        Returns (tuple):
          fclf (FeatureSelectionClassifier): complex classifier object schema
          cv (CrossValidation): crossvalidation object, expects Dataset argument
                                to actually trigger training
        """
        clf = LinearCSVMC(C = C)
        fsel = SensitivityBasedFeatureSelection(
                        OneWayAnova(),
                        FractionTailSelector(anova_selection,
                                             mode = 'select',
                                             tail = 'upper'))
        fclf = FeatureSelectionClassifier(clf, fsel)
        cv = CrossValidation(fclf,
                             NFoldPartitioner(attr = 'block'),
                             errorfx = lambda p, t: np.mean(p == t),
                             enable_ca=['stats'])
        return fclf,cv

################################################################################
# Monte-Carlo null hypothesis estimation
################################################################################

def null_cv(permutations, classes, weights, anova_selection = 1.0, C = -1.0):
        """Like train() but with permutation tests.
           Train and cross-validate SVM with optional voxel selection, estimate
           null accuracy distribution.

        Args:
          permutations (int): number of permutations for null model estimation
          classes (list): list of class/label strings
          weights (dict): sensitivity_maps() dict to read class combinations and
                          populate subject SVM-max distro and p-val map
          anova_selection (float): select top F-scoring voxels by this amount (0, 1]
          C (float): SVM regularization hyperparameter (default -1. See LinearCSVMC())

        Returns (tuple):
          fclf (FeatureSelectionClassifier): complex classifier object schema
          cv_mc (CrossValidation): crossvalidation object, expects Dataset argument
                                to actually trigger training
          PVALS (dict): uncorrected voxelwise p-val maps of SVM parameters (ndarray),
                        one per class combination key
          SVM_MAX (dict): lists of maximum observed SVM parameter at each
                          permutation, one per class combination key
        """
        global PERMUTATIONS, PERMUTED_FOLD, PVALS, SVM_MAX, CLASSES, WEIGHTS
        PERMUTATIONS = permutations
        CLASSES = classes
        WEIGHTS = weights
        PERMUTED_FOLD = 0
        PVALS = {k : np.zeros(weights[k].shape[0]) for k in weights.keys()}
        SVM_MAX = {k : [] for k in weights.keys()}
        repeater = Repeater(count = permutations)
        permutator = AttributePermutator('targets',
                                         limit={'partitions': 1}, count = 1)
        clf = LinearCSVMC(C = C)
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
                             postproc = mean_sample(),
                             callback = p_vals_map)
        distr_est = MCNullDist(repeater, tail = 'right',
                               measure = cv,
                               enable_ca = ['dist_samples'])
        cv_mc = CrossValidation(fclf,
                                partitioner,
                                errorfx = lambda p, t: np.mean(p == t),
                                postproc = mean_sample(),
                                null_dist = distr_est,
                                enable_ca = ['stats'])
        return fclf,cv_mc,PVALS,SVM_MAX

def p_vals_map(data, node, result):
    """update p-val maps and SVM_MAX lists upon permutation trigger. Side-effect
       -only function, no return value.

    Args:
      data (pyMVPA Dataset): automatically bound, see CrossValidation() argument
                             `postproc`
      node (FeatureSelectionClassifier): automatically bound, see
                                         CrossValidation() argument `postproc`
      result: automatically bound, see CrossValidation() argument `postproc`
    """
    global PERMUTATIONS, PERMUTED_FOLD, PVALS, SVM_MAX, CLASSES, WEIGHTS
    PERMUTED_FOLD += 1
    if PERMUTED_FOLD % max(data.sa.block) == 0:
        permutation = PERMUTED_FOLD / max(data.sa.block)
        print "permutation {0}/{1}\r".format(permutation, PERMUTATIONS),
        if permutation == PERMUTATIONS:
            print("")
    nullmap = sensitivity_maps(sensitivity_maps_aux(node.measure, data),
                               CLASSES)
    for k in WEIGHTS.keys():
        PVALS[k] += abs(nullmap[k]) >= abs(WEIGHTS[k])
        SVM_MAX[k].append(max(abs(nullmap[k])))

################################################################################
# sensitivity analysis
################################################################################

def sensitivity_maps_aux(model, ds):
        """Auxilary function for sensitivity_maps(). Extract sensitivity analyzer

        Args:
          model (FeatureSelectionClassifier): complex classifier object schema
          ds (pyMVPA Dataset): dataset the model was trained with

        Returns: model.get_sensitivity_analyzer() object
        """
        analyzer = model.get_sensitivity_analyzer(force_train = False)
        return analyzer(ds)

def sensitivity_maps(sens, classes):
        """Obtain SVM weights (activation maps/sensitivity maps) from model

        Args:
          sens: sensitivity analyzer object, as obtained from sensitivity_maps_aux()
          classes (list): list of class/label strings

        Returns: (dict)
        """
        masks = dict()
        for comb in powerset(classes):
                subset_pairs = []
                for i in range(0, len(sens.targets)):
                        if set(sens.targets[i]).issubset(comb):
                                subset_pairs.append(sens[i].samples[0])
                name = sanitize_mask_name(str(sorted(comb)))
                masks[name] = normalize_weights(np.array(subset_pairs))
        return masks


def normalize_weights(weight_lists, significance = 1):
        """L2-normalize SVM weight vector(s), add them together in case of
           non-binary classification. Optionally. return n most significant
           weights

        Args:
          weight_list (list): ndarrays with SVM weights
          significance (float): select top voxels by this amount (0, 1]

        Returns: (ndarray)
        """
        for i in range(0, len(weight_lists)):
                weight_lists[i] = l2_normed(weight_lists[i])
        if len(weight_lists) > 1:
                total = np.sum(weight_lists, axis = 0)
        else:
                total = weight_lists[0]
        ntile = np.sort(total)[-int(round(len(total) * significance))]
        return np.array([(0 if (x < ntile) else x) for x in total])

def powerset(iterable):
        """Output partial (i.e., nontrivial) powerset of list of classes/labels

        Args:
          iterable (iterable): list of classes/labels

        Returns: partial powerset
        """
        s = list(iterable)
        return chain.from_iterable(combinations(s, r) for r in range(2, len(s) + 1))

def sanitize_mask_name(string):
        """Extract class/label names from list-like literal string
        
        Args:
          string (string)

        Returns: (string) string without square brackets
        """
        return re.sub("[' \[\]]", '', string)

# percentage of voxels with non-zero weights
def non_empty_weights_proportion(weights):
        """Compute proportion of nonzero elements in array

        Args:
          weights (ndarray)

        Returns: (float)
        """
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
    np.savetxt(path + ".csv", pattern, header = ' '.join(labels), fmt='%10.5f')
    counts = np.unique(labels, return_counts = True)
    classes = counts[0]
    classes_counts = counts[1]
    plt.matshow(pattern)
    plt.xticks(np.cumsum(classes_counts), classes)
    plt.yticks(np.cumsum(classes_counts), classes)
    plt.colorbar()
    plt.savefig(path + ".svg")
    plt.close()
