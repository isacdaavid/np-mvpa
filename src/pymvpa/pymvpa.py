#!/usr/bin/python2

# author: Isaac David <isacdaavid@at@isacdaavid@dot@info>
# license: GPLv3 or later

from mvpa2.tutorial_suite import *
import matplotlib.pyplot as plt
import sys

################################################################################
# volume labeling
################################################################################

SLICE_TIMING_REFERENCE = +1000 # ms
OPTIMAL_HRF_DELAY = SLICE_TIMING_REFERENCE + 300 # ms

def center_attr_onset_times(attr, tr, elapsed_nvols = 0):
        onset_times = np.array(attr.onset_time) - attr.onset_time[0] + \
                (tr * elapsed_nvols)
        attr.pop('onset_time')
        attr['onset_time'] = np.ndarray.tolist(onset_times)

def label(ds, attr, slice_timing_reference, hrf_delay):
        # time-shift for eprime preambulus
        onset_times = np.array(attr.onset_time) - attr.onset_time[0]
        attr.pop('onset_time')
        attr['onset_time'] = np.ndarray.tolist(onset_times)

        # convert volume times to ms, as eprime, then time-slice-shift
        # (tr=2s, interleaved)
        ds.sa.time_coords = (ds.sa.time_coords * 1000) + slice_timing_reference

        indices = list()
        # exclude volumes prior to 1 hrf_delay, since they are not supposed to
        # to reflect a maximum of task-related HRF activity
        for vol_time in ds[ds.sa.time_coords >= hrf_delay].sa.time_coords:
                # select latest event to occur before (vol_time - hrf_delay)
                match_floor = \
                        max(onset_times[onset_times <= (vol_time - hrf_delay)])
                indices.append(attr.onset_time.index(match_floor))

        # excluded premature volumes will be assigned the attributes of first
        # non-excluded volume
        useless_vols = 0
        while(len(indices) < len(ds)):
                indices.insert(0, indices[0])
                useless_vols = useless_vols + 1

        # label volumes according to attribute sublist
        for at in attr.keys():
                ds.sa[at] = [attr[at][index] for index in indices]

        return ds[useless_vols:]

################################################################################
# extra preprocessing
################################################################################

def prepro(ds):
        # detrending
        detrender = PolyDetrendMapper(polyord = 1)
        ds2 = ds.get_mapped(detrender)
        # poly_detrend(ds, polyord = 1) # in-place alternative

        # intensity normalisation
        zscorer = ZScoreMapper(chunks_attr = None)
        zscorer.train(ds2)
        ds3 = ds2.get_mapped(zscorer)
        # zscore(ds, chunks_attr = None) # in-place alternative

        return ds3

################################################################################
# classification
################################################################################

def subsample(ds0):
        ds1 = ds0[{'emotion': ['happy', 'neutral', 'sad']}]
        # sample from independent blocks so as to avoid temporally-correlated
        # volumes in both training and test partitions. take earliest volume
        indep = []
        for b in np.unique(ds1.sa.block):
                indep.append(min(ds1[{'block': [b]}].sa.time_indices))
                ds2 = ds1[{'time_indices': indep}]

        # good machines are emotionally-balanced machines ;)
        max_samples = min(np.unique(ds2[{'emotion': ['happy',
                                                     'sad',
                                                     'neutral']}].sa.emotion,
                                    return_counts = True)[1])
        ds3 = ds2[{'emotion': ['happy']}][:max_samples]
        ds3 = vstack((ds3, ds2[{'emotion': ['sad']}][:max_samples]))
        ds3 = vstack((ds3, ds2[{'emotion': ['neutral']}][:max_samples]))

        return ds3

def train(ds):
        ds.sa['targets'] = ds.sa.emotion # this is the label/target attribute
        clf = LinearCSVMC()
        cvte = CrossValidation(clf, NFoldPartitioner(attr = 'block'),
                               errorfx = lambda p, t: np.mean(p == t),
                               enable_ca=['stats'])
        return clf,cvte

################################################################################
# sensitivity analysis
################################################################################

# outputs an "activation" map (rather, a sensibility map)
def sensibility_map(model, ds):
        sensana = model.get_sensitivity_analyzer()
        sens = sensana(ds)
        for i in range(0, len(sens.targets)):
                sensmap = sens.targets[i]
                # emotional vs neutral
                if str(sensmap) == "('neutral', 'happy')" or \
                   str(sensmap) == "('happy', 'neutral')":
                   i1_1 = i
                if str(sensmap) == "('neutral', 'sad')" or \
                   str(sensmap) == "('sad', 'neutral')":
                   i1_2 = i
                # happy vs sad
                if str(sensmap) == "('happy', 'sad')" or \
                   str(sensmap) == "('sad', 'happy')":
                   i2 = i

        all_weights = (sens[0].samples[0] + sens[1].samples[0] +
                       sens[2].samples[0]) / 3.0
        emo_vs_neu = (sens[i1_1].samples[0] + sens[i1_2].samples[0]) / 2.0
        hap_vs_sad = sens[i2].samples[0]

        return all_weights,emo_vs_neu,hap_vs_sad

# remove sign, take n most significant weights
def normalize_weights(weights, significance = .05):
        ntile = np.sort(abs(weights))[-int(round(len(weights) * significance))]
        return np.array([(0 if (x < ntile) else x) for x in abs(weights)])

################################################################################
# main
################################################################################

# load eprime events/design matrix (aka target attributes)
ATTR_FNAME = sys.argv[1]
attr = SampleAttributes(ATTR_FNAME,
                        header = ['onset_time', 'age', 'sex', 'handedness',
                                  'block', 'visual', 'face', 'face_gender',
                                  'emotion', 'gaze', 'target', 'response'])

# loading fmri data
BOLD_FNAME = sys.argv[2]

STEP = 200 # ms
result_dist = []
for delay in range(0, 600, STEP):
        ds = fmri_dataset(BOLD_FNAME)
        ds2 = label(ds, attr, SLICE_TIMING_REFERENCE, delay)
        ds3 = prepro(ds2)

        # nimg = map2nifti(ds, ds3) # use ds.mapper to preserve preprocessing
        # nimg.to_filename('./filtered_func_data_pymvpa.nii.gz')

        ds4 = subsample(ds3)
        model,validator = train(ds4)
        results = validator(ds4)

        # display results
        result_dist.append(np.mean(results))
        print(str(ds4.nsamples / 3) + '\tmean accuracy: ' + str(np.mean(results)))

OUTDIR = sys.argv[3]

fo = open(OUTDIR + "/result-time-series.txt", "w+")
fo.writelines("%s\n" % item for item in result_dist)
fo.close()

plt.plot(result_dist)
plt.savefig(OUTDIR + '/result-time-series.svg')
plt.close()
plt.hist(result_dist, bins=1000)
plt.savefig(OUTDIR + '/result-dist.svg')
plt.close()

# best result
optimal_delay = result_dist.index(max(result_dist)) * STEP
ds = fmri_dataset(BOLD_FNAME)
ds2 = label(ds, attr, SLICE_TIMING_REFERENCE, optimal_delay)
ds3 = prepro(ds2)
ds4 = subsample(ds3)
model,validator = train(ds4)
results = validator(ds4)

print(results.samples)
print(validator.ca.stats.as_string(description = True))
validator.ca.stats.plot()
plt.savefig(OUTDIR + '/conf-matrix.svg')
plt.close()

# sensibility map
all_weights,emo_vs_neu,hap_vs_sad = sensibility_map(model, ds4)
# percentage of voxels with non-zero weights
print(len(all_weights[all_weights != 0.0]) / float(len(all_weights)))
all_weights = normalize_weights(all_weights, significance = 1)
emo_vs_neu = normalize_weights(emo_vs_neu, significance = 1)
hap_vs_sad = normalize_weights(hap_vs_sad, significance = 1)
# distribution of non-zero weights, normalized to the maximum one
plt.hist(all_weights[all_weights != 0] / max(all_weights), bins=50)
plt.savefig(OUTDIR + '/weights-dist.svg')
plt.close()

# export sensitivity maps
nimg = map2nifti(ds, all_weights) # use ds.a.mapper to reverse flattening
nimg.to_filename(OUTDIR + '/all-weights-nn.nii.gz')
nimg = map2nifti(ds, emo_vs_neu)
nimg.to_filename(OUTDIR + '/emo-vs-neu-weights-nn.nii.gz')
nimg = map2nifti(ds, hap_vs_sad)
nimg.to_filename(OUTDIR + '/hap_vs_sad-weights-nn.nii.gz')
