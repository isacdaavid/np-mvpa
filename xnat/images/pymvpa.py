#!/usr/bin/python3

# author: Isaac David <isacdaavid@at@isacdaavid@dot@info>
# license: GPLv3 or later

from mvpa2.tutorial_suite import *
import matplotlib.pyplot as plt

# load target attributes (eprime events/design matrix)
attr_fname = './670-A.csv'
attr = SampleAttributes(attr_fname,
                        header = ['onset_time', 'age', 'sex', 'handedness',
                                  'visual', 'face', 'face_gender', 'emotion',
                                  'gaze', 'target', 'response'])

# loading fmri data
bold_fname = './filtered_func_data-volbrain-mask.nii.gz'
ds = fmri_dataset(bold_fname)

################################################################################
# volume labeling
################################################################################

SLICE_TIMING_REFERENCE = +1000 # ms
OPTIMAL_HRF_DELAY = SLICE_TIMING_REFERENCE + 300 # ms

def label(ds, attr, slice_timing_reference, hrf_delay):
        # time-shift for eprime preambulus
        onset_times = np.array(attr.onset_time) - attr.onset_time[0]
        attr.pop('onset_time')
        attr['onset_time'] = np.ndarray.tolist(onset_times)

        # convert volume times to ms, as eprime, then time-slice-shift
        # (tr=2s, interleaved)
        ds.sa.time_coords = (ds.sa.time_coords * 1000) + slice_timing_reference

        # select attributes that correspond to closest volume in the
        # present or future
        indices = list()
        for vol_time in ds[ds.sa.time_coords >= hrf_delay].sa.time_coords:
                match_floor = max(onset_times[onset_times <= (vol_time
                                                              - hrf_delay)])
                indices.append(attr.onset_time.index(match_floor))

        while(len(indices) < len(ds)):
                indices.append(indices[-1])

        # add
        for at in attr.keys():
                ds.sa[at] = [attr[at][index] for index in indices]

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

# balanced subsampling of volumes
def subsample(ds):
        max_samples_per_cat = min(np.unique(ds[{'emotion': ['happy', 'sad', 'neutral']}].sa.emotion,
                                            return_counts = True)[1])
        dssub = ds[{'emotion': ['happy']}][:max_samples_per_cat]
        dssub = vstack((dssub, ds[{'emotion': ['sad']}][:max_samples_per_cat]))
        dssub = vstack((dssub, ds[{'emotion': ['neutral']}][:max_samples_per_cat]))
        return dssub

def train(ds):
        # please crossvalidation()'s stubborness
        ds.sa['chunks'] = [1] * len(ds.samples) # all samples belong to 1 chunk
        ds.sa['targets'] = ds.sa.emotion        # label attribute
        parity = [0, 1]  * (len(ds.samples) / 2) # partition attribute
        parity.append(0) if ((len(ds) / 3) % 2 == 1) else True
        ds.sa['parity'] = parity

        clf = LinearCSVMC()
        cvte = CrossValidation(clf, HalfPartitioner(attr='parity'),
                               errorfx = lambda p, t: np.mean(p == t),
                               enable_ca=['stats'])
        return clf,cvte



################################################################################
# sensibility analysis
################################################################################

# outputs an "activation" map (rather, a sensibility map)
def sensibility_map(model, ds, significance = .05):
        sensana = model.get_sensitivity_analyzer()
        sens = sensana(ds)
        weights = (sens[0].samples[0] +
                   sens[1].samples[0] +
                   sens[2].samples[0]) / 3.0

        # n most important weights
        ntile = np.sort(abs(weights))[-int(round(len(weights) * significance))]
        weights_sub = np.array([(0 if (x < ntile) else x) for x in abs(weights)])
        return weights_sub

################################################################################
# main
################################################################################

result_dist = []
for delay in range(OPTIMAL_HRF_DELAY, OPTIMAL_HRF_DELAY + 1, 1):
        ds = fmri_dataset(bold_fname)
        label(ds, attr, SLICE_TIMING_REFERENCE, delay)
        ds2 = prepro(ds)
        # nimg = map2nifti(ds, ds2) # use ds.mapper to preserve preprocessing
        # nimg.to_filename('./filtered_func_data_pymvpa.nii.gz')
        ds2 = subsample(ds2)
        model,validator = train(ds2)
        results = validator(ds2)

        # display results
        result_dist.append(np.mean(results))
        print('mean accuracy: ' + str(np.mean(results)))

plt.plot(result_dist)
plt.hist(result_dist, bins=1000)

# best result
print(results.samples)
print(validator.ca.stats.as_string(description = True))
validator.ca.stats.plot()
plt.show()

# sensibility map
weights = sensibility_map(model, ds2, 1)
# percentage of voxels with non-zero weights
print(len(weights[weights != 0.0]) / float(len(weights))) # .1676
# distribution of non-zero weights, normalized to the maximum weight
plt.hist(weights[weights != 0] / max(weights), bins=50)

# export sensibility map
nimg = map2nifti(ds, weights) # use ds.a.mapper to reverse flattening
nimg.to_filename('./weights-nn.nii.gz')
