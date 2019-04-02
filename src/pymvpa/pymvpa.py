#!/usr/bin/python2

# author: Isaac David <isacdaavid@at@isacdaavid@dot@info>
# license: GPLv3 or later

from mvpa2.tutorial_suite import *
import matplotlib.pyplot as plt

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

# load eprime events/design matrix (aka target attributes)
attr_fname = 'out/eprime/670/A.txt.csv'
attr = SampleAttributes(attr_fname,
                        header = ['onset_time', 'age', 'sex', 'handedness',
                                  'block', 'visual', 'face', 'face_gender',
	                          'emotion', 'gaze', 'target', 'response'])

# loading fmri data
bold_fname = 'data/feat/670/scans/5-fMRI_GazeCueing_1.feat/filtered_func_data_brain.nii.gz'

result_dist = []
for delay in range(9600, 9601, 1):
        ds = fmri_dataset(bold_fname)
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

plt.plot(result_dist)
plt.hist(result_dist, bins=1000)

# best result
print(results.samples)
print(validator.ca.stats.as_string(description = True))
validator.ca.stats.plot()
plt.show()

# sensibility map
weights = sensibility_map(model, ds4, 1)
# percentage of voxels with non-zero weights
print(len(weights[weights != 0.0]) / float(len(weights)))
# distribution of non-zero weights, normalized to the maximum weight
plt.hist(weights[weights != 0] / max(weights), bins=50)

# export sensibility map
nimg = map2nifti(ds, weights) # use ds.a.mapper to reverse flattening
nimg.to_filename('./weights-nn.nii.gz')
