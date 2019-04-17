#!/usr/bin/python2

# author: Isaac David <isacdaavid@at@isacdaavid@dot@info>
# license: GPLv3 or later

from mvpa2.tutorial_suite import *
import sys

TR = 2000 # ms
NVOLS = 260 # number of volumes per fMRI sequence
HEADER = ['onset_time', 'age', 'sex', 'handedness', 'block', 'visual', 'face',
          'face_gender', 'emotion', 'gaze', 'target', 'response']

def center_onset_times(attr, tr, elapsed_nvols = 0):
        onset_times = np.array(attr.onset_time) - attr.onset_time[0] + \
                (tr * elapsed_nvols)
        attr.pop('onset_time')
        attr['onset_time'] = np.ndarray.tolist(onset_times)

datasets = []
for i in range(0, len(sys.argv[2:])):
	datasets.append(SampleAttributes(sys.argv[i + 2], header = HEADER))
	center_onset_times(datasets[i], TR, NVOLS * i)

attr = SampleAttributes(sys.argv[2], header = HEADER)
for at in attr.keys():
	series = []
	for ds in datasets:
		series += ds[at] # concatenate to series' current ending
	attr.pop(at)
	attr[at] = series

fo = open(sys.argv[1], "w+")
for row in range(0, len(attr[HEADER[0]])):
	for at in HEADER:
		fo.write(str(attr[at][row]))
		fo.write('\t') if at != HEADER[-1] else True
	fo.write('\n')
fo.close()
