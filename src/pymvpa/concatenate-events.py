#!/usr/bin/python2

# author: Isaac David <isacdaavid@at@isacdaavid@dot@info>
# license: GPLv3 or later

from mvpa2.tutorial_suite import *
import sys

TR = 2000 # ms
NVOLS = 260 # number of volumes per fMRI sequence
HEADER = ['onset_time', 'age', 'sex', 'handedness', 'block', 'visual', 'face',
          'face_gender', 'emotion', 'gaze', 'target', 'response']

def adjust_onset_times(attr, tr, elapsed_nvols = 0):
        onset_times = np.array(attr.onset_time) - attr.onset_time[0] + \
                (tr * elapsed_nvols)
        attr.pop('onset_time')
        attr['onset_time'] = np.ndarray.tolist(onset_times)

def adjust_block_count(attr, elapsed_blocks):
	blocks = np.array(attr.block) + elapsed_blocks
	attr.pop('block')
	attr['block'] = np.ndarray.tolist(blocks)

# read attribute files from sys.argv[2:] and adjust to follow each other
datasets = []
for i in range(0, len(sys.argv[2:])):
	datasets.append(SampleAttributes(sys.argv[i + 2], header = HEADER))
	adjust_onset_times(datasets[i], TR, NVOLS * i)
	if i > 0:
		adjust_block_count(datasets[i], datasets[i - 1]['block'][-1])

# concatenate attribute objects
attr = SampleAttributes(sys.argv[2], header = HEADER)
for at in attr.keys():
	series = []
	for ds in datasets:
		series += ds[at] # concatenate to series' current ending
	attr.pop(at)
	attr[at] = series

# write concatenated and adjusted file into sys.argv[1]
fo = open(sys.argv[1], "w+")
for row in range(0, len(attr[HEADER[0]])):
	for at in HEADER:
		fo.write(str(attr[at][row]))
		fo.write('\t') if at != HEADER[-1] else True
	fo.write('\n')
fo.close()

