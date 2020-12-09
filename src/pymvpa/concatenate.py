#!/usr/bin/python2

# author: Isaac David <isacdaavid@at@isacdaavid@dot@info>
# license: GPLv3 or later

import sys
# sys.path.append('/home/inb/lconcha/fmrilab_software/miniconda2/lib/python2.7/site-packages')
from mvpa2.suite import *

datasets = []
for filename in sys.argv[2:]:
	datasets.append(fmri_dataset(filename))

merged = vstack(datasets, a = 'drop_nonunique')
merged.a.mapper = datasets[0].a.mapper
merged.a.voxel_eldim = datasets[0].a.voxel_eldim
merged.a.voxel_dim = datasets[0].a.voxel_dim
merged.a.imghdr = datasets[0].a.imghdr
merged.a.imgtype = datasets[0].a.imgtype

nimg = map2nifti(merged, merged)
nimg.to_filename(sys.argv[1])
