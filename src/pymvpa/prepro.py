#!/usr/bin/python2

# author: Isaac David <isacdaavid@at@isacdaavid@dot@info>
# license: GPLv3 or later

from mvpa2.tutorial_suite import *
import sys

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

datasets = []
for filename in sys.argv[2:]:
	datasets.append(prepro(fmri_dataset(filename)))

merged = vstack(datasets, a = 'drop_nonunique')
merged.a.mapper = datasets[0].a.mapper[0]
merged.a.voxel_eldim = datasets[0].a.voxel_eldim
merged.a.voxel_dim = datasets[0].a.voxel_dim
merged.a.imghdr = datasets[0].a.imghdr
merged.a.imgtype = datasets[0].a.imgtype
merged.a.imghdr = datasets[0].a.imghdr

nimg = map2nifti(merged, merged)
nimg.to_filename(sys.argv[1])

