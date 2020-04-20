#!/usr/bin/python2

# author: Isaac David <isacdaavid@at@isacdaavid@dot@info>
# license: GPLv3 or later

import sys
sys.path.append('/home/inb/lconcha/fmrilab_software/miniconda2/lib/python2.7/site-packages')
from mvpa2.suite import *

VOLS_PER_ACQ = 185

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

ds = fmri_dataset(sys.argv[2])
datasets = []
for i in range(0, ds.a.imghdr['dim'][4] / VOLS_PER_ACQ):
        datasets.append(prepro(ds[i * VOLS_PER_ACQ:(i + 1) * VOLS_PER_ACQ]))

merged = vstack(datasets, a = 'drop_nonunique')
merged.a.mapper = ds.a.mapper
merged.a.voxel_eldim = ds.a.voxel_eldim
merged.a.voxel_dim = ds.a.voxel_dim
merged.a.imghdr = ds.a.imghdr
merged.a.imgtype = ds.a.imgtype

nimg = map2nifti(ds, merged)
nimg.to_filename(sys.argv[1])
