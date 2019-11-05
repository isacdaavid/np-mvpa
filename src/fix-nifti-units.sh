#!/usr/bin/env bash

# author: Isaac David <isacdaavid@at@isacdaavid@dot@info>
# license: GPLv3 or later

# this script serves as an editor passed to `fsledithd`.
# pymvpa mangles the units and dimensions of certain generated Niftis
# https://github.com/PyMVPA/PyMVPA/blob/8a372b0ccc7fe8d3e398009222b8b702839c773f/mvpa2/datasets/mri.py#L311
# we overcome this manually using `fsledithd`.

sed -Ei "s/(d)(x|y|z)( = '1')/\1\2 = '4'/ ; \
         s/dt = '1'/dt = '2'\n  xyz_units = '2'\n  time_units = '8'/" \
        "$1"
