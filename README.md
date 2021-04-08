# NP-MVPA: Non-programmers' Multi-Voxel Pattern Analysis

Collection of scripts to deploy a vanilla-yet-sound _classification_
and _RSA_ analysis pipeline for fMRI datasets, built on top of
[pyMVPA](http://www.pymvpa.org/) and R.

You should be able to adapt it with
minimal effort to any task-based (block paradigm) fMRI dataset, _no programming required_.
Think of it as a sort of barebones FSL's Feat substitute: from DICOM ->
NIFTI conversion and preprocessing all the way to group-level,
whole-brain or ROI activation maps.

Features:

 - [BIDS](https://bids.neuroimaging.io/)-compatible! (task TSVs coming soon)
 - flexible k-fold cross-validation
 - p-vals and effect sizes: rank-based hypothesis testing using Montecarlo permutations
 - simple Representational Symilarity Analysis (Pearson correlation and Euclidean distance matrices)
 - only linear SVM classifiers are currently supported
 - pretty plots (first-level and higher-level stats) and NIFTI sensitivity maps
 - simple CLI, because less is more
 - easy to parallelize for multiple subjects (multi-core thanks to `make`, optionally scale to SGE cluster for critical parts)
 - smart build-dependency tracking to resume processing where you left it, thanks to `make`

Use wisely! Automation is not a substitute for methodological and
statistical understanding.

## Getting Started

The following binaries are expected from your `$PATH`:

- typical GNU environment (`make`, `bash`, coreutils, `grep`, `sed`, `find`)
- `unzip`
- `dcm2niix` (optional, or just bring your own Nifti images)
- FSL:

   - `fslmaths`
   - `fsledithd`
   - `feat`
   - `flirt`
   - `fslstats`

- `python2`
- `Rscript`

Needed libraries:

- python2:

   - pyMVPA (tested with pyMVPA 2.5.0)
   - matplotlib

- R:

   - ggplot2
   - reshape2 (`acast()`)
   - effsize (`cohen.d()`)
   - [rainclouds](https://github.com/RainCloudPlots/RainCloudPlots)
     may ask for additional packages

## Running

1. Place your BIDS dataset inside `data/bids/` and fill in `TASKNAME`
   in the `Makefile`

2. obtain gray-matter segmentation masks (currently depends on SaaSS
   [volBrain](https://www.volbrain.upv.es/))

   ```
   make [ -j #CORES ] volbrain_tree
   ```

   process T1w images with volBrain and download ZIPs to
   `data/volbrain/$ID/`. Afterwards:

   ```
   make [ -j #CORES ] volbrain_unzip
   ```

   De-skull with:

   ```
   make [ -j #CORES ] t1w_brain_extraction
   ```

3. Preprocess with FSL (high-pass filter, slice-timing,
   motion-correction, fMRI <-> T1 registration):

   ```
   make [ -j #CORES ] concatenate
   make feat_prepro
   ```

4. produce gray-matter masks (volBrain T1 -> fMRI) and apply them to
   fMRI data

   ```
   make [ -j #CORES ] feat_masks
   make [ -j #CORES ] feat_brains
   ```

5. More prepro: adjust for the number of fMRI `RUNS` inside the
   `Makefile` (this currently **assumes all runs have the same number
   of volumes**), then linearly detrend and normalize signals:

   ```
   make [ -j #CORES ] detrend_normalize
   ```

6. specify the task layout before supervised learning can occur

   Add one `$ID.csv` per subject at `data/psychopy/` according to the
   following format (or just create symlinks to a single CSV if all
   subjects performed the same task):

   | block start time (ms) | label name | cross-validation fold # |
   |-|-|-|

7. Then specify which contrasts to try out from the label names. These
   are stored in `src/pymvpa/contrasts/`. You _must_ create at least
   one contrast category file there. See provided examples.

8. Train and validate classification models. Estimate accuracy
   distribution under the null hypothesis (currently set inside `Makefile` to
   5000 permutations and fsl_sub/SGE parallel processing cluster):

   __Note__: the 5000 Montecarlo simulations can take long!, but you still want
   at least a few thousand of them for accurate results.

   ```
   make pymvpa_whole
   ```

    First-level results will go to `out/pymvpa/`

9. transform resulting sensitivity maps back to T1w space, then
    denoise to improve spatial detection at group analysis:

    ```
    make [ -j #CORES ] register_results
    ```

10. Ask R to perform group-level inference and compose pretty plots
    (one test per contrast). Note this is currently set to run serially

    ```
    make poststats
    make group_level_mvpa
    ```

    Results will go to `out/poststats/`

11. (optional) make per-category summary CSVs with mean p-values and
    effect sizes:

    ```
    make postpoststats
    ```

## Citation

Please cite this paper when using the software:
<https://doi.org/10.1101/2021.04.04.438425>

## License

GNU General Public License (GPL) version 3, or (at your option) any
later version.
