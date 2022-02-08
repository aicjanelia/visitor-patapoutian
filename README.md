# Ardem Patapoutian -  iPALM

**Project:** Visualizing mechanically activated PIEZO channels

**Dates:** 2021-10-25 - 2021-11-05

**Hypothesis:** Inherently mechanosensitive mPIEZO1 expand and flatten upon channel activation and opening

**Aims:**
1. Determine optimal PIEZO1 constructs for iPALM imaging/Measure conformational parameters of mPIEZO1 without application of Yoda1 (closed	state).
2. Measure conformational parameters of mPIEZO1 with application of Yoda1 (open state).
3. Measure conformational parameters of mPIEZO1 with and without mechanical force (long-term goal).	

| Fluorophore | Target | Notes |
| :----: | :----: | :----: |
Alexa-647 | N-terminus of PIEZO1 | 'blades'
tdEOs | C-terminus of PIEZO1 | 'central pore'

## Analysis Workflow
The current analysis tools are designed to work in tandem with the software PeakSelector, which allows for filtering of iPALM localizations and pixelated image rendering of the underlying data. The bead removal approach can be used separately from any downstream analysis to clean up localization data for better visualization in PeakSelector. The candidate PIEZO segmentation approach looks for three nearby peaks in a pixelated/rendered image of the data, and uses this information to segment the underlying localization data into candidate PIEZO 'particles.'

**Example Data:** Initial testing of the analysis pipeline was conducted on the file _Run1-561_c123_sum_X14_processed_overlay_Fiducial_transform_complete_IDL.sav_ from the folder _21.1102-5/Run1-561/_

**Acknowledgments:** The subfunctions `bpass`, `cntrd`, and `pkfnd` were adapted for MATLAB by [Daniel Blair and Eric Dufresne](https://site.physics.georgetown.edu/matlab/code.html) from the IDL Particle Tracking software developed by David Grier, John Crocker, and Eric Weeks.

### Bead Removal
Bead removal takes as input an image and txt file generated by PeakSelector and outputs a new txt file with localizations corresponding to beads removed. This txt file can be reloaded into PeakSelector for further exploration and visualization of the data.

To use this approach:
1. Load the data of interest into PeakSelector and ensure it is properly filtered. For example, set quaulity control settings such as sigma rtNPh < 0.06 or set bounds on the unwrapped z value.
2. Export the PeakSelector data as an ASCII file with the option of xy values in pixels.
3. Open the 'Cust. Tiff' menu and export the Total Raw Data as a tiff file with 133.33 nm per pixel.
4. Fill in appropriate filenames for the ASCII and tiff files under the USER PARAMETERS section of `beadRemoval_v0_anisotropic.m`.
5. Optionally adjust the parameters for bead removal:
    - `rRemoveX` and `rRemoveY` indicate the size of the region in pixels (133.33 nm/pixel) removed around each bead, in the X & Y directions respectively.
    - `rParticle` indicates the approximate size of beads in the total raw data image, in pixels.
    - `beadThresh` is a threshold for distinguishing beads from the background.
6. Run `beadRemoval_v0_anisotropic.m`. Figures will provide a check on the regions of data that were removed, and a new txt file will be generated without localizations from beads.

The bead-removed txt file can be reloaded into PeakSelctor through the menu option _Import User ASCII_ with xy coords in pixels, the box for headers checked, and the values 0 to 48 for columns (see the file _peakSelector_columnIndices_ for an easy list of values to copy and paste).

### Segmentation of Candidate PIEZOs
Segmentation of the candidate PIEZO 'particles' is based on the assumption the PIEZOs will appear in a rendered image as 3 neighboring peaks. Data is stored as a .mat file in a format that can be used for further processing a single particle averaging routine developed by [Heydarian, et. al.](https://github.com/imphys/smlm_datafusion3d). Data in the PeakSelector format is also included in the particle structure, and thus the data can also be pulled back into PeakSelector for further rendering.

To use this approach:
1. First run `beadRemoval_v0_anisotropic.m` (described in the Bead Removal section) or otherwise remove beads from the data.
2. Load the bead-removed localization data into PeakSelector.
3. Select the option: PeakSelector > SpecialFunctions > SwapZ with Unwrapped Z. (This makes sure any filtering on unwrapped Z is properly used in the rendering.)
4. Open the 'Cust. Tiff' window
    - Make sure "Filter" is set to "Frame Peaks"
    - Set "NM per Image Pixel" to 2
5. Click "Render." Note that this process could be slow due to the size of the image.
6. Click "Save TIFF float" and save the file to the appropriate directory.
7. Fill in the USER PARAMETERS for `piezoSegment_beadsRemoved.m`. The data directory and file basename are required. Other parameters should be held constant across images.
    - `dataDir` is the folder containing the txt data and rendered images.
    - `saveTag` is the basename for the files and will be used for saving figures and data.
    - `xyScale` is the scale of pixels in the txt file, which should be 133.33 nm for iPALM data.
    - `nmPix` is the scale of pixels in the rendered image, which should be 2 for clear visualization of peaks in the rendered image.
    - `rBlade` is the scale in pixels of a blade-blob in the rendered image.
    - `peakThresh` is a threshold on the bandpassed image used to find blade-blobs.
    - `rThresh` is used to exclude peaks that are too close together.
    - `neighborNumber` is the number of neighbors each candidate PIEZO blade should have (i.e., each particle in a group of 3 has 2 neighbors)
    - `neighborDist` is the maximum distance for two blobs to be considered neighbors.
    - `minDist` is the minimum distance for two blobs to be considered neighbors; this parameter gets rid of self comparisons and over-finding artifacts
    - `bord` is the size of the region around each candidate PIEZO to include in the localization segmentation.
8. Run `piezoSegment_beadsRemoved.m`. Outputs include sanity-check figures, a .mat file containing the entire workspace (for current troubleshooting purposes), and a particles.mat file containing the data for later averaging.
