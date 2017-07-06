# nmr
Python Quadrupolar NMR Visualiser

This software simulates NMR spectra using peak and splitting data. Cauchy-Lorentz distributions are generated for each peak which are then convolved against a splitting pattern calculated using spin, nuclei and J coupling.

## Usage

### Machine

Operating Frequency: This is the frequency (in MHz) that the machine runs on.

Noise: This will generate a small amount of white noise to add to the plot.

Resolution: The space between each calculated point.

Range: The region of the NMR plot to simulate.

### Peaks

Add Peak: Append a peak to the plot. By default, this peak will be from one nucleus, have a chemical shift of 7, and a half-width at half-maximum of 0.01 ppm.

Remove Peak: Removes the current peak (and all its splittings).

Peaks: A dropdown menu of all the peaks in the plot. Use this to change the current peak.

Shift: The chemical shift of the current peak.

Nuclei: The number of nuclei that contribute to the peak.

Half Width Half Maximum: The HWHM of the peak (from 1 nucleus and when there is no splitting).

### Splitting Nuclei

Add Splitting: Attach a splitting pattern to the current peak. By default, the splitting pattern will consist of 1 nucleus with half spin, with a J coupling value of 20Hz and 100% abundace. In other words, this is a doublet.

Remove Splitting: Remove the current splitting from the current peak.

J Coupling: The amount in Hz that the splitting pattern will split the current peak. Corresponds to the machine frequency divided by J value in ppm.

Nuclei: The number of nuclei that create the splitting pattern.

Spin: The spin of all the nuclei in the splitting pattern.

Abundance: The isotopic abundance of all the nuclei in the splitting pattern. Used to simulate satellites.

### Controls

Update: Update the plot.

From .log: Load data from a Gaussian log file. Note that this file must include the keywords NMR=(GIAO,spinspin). See below.

### Loading from a Gaussian log file.

Provided the calculation has been performed with the keywords NMR(GIAO,spinspin), data can be imported from a Gaussian log file. Once a file is chosen, the following options will appear:

Select Element: Select the element that is used for the plot.

Reference Shift: The chemical shift from the reference peak in ppm. This can be from the default list in Gaussian, or calculated.

Degeneracy Threshold: When J coupling values are imported, any that differ by less than this amount will be merged (effectively increasing the number of nuclei in a splitting pattern. This saves time when generating the spectrum.

Decouple Elements: If true, only homonuclear coupling is considered.

Go: Import the data to the GUI.

