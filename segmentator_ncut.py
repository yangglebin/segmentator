#!/usr/bin/env python

# Part of the Segmentator library
# Copyright (C) 2016  Omer Faruk Gulban and Marian Schneider
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from nibabel import load
from matplotlib.colors import LogNorm
from matplotlib.widgets import Slider, Button
from utils import Ima2VolHistMapping
from segmentator_functions import responsiveObj

# %%
"""Load Data"""
#
import segmentator
nii = load(segmentator.args.filename)
#
#nii = load('/media/sf_D_DRIVE/Segmentator/Segmentator_local/#nii/sub-11_T1wDivPD_bet_cdfMatch.nii.gz')

ncut_labels = np.load(segmentator.args.ncut)
#ncut_labels = np.load('/media/sf_D_DRIVE/Segmentator/Segmentator_local/nii/multi-level_ncut_output_6recursion_S11.npy')

# transpose the labels
ncut_labels = np.transpose(ncut_labels, (1, 0, 2))
nrTotal_labels = sum([2**x for x in range(ncut_labels.shape[2])])
total_labels = np.arange(nrTotal_labels)
total_labels[1::2] = total_labels[-2:0:-2]

# relabel the labels from ncut, assign ascending integers starting with counter
counter = 0
for ind in np.arange(ncut_labels.shape[2]):
    tmp = ncut_labels[:, :, ind]
    uniqueVals = np.unique(tmp)
    nrUniqueVals = len(uniqueVals)
    newInd = np.arange(nrUniqueVals) + counter
    newVals = total_labels[newInd]
    tmp2 = np.zeros((tmp.shape))
    for ind2, val in enumerate(uniqueVals):
        tmp2[tmp == val] = newVals[ind2]
    counter = counter + nrUniqueVals
    ncut_labels[:, :, ind] = tmp2
lMax = np.max(ncut_labels)

orig_ncut_labels = ncut_labels
ima_ncut_labels = ncut_labels.copy()


# %%
"""Data Pre-Processing"""
orig = np.squeeze(nii.get_data())
# truncate too low and too high values
percDataMin = np.percentile(orig, 0.01)
orig[orig < percDataMin] = percDataMin
percDataMax = np.percentile(orig, 99.9)
orig[orig > percDataMax] = percDataMax
# auto-scaling for faster interface (0-500 or 600 seems fine)
scaleFactor = 500
orig = orig - orig.min()
orig = scaleFactor/orig.max() * orig
# define dataMin and dataMax for later use
dataMin = np.round(orig.min())
dataMax = np.round(orig.max())

# copy intensity data so we can flatten the copy and leave original intact
ima = orig.copy()
# calculate gradient magnitude (using L2 norm of the vector)
gra = np.gradient(ima)
gra = np.sqrt(np.power(gra[0], 2) + np.power(gra[1], 2) + np.power(gra[2], 2))

# reshape ima (more intuitive for voxel-wise operations)
ima = np.ndarray.flatten(ima)
gra = np.ndarray.flatten(gra)


# fill ncut_labels up with zeros if needed
binEdges = np.arange(dataMin, dataMax+1)
nrBins = len(binEdges)-1

if ncut_labels.shape[0] < nrBins:
    dif1 = nrBins - ncut_labels.shape[0]
    ncut_labels = np.append(ncut_labels,
                            np.zeros((dif1,
                                      ncut_labels.shape[1],
                                      ncut_labels.shape[2])),
                            axis=0)

if ncut_labels.shape[1] < nrBins:
    dif2 = nrBins - ncut_labels.shape[1]
    ncut_labels = np.append(ncut_labels,
                            np.zeros((ncut_labels.shape[0],
                                      dif2,
                                      ncut_labels.shape[2])),
                            axis=1)
# set number of free labels that the user can set
varNumAddLabel = 5

# %%
"""Plots"""
# Plot 2D histogram
fig = plt.figure()
ax = fig.add_subplot(121)
binEdges = np.arange(dataMin, dataMax+1)
nrBins = len(binEdges)-1
counts, xedges, yedges, volHistH = plt.hist2d(ima, gra,
                                              bins=binEdges,
                                              cmap='Greys'
                                              )

ax.set_xlim(dataMin, dataMax)
ax.set_ylim(0, dataMax)
ax.set_xlabel("Intensity f(x)")
ax.set_ylabel("Gradient Magnitude f'(x)")
ax.set_title("2D Histogram")

# plot colorbar for 2d hist
volHistH.set_norm(LogNorm(vmax=1000))
plt.colorbar(volHistH)

# plot hist mask (with ncut labels)
volHistMask = np.squeeze(ncut_labels[:, :, 0])
volHistMaskH = ax.imshow(volHistMask, interpolation='none',
                         alpha=0.2, cmap=plt.cm.gist_rainbow,
                         vmin=np.min(ncut_labels),
                         vmax=lMax+varNumAddLabel,
                         extent=[0, nrBins, nrBins, 0])

# plot 3D ima by default
ax2 = fig.add_subplot(122)
slcH = ax2.imshow(orig[:, :, int(orig.shape[2]/2)], cmap=plt.cm.gray,
                  vmin=ima.min(), vmax=ima.max(), interpolation='none',
                  extent=[0, orig.shape[1], orig.shape[0], 0])
imaMask = np.zeros(orig.shape[0:2])*total_labels[1]
imaMaskH = ax2.imshow(imaMask, interpolation='none', alpha=0.5,
                      cmap=plt.cm.gist_rainbow, vmin=np.min(ncut_labels),
                      vmax=lMax+varNumAddLabel,
                      extent=[0, orig.shape[1], orig.shape[0], 0])

# adjust subplots on figure
bottom = 0.30
fig.subplots_adjust(bottom=bottom)
plt.axis('off')


# %%
"""Initialisation"""
segmType = 'ncut'
# initiate a flexible figure object, pass to it usefull properties
flexFig = responsiveObj(figure=ax.figure,
                        axes=ax.axes,
                        axes2=ax2.axes,
                        segmType=segmType,
                        orig=orig,
                        nii=nii,
                        nrBins=nrBins,
                        sliceNr=int(0.5*orig.shape[2]),
                        slcH=slcH,
                        imaMask=imaMask,
                        imaMaskH=imaMaskH,
                        volHistMask=volHistMask,
                        volHistMaskH=volHistMaskH,
                        counterField=np.zeros((nrBins, nrBins)),
                        orig_ncut_labels=orig_ncut_labels,
                        ima_ncut_labels=ima_ncut_labels,
                        )

# make the figure responsive to clicks
flexFig.connect()
ima2volHistMap = Ima2VolHistMapping(xinput=ima, yinput=gra, binsArray=binEdges)
flexFig.invHistVolume = np.reshape(ima2volHistMap, orig.shape)

# %%
"""Sliders and Buttons"""
# colorbar slider
axcolor = 'lightgoldenrodyellow'
axHistC = plt.axes([0.15, bottom-0.15, 0.25, 0.025], axisbg=axcolor)
flexFig.sHistC = Slider(axHistC, 'Colorbar', 1, 5, valinit=3, valfmt='%0.1f')

# ima browser slider
axSliceNr = plt.axes([0.6, bottom-0.15, 0.25, 0.025], axisbg=axcolor)
flexFig.sSliceNr = Slider(axSliceNr, 'Slice', 0, 0.999,
                          valinit=0.5, valfmt='%0.3f')

# label slider
axLabels = plt.axes([0.15, bottom-0.2, 0.25, 0.025], axisbg=axcolor)
flexFig.sLabelNr = Slider(axLabels, 'Labels', 0, lMax+varNumAddLabel,
                          valinit=lMax, valfmt='%i')

# cycle button
cycleax = plt.axes([0.55, bottom-0.285, 0.075, 0.075])
flexFig.bCycle = Button(cycleax, 'Cycle\nView',
                        color=axcolor, hovercolor='0.975')
flexFig.cycleCount = 0

# export nii button
exportax = plt.axes([0.75, bottom-0.285, 0.075, 0.075])
flexFig.bExport = Button(exportax, 'Export\nNifti',
                         color=axcolor, hovercolor='0.975')

# export nyp button
exportax = plt.axes([0.85, bottom-0.285, 0.075, 0.075])
flexFig.bExportNyp = Button(exportax, 'Export\nNyp',
                            color=axcolor, hovercolor='0.975')

# reset button
resetax = plt.axes([0.65, bottom-0.285, 0.075, 0.075])
flexFig.bReset = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')

# %%
"""Updates"""
flexFig.sHistC.on_changed(flexFig.updateColorBar)
flexFig.sSliceNr.on_changed(flexFig.updateImaBrowser)
flexFig.sLabelNr.on_changed(flexFig.updateLabels)
flexFig.bCycle.on_clicked(flexFig.cycleView)
flexFig.bExport.on_clicked(flexFig.exportNifti)
flexFig.bExportNyp.on_clicked(flexFig.exportNyp)
flexFig.bReset.on_clicked(flexFig.resetGlobal)

plt.show()