"""
Processing input and plotting.

Segmentator main.
"""

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
import config as cfg
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.widgets import Slider, Button, LassoSelector
from matplotlib import path
from nibabel import load
from segmentator.utils import map_ima_to_2D_hist, prep_2D_hist
from segmentator.utils import truncate_range, scale_range
from segmentator.utils import set_gradient_magnitude
from gui_utils import sector_mask, responsiveObj

"""Load Data"""
nii = load(cfg.filename)

#
"""Data Processing"""
orig = np.squeeze(nii.get_data())
dims = orig.shape
orig = truncate_range(orig, percMin=cfg.perc_min, percMax=cfg.perc_max)
orig = scale_range(orig, scale_factor=cfg.scale, delta=0.0001)
gra = set_gradient_magnitude(orig, cfg.gramag)

# reshape ima (more intuitive for voxel-wise operations)
ima = np.copy(orig.flatten())
gra = gra.flatten()

#
"""Plots"""
# Set up a colormap:
palette = plt.cm.Reds
palette.set_over('r', 1.0)
palette.set_under('w', 0)
palette.set_bad('m', 1.0)

# Plot 2D histogram
fig = plt.figure(facecolor='0.775')
ax = fig.add_subplot(121)

counts, volHistH, d_min, d_max, nr_bins, bin_edges \
    = prep_2D_hist(ima, gra, discard_zeros=cfg.discard_zeros)

ax.set_xlim(d_min, d_max)
ax.set_ylim(d_min, d_max)
ax.set_xlabel("Intensity f(x)")
ax.set_ylabel("Gradient Magnitude f'(x)")
ax.set_title("2D Histogram")

# plot colorbar for 2d hist
volHistH.set_norm(LogNorm(vmax=np.power(10, cfg.cbar_init)))
plt.colorbar(volHistH)

# plot 3D ima by default
ax2 = fig.add_subplot(122)
sliceNr = int(0.5*dims[2])
imaSlcH = ax2.imshow(orig[:, :, sliceNr], cmap=plt.cm.gray, vmin=ima.min(),
                     vmax=ima.max(), interpolation='none',
                     extent=[0, dims[1], dims[0], 0])

imaSlcMsk = np.ones(dims[0:2])
imaSlcMskH = ax2.imshow(imaSlcMsk, cmap=palette, vmin=0.1,
                        interpolation='none', alpha=0.5,
                        extent=[0, dims[1], dims[0], 0])

# adjust subplots on figure
bottom = 0.30
fig.subplots_adjust(bottom=bottom)
fig.canvas.set_window_title(nii.get_filename())
plt.axis('off')

#
"""Initialisation"""
# create first instance of sector mask
sectorObj = sector_mask((nr_bins, nr_bins), cfg.init_centre, cfg.init_radius,
                        cfg.init_theta)

# draw sector mask for the first time
volHistMaskH, volHistMask = sectorObj.draw(ax, cmap='Reds', alpha=0.2,
                                           vmin=0.1, interpolation='nearest',
                                           origin='lower',
                                           extent=[0, nr_bins, 0, nr_bins])

# initiate a flexible figure object, pass to it usefull properties
idxLasso = np.zeros(nr_bins*nr_bins, dtype=bool)
lassoSwitchCount = 0
flexFig = responsiveObj(figure=ax.figure, axes=ax.axes, axes2=ax2.axes,
                        segmType='main', orig=orig, nii=nii,
                        sectorObj=sectorObj,
                        nrBins=nr_bins,
                        sliceNr=sliceNr,
                        imaSlcH=imaSlcH,
                        imaSlcMsk=imaSlcMsk, imaSlcMskH=imaSlcMskH,
                        volHistMask=volHistMask, volHistMaskH=volHistMaskH,
                        contains=volHistMaskH.contains,
                        counts=counts,
                        idxLasso=idxLasso,
                        initTpl=(cfg.perc_min, cfg.perc_max, cfg.scale),
                        lassoSwitchCount=lassoSwitchCount)

# make the figure responsive to clicks
flexFig.connect()
ima2volHistMap = map_ima_to_2D_hist(xinput=ima, yinput=gra, bins_arr=bin_edges)
flexFig.invHistVolume = np.reshape(ima2volHistMap, dims)

#
"""Sliders and Buttons"""
# colorbar slider
axcolor, hovcolor = '0.875', '0.975'
axHistC = plt.axes([0.15, bottom-0.20, 0.25, 0.025], facecolor=axcolor)
flexFig.sHistC = Slider(axHistC, 'Colorbar', 1, cfg.cbar_max,
                        valinit=cfg.cbar_init, valfmt='%0.1f')

# ima browser slider
axSliceNr = plt.axes([0.6, bottom-0.15, 0.25, 0.025], facecolor=axcolor)
flexFig.sSliceNr = Slider(axSliceNr, 'Slice', 0, 0.999, valinit=0.5,
                          valfmt='%0.2f')

# theta sliders
aThetaMin = plt.axes([0.15, bottom-0.10, 0.25, 0.025], facecolor=axcolor)
flexFig.sThetaMin = Slider(aThetaMin, 'ThetaMin', 0, 359.9,
                           valinit=cfg.init_theta[0], valfmt='%0.1f')
aThetaMax = plt.axes([0.15, bottom-0.15, 0.25, 0.025], facecolor=axcolor)
flexFig.sThetaMax = Slider(aThetaMax, 'ThetaMax', 0, 359.9,
                           valinit=cfg.init_theta[1]-0.1, valfmt='%0.1f')

# cycle button
cycleax = plt.axes([0.55, bottom-0.2475, 0.075, 0.0375])
flexFig.bCycle = Button(cycleax, 'Cycle',
                        color=axcolor, hovercolor=hovcolor)

# rotate button
rotateax = plt.axes([0.55, bottom-0.285, 0.075, 0.0375])
flexFig.bRotate = Button(rotateax, 'Rotate',
                         color=axcolor, hovercolor=hovcolor)

# export nii button
exportax = plt.axes([0.75, bottom-0.285, 0.075, 0.075])
flexFig.bExport = Button(exportax, 'Export\nNifti',
                         color=axcolor, hovercolor=hovcolor)

# export nyp button
exportax = plt.axes([0.85, bottom-0.285, 0.075, 0.075])
flexFig.bExportNyp = Button(exportax, 'Export\nHist',
                            color=axcolor, hovercolor=hovcolor)

# reset button
resetax = plt.axes([0.65, bottom-0.285, 0.075, 0.075])
flexFig.bReset = Button(resetax, 'Reset', color=axcolor, hovercolor=hovcolor)


#
"""Updates"""
flexFig.sHistC.on_changed(flexFig.updateColorBar)
flexFig.sSliceNr.on_changed(flexFig.updateImaBrowser)
flexFig.sThetaMin.on_changed(flexFig.updateThetaMin)
flexFig.sThetaMax.on_changed(flexFig.updateThetaMax)
flexFig.bCycle.on_clicked(flexFig.cycleView)
flexFig.bRotate.on_clicked(flexFig.changeRotation)
flexFig.bExport.on_clicked(flexFig.exportNifti)
flexFig.bExportNyp.on_clicked(flexFig.exportNyp)
flexFig.bReset.on_clicked(flexFig.resetGlobal)


#
"""New stuff: Lasso (Experimental)"""
# Lasso button
lassoax = plt.axes([0.15, bottom-0.285, 0.075, 0.075])
bLasso = Button(lassoax, 'Lasso\nOff', color=axcolor, hovercolor=hovcolor)


def lassoSwitch(event):
    """Enable disable lasso tool."""
    global lasso
    lasso = []
    flexFig.lassoSwitchCount = (flexFig.lassoSwitchCount+1) % 2
    if flexFig.lassoSwitchCount == 1:  # enable lasso
        flexFig.disconnect()  # disable drag function of sector mask
        lasso = LassoSelector(ax, onselect)
        bLasso.label.set_text("Lasso\nOn")
    else:  # disable lasso
        # lasso = []  # I am not sure we want to reset lasso with this button
        flexFig.connect()  # enable drag function of sector mask
        bLasso.label.set_text("Lasso\nOff")


# Pixel coordinates
pix = np.arange(nr_bins)
xv, yv = np.meshgrid(pix, pix)
pix = np.vstack((xv.flatten(), yv.flatten())).T


def onselect(verts):
    """Lasso related."""
    global pix
    p = path.Path(verts)
    newLasIdx = p.contains_points(pix, radius=1.5)  # new lasso indices
    flexFig.idxLasso[newLasIdx] = True  # updated old lasso indices
    # update volume histogram mask
    flexFig.remapMsks()
    flexFig.updatePanels(update_slice=False, update_rotation=True,
                         update_extent=True)


bLasso.on_clicked(lassoSwitch)
flexFig.remapMsks()
flexFig.updatePanels(update_slice=True, update_rotation=False,
                     update_extent=False)
plt.show()
