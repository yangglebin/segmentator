"""Functions covering the user interaction with the GUI."""

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
import os
import numpy as np
import matplotlib.pyplot as plt
from utils import map_2D_hist_to_ima
from nibabel import save, Nifti1Image
import config as cfg


class responsiveObj:
    """Stuff to interact in the user interface."""

    def __init__(self, **kwargs):
        if kwargs is not None:
            for key, value in kwargs.iteritems():
                setattr(self, key, value)
        self.basename = self.nii.get_filename().split(os.extsep, 1)[0]
        self.press = None
        self.ctrlHeld = False
        self.labelNr = 0
        self.imaSlcMskSwitch, self.volHistHighlightSwitch = 0, 0
        self.TranspVal = 0.5
        self.nrExports = 0
        self.entropWin = 0
        self.borderSwitch = 0
        self.imaSlc = self.orig[:, :, self.sliceNr]  # selected slice
        self.cycleCount = 0
        self.cycRotHistory = [[0, 0], [0, 0], [0, 0]]
        self.highlights = [[], []]  # to hold image to histogram circles

    def remapMsks(self, remap_slice=True):
        """Update volume histogram to image mapping.

        Parameters
        ----------
        remap_slice : bool
            Do histogram to image mapping. Used to map displayed slice mask.

        """
        if self.segmType == 'main':
            self.volHistMask = self.sectorObj.binaryMask()
            self.volHistMask = self.lassoArr(self.volHistMask,
                                             self.idxLasso)
            self.volHistMaskH.set_data(self.volHistMask)
        elif self.segmType == 'ncut':
            self.labelContours()
            self.volHistMaskH.set_data(self.volHistMask)
            self.volHistMaskH.set_extent((0, self.nrBins, self.nrBins, 0))
        # histogram to image mapping
        if remap_slice:
            self.imaSlcMsk = map_2D_hist_to_ima(
                self.invHistVolume[:, :, self.sliceNr], self.volHistMask)
            if self.borderSwitch == 1:
                self.imaSlcMsk = self.calcImaMaskBrd()

    def updatePanels(self, update_slice=True, update_rotation=False,
                     update_extent=False):
        """Update histogram and image panels."""
        if update_rotation:
            self.checkRotation()
        if update_extent:
            self.updateImaExtent()
        if update_slice:
            self.imaSlcH.set_data(self.imaSlc)
        self.imaSlcMskH.set_data(self.imaSlcMsk)
        self.figure.canvas.draw()

    def connect(self):
        """Make the object responsive."""
        self.cidpress = self.figure.canvas.mpl_connect(
            'button_press_event', self.on_press)
        self.cidrelease = self.figure.canvas.mpl_connect(
            'button_release_event', self.on_release)
        self.cidmotion = self.figure.canvas.mpl_connect(
            'motion_notify_event', self.on_motion)
        self.cidkeypress = self.figure.canvas.mpl_connect(
            'key_press_event', self.on_key_press)
        self.cidkeyrelease = self.figure.canvas.mpl_connect(
            'key_release_event', self.on_key_release)

    def on_key_press(self, event):
        """Determine what happens if key is pressed."""
        if event.key == 'control':
            self.ctrlHeld = True
        elif event.key == 'q':
            self.imaSlcMskIncr(-0.1)
        elif event.key == 'w':
            self.imaSlcMskTransSwitch()
        elif event.key == 'h':
            self.volHistHighlightTransSwitch()
        elif event.key == 'e':
            self.imaSlcMskIncr(0.1)
        elif event.key == '1':
            self.borderSwitch = (self.borderSwitch + 1) % 2
            self.remapMsks()
            self.updatePanels(update_slice=False, update_rotation=False,
                              update_extent=False)

        if self.segmType == 'main':
            if event.key == 'up':
                self.sectorObj.scale_r(1.05)
                self.remapMsks()
                self.updatePanels(update_slice=False, update_rotation=True,
                                  update_extent=False)
            elif event.key == 'down':
                self.sectorObj.scale_r(0.95)
                self.remapMsks()
                self.updatePanels(update_slice=False, update_rotation=True,
                                  update_extent=False)
            elif event.key == 'right':
                self.sectorObj.rotate(-10.0)
                self.remapMsks()
                self.updatePanels(update_slice=False, update_rotation=True,
                                  update_extent=False)
            elif event.key == 'left':
                self.sectorObj.rotate(10.0)
                self.remapMsks()
                self.updatePanels(update_slice=False, update_rotation=True,
                                  update_extent=False)
            else:
                return

    def on_key_release(self, event):
        """Determine what happens if key is released."""
        if event.key == 'control':
            self.ctrlHeld = False

    def findVoxInHist(self, event):
        """Find voxel's location in histogram."""
        self.press = event.xdata, event.ydata
        pixel_x = int(np.floor(event.xdata))
        pixel_y = int(np.floor(event.ydata))
        aoi = self.invHistVolume[:, :, self.sliceNr]  # array of interest
        # Check rotation (TODO: code repetition!)
        cyc_rot = self.cycRotHistory[self.cycleCount][1]
        if cyc_rot == 1:  # 90
            aoi = np.rot90(aoi, axes=(0, 1))
        elif cyc_rot == 2:  # 180
            aoi = aoi[::-1, ::-1]
        elif cyc_rot == 3:  # 270
            aoi = np.rot90(aoi, axes=(1, 0))
        # Switch x and y voxel to get linear index since not Cartesian!!!
        pixelLin = aoi[pixel_y, pixel_x]
        # ind2sub
        xpix = (pixelLin / self.nrBins)
        ypix = (pixelLin % self.nrBins)
        # Switch x and y for circle centre since back to Cartesian.
        circle_colors = [np.array([8, 48, 107, 255])/255,
                         np.array([33, 113, 181, 255])/255]
        self.highlights[0].append(plt.Circle((ypix, xpix), radius=1,
                                  edgecolor=None, color=circle_colors[0]))
        self.highlights[1].append(plt.Circle((ypix, xpix), radius=5,
                                  edgecolor=None, color=circle_colors[1]))
        self.axes.add_artist(self.highlights[0][-1])  # small circle
        self.axes.add_artist(self.highlights[1][-1])  # large circle
        self.figure.canvas.draw()

    def on_press(self, event):
        """Determine what happens if mouse button is clicked."""
        if self.segmType == 'main':
            if event.button == 1:  # left button
                if event.inaxes == self.axes:  # cursor in left plot (hist)
                    if self.ctrlHeld is False:  # ctrl no
                        contains = self.contains(event)
                        if not contains:
                            print 'cursor outside circle mask'
                        if not contains:
                            return
                        # get sector centre x and y positions
                        x0 = self.sectorObj.cx
                        y0 = self.sectorObj.cy
                        # also get cursor x and y position and safe to press
                        self.press = x0, y0, event.xdata, event.ydata
                elif event.inaxes == self.axes2:  # cursor in right plot (brow)
                    self.findVoxInHist(event)
                else:
                    return
            elif event.button == 2:  # scroll button
                if event.inaxes != self.axes:  # outside axes
                    return
                # increase/decrease radius of the sector mask
                if self.ctrlHeld is False:  # ctrl no
                    self.sectorObj.scale_r(1.05)
                    self.remapMsks()
                    self.updatePanels(update_slice=False, update_rotation=True,
                                      update_extent=False)
                elif self.ctrlHeld is True:  # ctrl yes
                    self.sectorObj.rotate(10.0)
                    self.remapMsks()
                    self.updatePanels(update_slice=False, update_rotation=True,
                                      update_extent=False)
                else:
                    return
            elif event.button == 3:  # right button
                if event.inaxes != self.axes:
                    return
                # rotate the sector mask
                if self.ctrlHeld is False:  # ctrl no
                    self.sectorObj.scale_r(0.95)
                    self.remapMsks()
                    self.updatePanels(update_slice=False, update_rotation=True,
                                      update_extent=False)
                elif self.ctrlHeld is True:  # ctrl yes
                    self.sectorObj.rotate(-10.0)
                    self.remapMsks()
                    self.updatePanels(update_slice=False, update_rotation=True,
                                      update_extent=False)
                else:
                    return
        elif self.segmType == 'ncut':
            if event.button == 1:  # left button
                if event.inaxes == self.axes:  # cursor in left plot (hist)
                    xbin = int(np.floor(event.xdata))
                    ybin = int(np.floor(event.ydata))
                    val = self.volHistMask[ybin][xbin]
                    # increment counterField for values in clicked subfield, at
                    # the first click the entire field constitutes the subfield
                    counter = int(self.counterField[ybin][xbin])
                    if counter+1 >= self.ima_ncut_labels.shape[2]:
                        print "already at maximum ncut dimension"
                        return
                    self.counterField[(
                        self.ima_ncut_labels[:, :, counter] ==
                        self.ima_ncut_labels[[ybin], [xbin], counter])] += 1
                    print "counter:" + str(counter+1)
                    # define arrays with old and new labels for later indexing
                    oLabels = self.ima_ncut_labels[:, :, counter]
                    nLabels = self.ima_ncut_labels[:, :, counter+1]
                    # replace old values with new values (in clicked subfield)
                    self.volHistMask[oLabels == val] = np.copy(
                        nLabels[oLabels == val])
                    self.remapMsks()
                    self.updatePanels(update_slice=False, update_rotation=True,
                                      update_extent=False)

                elif event.inaxes == self.axes2:  # cursor in right plot (brow)
                    self.findVoxInHist(event)
                else:
                    return
            elif event.button == 3:  # right button
                if event.inaxes == self.axes:  # cursor in left plot (hist)
                    xbin = int(np.floor(event.xdata))
                    ybin = int(np.floor(event.ydata))
                    val = self.volHistMask[ybin][xbin]
                    # fetch the slider value to get label nr
                    self.volHistMask[self.volHistMask == val] = \
                        np.copy(self.labelNr)
                    self.remapMsks()
                    self.updatePanels(update_slice=False, update_rotation=True,
                                      update_extent=False)

    def on_motion(self, event):
        """Determine what happens if mouse button moves."""
        if self.segmType == 'main':
            # ... button is pressed
            if self.press is None:
                return
            # ... cursor is in left plot
            if event.inaxes != self.axes:
                return
            # get former sector centre x and y positions,
            # cursor x and y positions
            x0, y0, xpress, ypress = self.press
            # calculate difference betw cursor pos on click
            # and new pos dur motion
            # switch x0 & y0 cause volHistMask not Cart
            dy = event.xdata - xpress
            dx = event.ydata - ypress
            # update x and y position of sector,
            # based on past motion of cursor
            self.sectorObj.set_x(x0 + dx)
            self.sectorObj.set_y(y0 + dy)
            # update masks
            self.remapMsks()
            self.updatePanels(update_slice=False, update_rotation=True,
                              update_extent=False)
        else:
            return

    def on_release(self, event):
        """Determine what happens if mouse button is released."""
        self.press = None
        # Remove highlight circle
        if self.highlights[1]:
            self.highlights[1][-1].set_visible(False)
        self.figure.canvas.draw()

    def disconnect(self):
        """Make the object unresponsive."""
        self.figure.canvas.mpl_disconnect(self.cidpress)
        self.figure.canvas.mpl_disconnect(self.cidrelease)
        self.figure.canvas.mpl_disconnect(self.cidmotion)

    def updateColorBar(self, val):
        """Update slider for scaling log colorbar in 2D hist."""
        histVMax = np.power(10, self.sHistC.val)
        plt.clim(vmax=histVMax)

    def updateSliceNr(self):
        """Update slice number and the selected slice."""
        self.sliceNr = int(self.sSliceNr.val*self.orig.shape[2])
        self.imaSlc = self.orig[:, :, self.sliceNr]

    def updateImaBrowser(self, val):
        """Update image browse."""
        # scale slider value [0,1) to dimension index
        self.updateSliceNr()
        self.remapMsks()
        self.updatePanels(update_slice=True, update_rotation=True,
                          update_extent=True)

    def updateImaExtent(self):
        """Update both image and mask extent in image browser."""
        self.imaSlcH.set_extent((0, self.imaSlc.shape[1],
                                 self.imaSlc.shape[0], 0))
        self.imaSlcMskH.set_extent((0, self.imaSlc.shape[1],
                                    self.imaSlc.shape[0], 0))

    def cycleView(self, event):
        """Cycle through views."""
        self.cycleCount = (self.cycleCount + 1) % 3
        # transpose data
        self.orig = np.transpose(self.orig, (2, 0, 1))
        # transpose ima2volHistMap
        self.invHistVolume = np.transpose(self.invHistVolume, (2, 0, 1))
        # updates
        self.updateSliceNr()
        self.remapMsks()
        self.updatePanels(update_slice=True, update_rotation=True,
                          update_extent=True)

    def rotateIma90(self, axes=(0, 1)):
        """Rotate image slice 90 degrees."""
        self.imaSlc = np.rot90(self.imaSlc, axes=axes)
        self.imaSlcMsk = np.rot90(self.imaSlcMsk, axes=axes)

    def changeRotation(self, event):
        """Change rotation of image after clicking the button."""
        self.cycRotHistory[self.cycleCount][1] += 1
        self.cycRotHistory[self.cycleCount][1] %= 4
        self.rotateIma90()
        self.updatePanels(update_slice=True, update_rotation=False,
                          update_extent=True)

    def checkRotation(self):
        """Check rotation update if changed."""
        cyc_rot = self.cycRotHistory[self.cycleCount][1]
        if cyc_rot == 1:  # 90
            self.rotateIma90(axes=(0, 1))
        elif cyc_rot == 2:  # 180
            self.imaSlc = self.imaSlc[::-1, ::-1]
            self.imaSlcMsk = self.imaSlcMsk[::-1, ::-1]
        elif cyc_rot == 3:  # 270
            self.rotateIma90(axes=(1, 0))

    def exportNifti(self, event):
        """Export labels in the image browser as a nifti file."""
        print "start exporting labels..."
        # put the permuted indices back to their original format
        cycBackPerm = (self.cycleCount, (self.cycleCount+1) % 3,
                       (self.cycleCount+2) % 3)
        # assing unique integers (for ncut labels)
        out_volHistMask = np.copy(self.volHistMask)
        labels = np.unique(self.volHistMask)
        intLabels = [i for i in range(labels.size)]
        for label, newLabel in zip(labels, intLabels):
            out_volHistMask[out_volHistMask == label] = intLabels[newLabel]
        # get 3D brain mask
        temp = np.transpose(self.invHistVolume, cycBackPerm)
        outNii = map_2D_hist_to_ima(temp, out_volHistMask)
        outNii = outNii.reshape(temp.shape)
        # save mask image as nii
        new_image = Nifti1Image(outNii, header=self.nii.get_header(),
                                affine=self.nii.get_affine())
        # get new flex file name and check for overwriting
        self.nrExports = 0
        self.flexfilename = '_labels_' + str(self.nrExports) + '.nii.gz'
        while os.path.isfile(self.basename + self.flexfilename):
            self.nrExports += 1
            self.flexfilename = '_labels_' + str(self.nrExports) + '.nii.gz'
        save(new_image, self.basename + self.flexfilename)
        print "successfully exported image labels as: \n" + \
            self.basename + self.flexfilename

    def clearOverlays(self):
        """Clear overlaid items such as circle highlights."""
        if self.highlights[0]:
            {h.remove() for h in self.highlights[0]}
            {h.remove() for h in self.highlights[1]}
        self.highlights[0] = []

    def resetGlobal(self, event):
        """Reset stuff."""
        # reset highlights
        self.clearOverlays()
        # reset color bar
        self.sHistC.reset()
        # reset transparency
        self.TranspVal = 0.5
        if self.segmType == 'main':
            if self.lassoSwitchCount == 1:  # reset only lasso drawing
                self.idxLasso = np.zeros(self.nrBins*self.nrBins, dtype=bool)
            else:
                # reset theta sliders
                self.sThetaMin.reset()
                self.sThetaMax.reset()
                # reset values for mask
                self.sectorObj.set_x(cfg.init_centre[0])
                self.sectorObj.set_y(cfg.init_centre[1])
                self.sectorObj.set_r(cfg.init_radius)
                self.sectorObj.tmin, self.sectorObj.tmax = np.deg2rad(
                    cfg.init_theta)

        elif self.segmType == 'ncut':
            self.sLabelNr.reset()
            # reset ncut labels
            self.ima_ncut_labels = np.copy(self.orig_ncut_labels)
            # reset values for volHistMask
            self.volHistMask = self.ima_ncut_labels[:, :, 0].reshape(
                (self.nrBins, self.nrBins))
            # reset counter field
            self.counterField = np.zeros((self.nrBins, self.nrBins))
            # reset political borders
            self.pltMap = np.zeros((self.nrBins, self.nrBins))
            self.pltMapH.set_data(self.pltMap)
        self.updateSliceNr()
        self.remapMsks()
        self.updatePanels(update_slice=False, update_rotation=True,
                          update_extent=False)

    def updateThetaMin(self, val):
        """Update theta (min) in volume histogram mask."""
        if self.segmType == 'main':
            theta_val = self.sThetaMin.val  # get theta value from slider
            self.sectorObj.theta_min(theta_val)
            self.remapMsks()
            self.updatePanels(update_slice=False, update_rotation=True,
                              update_extent=False)
        else:
            return

    def updateThetaMax(self, val):
        """Update theta(max) in volume histogram mask."""
        if self.segmType == 'main':
            theta_val = self.sThetaMax.val  # get theta value from slider
            self.sectorObj.theta_max(theta_val)
            self.remapMsks()
            self.updatePanels(update_slice=False, update_rotation=True,
                              update_extent=False)
        else:
            return

    def exportNyp(self, event):
        """Export histogram counts as a numpy array."""
        outFileName = self.basename + '_identifier' \
            + '_pcMax' + str(self.initTpl[0]) \
            + '_pcMin' + str(self.initTpl[1]) \
            + '_sc' + str(int(self.initTpl[2]))
        if self.segmType == 'ncut':
            outFileName = outFileName.replace('identifier', 'volHistLabels')
            outFileName = outFileName.replace('.', 'pt')
            np.save(outFileName, self.volHistMask)
            print "successfully exported histogram colors as: \n" + \
                outFileName
        elif self.segmType == 'main':
            outFileName = outFileName.replace('identifier', 'volHist')
            outFileName = outFileName.replace('.', 'pt')
            np.save(outFileName, self.counts)
            print "successfully exported histogram counts as: \n" + \
                outFileName
        else:
            return

    def updateLabels(self, val):
        """Update labels in volume histogram with slider."""
        if self.segmType == 'ncut':
            self.labelNr = self.sLabelNr.val
        else:
            return

    def imaSlcMskIncr(self, incr):
        """Update transparency of image mask by increment."""
        if (self.TranspVal + incr >= 0) & (self.TranspVal + incr <= 1):
            self.TranspVal += incr
        self.imaSlcMskH.set_alpha(self.TranspVal)
        self.figure.canvas.draw()

    def imaSlcMskTransSwitch(self):
        """Update transparency of image mask to toggle transparency."""
        self.imaSlcMskSwitch = (self.imaSlcMskSwitch+1) % 2
        if self.imaSlcMskSwitch == 1:  # set imaSlcMsk transp
            self.imaSlcMskH.set_alpha(0)
        else:  # set imaSlcMsk opaque
            self.imaSlcMskH.set_alpha(self.TranspVal)
        self.figure.canvas.draw()

    def volHistHighlightTransSwitch(self):
        """Update transparency of highlights to toggle transparency."""
        self.volHistHighlightSwitch = (self.volHistHighlightSwitch+1) % 2
        if self.volHistHighlightSwitch == 1 and self.highlights[0]:
            if self.highlights[0]:
                {h.set_visible(False) for h in self.highlights[0]}
        elif self.volHistHighlightSwitch == 0 and self.highlights[0]:
                {h.set_visible(True) for h in self.highlights[0]}
        self.figure.canvas.draw()

    def updateLabelsRadio(self, val):
        """Update labels with radio buttons."""
        labelScale = self.lMax / 6.  # nr of non-zero radio buttons
        self.labelNr = int(float(val) * labelScale)

    def labelContours(self):
        """Plot political borders used in ncut version."""
        grad = np.gradient(self.volHistMask)
        self.pltMap = np.greater(np.sqrt(np.power(grad[0], 2) +
                                         np.power(grad[1], 2)), 0)
        self.pltMapH.set_data(self.pltMap)
        self.pltMapH.set_extent((0, self.nrBins, self.nrBins, 0))

    def lassoArr(self, array, indices):
        """Lasso related."""
        lin = np.arange(array.size)
        newArray = array.flatten()
        newArray[lin[indices]] = 1
        return newArray.reshape(array.shape)

    def calcImaMaskBrd(self):
        """Calculate borders of image mask slice."""
        grad = np.gradient(self.imaSlcMsk)
        return np.greater(np.abs(grad[0], 2) + np.abs(grad[1], 2), 0)


class sector_mask:
    """A pacman-like shape with useful parameters.

    Disclaimer
    ----------
    This script is adapted from a stackoverflow  post by user ali_m:
    [1] http://stackoverflow.com/questions/18352973/mask-a-circular-sector-in-a-numpy-array

    """

    def __init__(self, shape, centre, radius, angle_range):
        self.radius = radius
        self.shape = shape
        self.x, self.y = np.ogrid[:shape[0], :shape[1]]
        self.cx, self.cy = centre
        self.tmin, self.tmax = np.deg2rad(angle_range)
        # ensure stop angle > start angle
        if self.tmax < self.tmin:
            self.tmax += 2*np.pi
        # convert cartesian --> polar coordinates
        self.r2 = (self.x-self.cx)*(self.x-self.cx) + (
            self.y-self.cy)*(self.y-self.cy)
        self.theta = np.arctan2(self.x-self.cx, self.y-self.cy) - self.tmin
        # wrap angles between 0 and 2*pi
        self.theta %= (2*np.pi)

    def set_polCrd(self):
        """Convert cartesian to polar coordinates."""
        self.r2 = (self.x-self.cx)*(self.x-self.cx) + (
            self.y-self.cy)*(self.y-self.cy)
        self.theta = np.arctan2(self.x-self.cx, self.y-self.cy) - self.tmin
        # wrap angles between 0 and 2*pi
        self.theta %= (2*np.pi)

    def set_x(self, x):
        """Set x axis value."""
        self.cx = x
        self.set_polCrd()  # update polar coordinates

    def set_y(self, y):
        """Set y axis value."""
        self.cy = y
        self.set_polCrd()  # update polar coordinates

    def set_r(self, radius):
        """Set radius of the circle."""
        self.radius = radius

    def scale_r(self, scale):
        """Scale (multiply) the radius."""
        self.radius = self.radius * scale

    def rotate(self, degree):
        """Rotate shape."""
        rad = np.deg2rad(degree)
        self.tmin += rad
        self.tmax += rad
        self.set_polCrd()  # update polar coordinates

    def theta_min(self, degree):
        """Angle to determine one the cut out piece in circular mask."""
        rad = np.deg2rad(degree)
        self.tmin = rad
        # ensure stop angle > start angle
        if self.tmax <= self.tmin:
            self.tmax += 2*np.pi
        # ensure stop angle- 2*np.pi NOT > start angle
        if self.tmax - 2*np.pi >= self.tmin:
            self.tmax -= 2*np.pi
        # update polar coordinates
        self.set_polCrd()

    def theta_max(self, degree):
        """Angle to determine one the cut out piece in circular mask."""
        rad = np.deg2rad(degree)
        self.tmax = rad
        # ensure stop angle > start angle
        if self.tmax <= self.tmin:
            self.tmax += 2*np.pi
        # ensure stop angle- 2*np.pi NOT > start angle
        if self.tmax - 2*np.pi >= self.tmin:
            self.tmax -= 2*np.pi
        # update polar coordinates
        self.set_polCrd()

    def binaryMask(self):
        """Return a boolean mask for a circular sector."""
        # circular mask
        self.circmask = self.r2 <= self.radius*self.radius
        # angular mask
        self.anglemask = self.theta <= (self.tmax-self.tmin)
        # return binary mask
        return self.circmask*self.anglemask

    def contains(self, event):
        """Check if a cursor pointer is inside the sector mask."""
        xbin = np.floor(event.xdata)
        ybin = np.floor(event.ydata)
        Mask = self.binaryMask()
        # the next line doesn't follow pep 8 (otherwise it fails)
        if Mask[ybin][xbin] is True:  # switch x and ybin, volHistMask not Cart
            return True
        else:
            return False

    def draw(self, ax, cmap='Reds', alpha=0.2, vmin=0.1,
             interpolation='nearest', origin='lower', extent=[0, 100, 0, 100]):
        """Draw stuff."""
        BinMask = self.binaryMask()
        FigObj = ax.imshow(
            BinMask,
            cmap=cmap,
            alpha=alpha,
            vmin=vmin,
            interpolation=interpolation,
            origin=origin,
            extent=extent)
        return (FigObj, BinMask)
