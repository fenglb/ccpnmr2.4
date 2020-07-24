
"""
======================COPYRIGHT/LICENSE START==========================

NmrPipePseudo.py: Part of the CcpNmr Analysis program

Copyright (C) 2003-2010 Wayne Boucher and Tim Stevens (University of Cambridge)

=======================================================================

The CCPN license can be found in ../../../../license/CCPN.license.

======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)

- email: ccpn@bioc.cam.ac.uk

- contact the authors: wb104@bioc.cam.ac.uk, tjs23@cam.ac.uk
=======================================================================

If you are using this software for academic purposes, we suggest
quoting the following references:

===========================REFERENCE START=============================
R. Fogh, J. Ionides, E. Ulrich, W. Boucher, W. Vranken, J.P. Linge, M.
Habeck, W. Rieping, T.N. Bhat, J. Westbrook, K. Henrick, G. Gilliland,
H. Berman, J. Thornton, M. Nilges, J. Markley and E. Laue (2002). The
CCPN project: An interim report on a data model for the NMR community
(Progress report). Nature Struct. Biol. 9, 416-418.

Wim F. Vranken, Wayne Boucher, Tim J. Stevens, Rasmus
H. Fogh, Anne Pajon, Miguel Llinas, Eldon L. Ulrich, John L. Markley, John
Ionides and Ernest D. Laue (2005). The CCPN Data Model for NMR Spectroscopy:
Development of a Software Pipeline. Proteins 59, 687 - 696.

===========================REFERENCE END===============================

"""

import os
import re

from memops.general.Implementation import ApiError

from memops.gui.Button import Button
from memops.gui.ButtonList import UtilityButtonList
from memops.gui.Entry import Entry
from memops.gui.FileSelectPopup import FileSelectPopup
from memops.gui.FloatEntry import FloatEntry
from memops.gui.Frame import Frame
from memops.gui.IntEntry import IntEntry
from memops.gui.Label import Label
from memops.gui.MessageReporter import showError
from memops.gui.RadioButtons import RadioButtons

from ccpnmr.analysis.popups.BasePopup import BasePopup

from ccp.format.spectra.params.NmrPipeParams import getSampledValues

template_re = re.compile('^(^\D+)(\d+)(\.ft\d)$')

class NmrPipePseudoPopup(BasePopup):

  pseudoEntries = ('Is Pseudo Expt', 'Is Not Pseudo Expt')

  def __init__(self, parent, params, dim, fileName='', *args, **kw):

    self.dim = dim
    self.params = params
    self.fileName = fileName

    m = template_re.match(fileName)
    if m:
      n = len(m.groups(2))
      ss = '%%0%dd' % n
      template = re.sub(template_re, r'\1%s\3' % ss, fileName)
    else:
      template = fileName
    self.template = template

    BasePopup.__init__(self, parent=parent, title='NMRPipe Pseudo Data', modal=True, **kw)

  def body(self, main):

    fileName = self.fileName
    directory = os.path.dirname(fileName)
    if not directory:
      directory = os.getcwd()
    fileName = os.path.basename(fileName)

    m = template_re.match(fileName)
    if m:
      n = len(m.groups(2))
      ss = '%%0%dd' % n
      template = re.sub(template_re, r'\1%s\3' % ss, fileName)
    else:
      template = fileName

    main.rowconfigure(0, weight=1)
    main.rowconfigure(1, weight=1)
    main.columnconfigure(1, weight=1)

    tipTexts = ['The experiment is pseudo-N dimensional, with a sampled axis',
                'The experiment is the regular kind with only NMR frequency axes']
    self.pseudoButton = RadioButtons(main, entries=self.pseudoEntries,
                                     select_callback=self.changedPseudoMode,
                                     grid=(0,0), sticky='nw', tipTexts=tipTexts)

    frame = self.pseudoFrame = Frame(main)
    self.pseudoFrame.grid(row=1, column=0, sticky='nsew')

    row = 0
    npts = self.params.npts[self.dim]
    tipText = 'Number of data points (planes) along sampled axis'
    label = Label(frame, text='Number of points: ')
    label.grid(row=row, column=0, sticky='e')
    self.nptsEntry = IntEntry(frame, text=npts, tipText=tipText, width=8, grid=(row,1))

    tipText = 'Load the values for the sampled axis from a text file containing a list of numeric values'
    Button(frame, text='Load values from text file', command=self.loadValues,
           tipText=tipText, grid=(row,2), sticky='ew')

    row = row + 1
    tipText = 'The values (e.g. T1, T2) corresponding to each data point (plane) along sampled axis'
    label = Label(frame, text='Point values: ')
    label.grid(row=row, column=0, sticky='e')
    self.valueEntry = FloatEntry(frame, isArray=True, tipText=tipText)
    self.valueEntry.grid(row=row, column=1, columnspan=2, sticky='ew')

    row = row + 1
    tipText = 'Fetch the Point values from the files given by the NMRPipe template'
    button = Button(frame, text='Fetch values from file(s) specified by template below', command=self.fetchValues, tipText=tipText)
    button.grid(row=row, column=1, columnspan=2, sticky='w')

    row = row + 1
    tipText = 'The directory where the data files reside'
    button = Button(frame, text='Data directory: ', command=self.chooseDirectory)
    button.grid(row=row, column=0, sticky='e')
    self.directoryEntry = Entry(frame, text=directory, width=40, tipText=tipText)
    self.directoryEntry.grid(row=row, column=1, columnspan=2, sticky='ew')

    row = row + 1
    tipText = 'The NMRPipe template for the data files, if you want to use these to fetch the point values from'
    button = Button(frame, text='NMRPipe template: ', command=self.chooseFile)
    button.grid(row=row, column=0, sticky='e')
    self.templateEntry = Entry(frame, text=template, tipText=tipText)
    self.templateEntry.grid(row=row, column=1, columnspan=2, sticky='ew')

    for n in range(row):
      frame.rowconfigure(n, weight=1)
    frame.columnconfigure(1, weight=1)

    buttons = UtilityButtonList(main, closeText='Ok', doClone=False, 
                                closeCmd=self.updateParams, helpUrl=self.help_url)
    buttons.grid(row=2, column=0, sticky='ew')

  def loadValues(self):

    directory = self.parent.fileSelect.getDirectory()
    fileSelectPopup = FileSelectPopup(self, title='Select Sampled Data File',
                                      dismiss_text='Cancel',
                                      selected_file_must_exist=True,
                                      multiSelect=False,
                                      directory=directory)


    fileName = fileSelectPopup.file_select.getFile()

    if not fileName:
      return

    fileObj = open(fileName, 'rU')

    data = ''
    line = fileObj.readline()
    while line:
      data += line
      line = fileObj.readline()

    fileObj.close()

    data = re.sub(',\s+', ',', data)
    data = re.sub('\s+', ',', data)
    data = re.sub(',,', ',', data)
    data = re.sub('[^0-9,.\-+eE]', '', data)

    self.valueEntry.set(data)

  def chooseDirectory(self):

    directory = os.path.dirname(self.fileName)
    if not directory:
      directory = os.getcwd()
    popup = FileSelectPopup(self, directory=directory, show_file=False)
    directory = popup.getDirectory()
    popup.destroy()

    if directory:
      self.directoryEntry.set(directory)

  def chooseFile(self):

    directory = self.directoryEntry.get()
    if not directory:
      directory = os.getcwd()
    popup = FileSelectPopup(self, directory=directory)
    file = popup.getFile()
    popup.destroy()

    if file:
      template = os.path.basename(file)
      self.templateEntry.set(template)

  def updateParams(self):

    params = self.params
    if self.pseudoButton.get() == self.pseudoEntries[0]:
      npts = self.nptsEntry.get()
      params.npts[self.dim] = npts
      values = self.valueEntry.get()
      try:
        params.setSampledDim(self.dim, values)
      except ApiError, e:
        showError('Set Sampled Dim', e.error_msg, parent=self)
        return
      params.fixNullDims(ignoreDim=self.dim)
    else:
      params.fixNullDims()

    self.close()

  def fetchValues(self):

    directory = self.directoryEntry.get()
    template = self.templateEntry.get()
    try:
      values = getSampledValues(directory, template)
    except ApiError, e:
      showError('Fetch Values', e.error_msg, parent=self)
      return
    self.nptsEntry.set(len(values))
    self.valueEntry.set(values)

  def changedPseudoMode(self, option):

    if option == self.pseudoEntries[0]:
      self.pseudoFrame.grid(row=1, column=0, sticky='nsew')
    else:
      self.pseudoFrame.grid_forget()

