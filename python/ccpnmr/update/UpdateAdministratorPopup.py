"""
======================COPYRIGHT/LICENSE START==========================

UpdateAdministratorPopup.py: Part of the CcpNmr Update program

Copyright (C) 2003-2010 Wayne Boucher and Tim Stevens (University of Cambridge)

=======================================================================

The CCPN license can be found in ../../../license/CCPN.license.

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
import os, sys
import Tkinter

from memops.gui.BasePopup       import BasePopup
from memops.gui.Entry           import Entry
from memops.gui.FileSelect      import FileType
from memops.gui.FileSelectPopup import FileSelectPopup
from memops.gui.IntEntry        import IntEntry
from memops.gui.Label           import Label
from memops.gui.MessageReporter import showWarning
from memops.gui.ScrolledMatrix  import ScrolledMatrix

from memops.gui.ButtonList      import ButtonList

from memops.universal.Io import splitPath

from ccpnmr.update.UpdateAgent import UpdateAgent, FileUpdate, UPDATE_SERVER_LOCATION, UPDATE_DIRECTORY, UPDATE_DATABASE_FILE

class UpdateAdministratorPopup(BasePopup, UpdateAgent):

  def __init__(self, parent, 
               serverLocation=UPDATE_SERVER_LOCATION,
               serverDirectory=UPDATE_DIRECTORY,
               dataFile=UPDATE_DATABASE_FILE):
    
    UpdateAgent.__init__(self, serverLocation, serverDirectory, dataFile, admin=1)

    self.fileTypes = [  FileType('Python', ['*.py']), FileType('C', ['*.c']), FileType('All', ['*'])]
    self.fileUpdate = None

    BasePopup.__init__(self, parent=parent, title='CcpNmr Update Administrator', quitFunc=self.quit)

  def body(self, guiParent):
  
    guiParent.grid_columnconfigure(3, weight=1)
    
    self.commentsEntry = Entry(self)
    self.priorityEntry = IntEntry(self)
    
    
    row = 0
    label = Label(guiParent, text='Server location:')
    label.grid(row=row, column=0, sticky='w')

    location = ''
    uid      = ''
    httpDir   = ''
    subDir   = ''
    version  = 'None'
    
    if self.server:
      location, uid, httpDir, subDir = self.server.identity
      version = self.server.version or 'None'
    
    self.serverEntry = Entry(guiParent, text=location)
    self.serverEntry.grid(row=row, column=1, stick='w')

    label = Label(guiParent, text='User ID:')
    label.grid(row=row, column=2, sticky='w')
    self.uidEntry = Entry(guiParent, text=uid)
    self.uidEntry.grid(row=row, column=3, stick='w')

    row += 1
    label = Label(guiParent, text='HTTP directory:')
    label.grid(row=row, column=0, sticky='w')
    self.httpDirEntry = Entry(guiParent, text=httpDir)
    self.httpDirEntry.grid(row=row, column=1, stick='w')

    label = Label(guiParent, text='Sub-directory:')
    label.grid(row=row, column=2, sticky='w')
    self.subDirEntry = Entry(guiParent, text=subDir)
    self.subDirEntry.grid(row=row, column=3, stick='w')

    row += 1
    self.localVerLabel = Label(guiParent, text='Local version: %s' % self.version)
    self.localVerLabel.grid(row=row, column=0, sticky='w')

    self.serverLabel = Label(guiParent, text='Server version: %s' % version)
    self.serverLabel.grid(row=row, column=2, sticky='w')

    row += 1
    guiParent.grid_rowconfigure(row, weight=1)
    headingList = ['File','Location','Date','Priority','Comments','StoredAs',]

    editWidgets      = [ None, None, None, self.priorityEntry, self.commentsEntry ]
    editGetCallbacks = [ None, None, None, self.getPriority,   self.getComments   ]
    editSetCallbacks = [ None, None, None, self.setPriority,   self.setComments   ]
    self.scrolledMatrix = ScrolledMatrix(guiParent, headingList=headingList,
                                         multiSelect=True,
                                         editWidgets=editWidgets, callback=self.selectCell,
                                         editGetCallbacks=editGetCallbacks,
                                         editSetCallbacks=editSetCallbacks)
    self.scrolledMatrix.grid(row=row, column=0, columnspan=4, sticky='nsew')

    row += 1
    texts = ['Add\nFiles','Remove\nFiles','Remove\nAll','Query\nServer','Commit\nSelected','Synchronise\nAll','Commit\nNew', 'Quit']
    commands = [self.addFile, self.removeFile, self.removeAll, self.queryFiles, self.synchroniseSelected, self.synchroniseServer, self.updateServer, self.quit]
    self.buttonList = ButtonList(guiParent, texts=texts, commands=commands, expands=1)
    self.buttonList.grid(row=row, column=0, columnspan=4, sticky='ew')

    self.update()

  def getPriority(self, fileUpdate):
  
    if fileUpdate:
      self.priorityEntry.set(fileUpdate.priority)
  
  def setPriority(self, event):

    i = self.priorityEntry.get()
    if self.fileUpdate:
      self.fileUpdate.priority = i
    
    self.update()

  def getComments(self, fileUpdate):

    if fileUpdate:
      self.commentsEntry.set(fileUpdate.details)
  
  def setComments(self, event):

    text = self.commentsEntry.get()
    if self.fileUpdate:
      self.fileUpdate.details = text
    
    self.update()

  def quit(self):
  
    self.close()
    self.destroy()
    sys.exit()


  def updateServer(self):
  
    if self.server:
      # Unly upload new updates
      success = self.server.setFileUpdates()
      self.serverLabel.set('Server version: %s' % self.server.version)
      if success:
        self.queryFiles()

  def synchroniseSelected(self):
  
    if self.server:
      selectedUpdates = self.scrolledMatrix.currentObjects
      if selectedUpdates:
        # Refresh all selected updates
        success = self.server.setFileUpdates(fileUpdates=selectedUpdates, refresh=True)
        self.serverLabel.set('Server version: %s' % self.server.version)
        if success:
          self.queryFiles()

  def synchroniseServer(self):
  
    if self.server:
      # Refresh all avaiable updates
      success = self.server.setFileUpdates(refresh=True)
      self.serverLabel.set('Server version: %s' % self.server.version)
      if success:
        self.queryFiles()

  def queryFiles(self):
  
    if self.server:
      self.server.getFileUpdates()
        
    self.update()

  def addFile(self):
  
    if self.server:
      fileSelectPopup = FileSelectPopup(self, title='Select Source File', dismiss_text='Cancel', 
                                        file_types=self.fileTypes, multiSelect=True, 
                                        selected_file_must_exist=True)
      
      filePaths = fileSelectPopup.file_select.getFiles()
      n = len(self.installRoot)

      for filePath in filePaths:
        if self.installRoot != filePath[:n]:
          showWarning('Warning','Install root %s not found in file path %s' % (self.installRoot, filePath))
          continue
 
        filePath = filePath[n+1:]
 
        if filePath:
          dirName, fileName = splitPath(filePath)
 
          if fileName[-3:] == '.py':
            language = 'python'
          else:
            language = 'None'
 
          fileUpdate = FileUpdate(self.server, fileName, dirName, language, isNew=True)
 
      self.update()
     

  def removeFile(self):
  
    if self.fileUpdate:
      self.fileUpdate.delete()
     
      self.update() 

  def removeAll(self):
  
    if self.server:
      for fileUpdate in list(self.server.fileUpdates):
        fileUpdate.delete()

      self.update() 
 
  def selectAll(self):
  
    if self.server:
      for fileUpdate in self.server.fileUpdates:
        fileUpdate.isSelected = 1
        
    self.update()

  def selectCell(self, object, row, col):
  
    self.fileUpdate = object
    self.updateButtons()

  def updateButtons(self):
  
    buttons = self.buttonList.buttons
    
    if self.server:
      buttons[0].enable()
      buttons[3].enable()
      buttons[4].enable()
    else:
      buttons[0].disable()  
      buttons[3].disable()  
      buttons[4].disable()  
    
    if self.server and self.server.fileUpdates:
      buttons[2].enable()
    else:
      buttons[2].disable()  

    if self.fileUpdate:
      buttons[1].enable()
    else:
      buttons[1].disable()  
      
  def update(self):
    
    location = self.serverEntry.get()
    uid      = self.uidEntry.get()
    httpDir   = self.httpDirEntry.get()
    subDir   = self.subDirEntry.get()
    if self.server:
      if (location,uid,httpDir,subDir) != self.server.identity:
        self.setServer(location)
    
    self.updateButtons()
    
    self.fileUpdate = None
    textMatrix  = []
    objectList  = []
    colorMatrix = []
    
    if self.server:
      
      for fileUpdate in self.server.fileUpdates:
                
        datum = []
        datum.append(fileUpdate.fileName)
        datum.append(fileUpdate.filePath)
        datum.append(fileUpdate.date)
        datum.append(fileUpdate.priority)
        datum.append(fileUpdate.details)
        datum.append(fileUpdate.storedAs)
        
        textMatrix.append(datum)
        objectList.append(fileUpdate)
        if fileUpdate.isNew:
          colorMatrix.append(5 * ['#B0FFB0'])
        
        elif not fileUpdate.getIsUpToDate():
          colorMatrix.append(5 * ['#FFB0B0'])
        
        else:
          colorMatrix.append(5 * [None])
    
    
    self.scrolledMatrix.update(textMatrix=textMatrix, objectList=objectList, colorMatrix=colorMatrix)
    

if __name__ == '__main__':

  root = Tkinter.Tk()
  root.withdraw()
  # server location, temp dir, dataFile
  top = UpdateAdministratorPopup(root)
  root.mainloop()
