#!/usr/bin/env python

"""
=========================================================================
Package:    - Code for the Graphical User Interface fontend to the 
              Haddock model package in the CCPN data model. 
            - Code for the export of a Haddock compatible project. A 
              Haddock compatible project can either be a parameter file
              ready for submission to the Haddock webserver or a
              directory structure with necessary files for use with a 
              localy installed version of Haddock.

Dependencies: The CCPN Haddock package requires CCPN data model version
              2.0 or higher. The export of a webserver compatible 
              parameter file requires Haddock webserver version 2.1 or 
              higher and a valid user account. The export of a 'classic' 
              Haddock project requires Haddock version 2.0 or higher.

Copyright and License information:
              The Haddock data model as implemented in the CCPN data
              model as well as the use of CCPN GUI code elements is 
              licenced to the CCPN Projects (Copyright (C) 2008) and
              distributed under the terms of the GNU Lesser General
              Public License.
            
              The Haddock project export code as well as the use of 
              Haddock software is covert in the Haddock License
              agreement (Copyright (C) 2008 Haddock Project, Bijvoet
              Center for Biomolecular Research, Utrecht University,
              The Netherlands).

GNU LGPL:        This library is free software; you can redistribute it 
              and/or modify it under the terms of the GNU Lesser General 
              Public License as published by the Free Software 
              Foundation; either version 2.1 of the License, or (at 
              your option) any later version.
 
              A copy of this license can be found in LGPL.license
 
              This library is distributed in the hope that it will be 
              useful, but WITHOUT ANY WARRANTY; without even the implied 
              warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
              PURPOSE. See the GNU Lesser General Public License for 
              more details.
 
              You should have received a copy of the GNU Lesser General 
              Public License along with this library; if not, write to 
              the Free Software Foundation, Inc., 59 Temple Place, Suite 
              330, Boston, MA 02111-1307 USA.

Information:  For further information regarding CCPN, please contact:
              - CCPN website (http://www.ccpn.ac.uk/)
              - email: ccpn@bioc.cam.ac.uk
              
              For further information regarding Haddock, please contact
              Alexandre M.J.J. Bonvin:
              - http://haddock.chem.uu.nl
              - email: a.m.j.j.bonvin@uu.nl    

Citing:          If you are using this software for academic purposes, we 
                suggest quoting the following references:

              For CCPN:    
              Rasmus H. Fogh, Wayne Boucher, Wim F. Vranken, Anne
              Pajon, Tim J. Stevens, T.N. Bhat, John Westbrook, John M.C. 
              Ionides and Ernest D. Laue (2005). A framework for 
              scientific data modeling and automated software development. 
              Bioinformatics 21, 1678-1684.
            
              For Haddock:
              Cyril Dominguez, Rolf Boelens and Alexandre M.J.J. Bonvin 
              (2003). HADDOCK: a protein-protein docking approach based 
              on biochemical and/or biophysical information. 
              J. Am. Chem. Soc. 125, 1731-1737.
            
              S.J. de Vries, A.D.J. van Dijk, M. Krzeminski, M. van Dijk, 
              A. Thureau, V. Hsu, T. Wassenaar and A.M.J.J. Bonvin (2007) 
              HADDOCK versus HADDOCK: New features and performance of 
              HADDOCK2.0 on the CAPRI targets. 
              Proteins: Struc. Funct. & Bioinformatic 69, 726-733.    
=========================================================================
"""

from memops.gui.ButtonList      import ButtonList
from memops.gui.CheckButton     import CheckButton
from memops.gui.IntEntry        import IntEntry
from memops.gui.LabelFrame      import LabelFrame
from memops.gui.Label           import Label
from memops.gui.MultiWidget     import MultiWidget
from memops.gui.PulldownMenu    import PulldownMenu
from memops.gui.ScrolledMatrix  import ScrolledMatrix
from memops.gui.MessageReporter import showWarning

from memops.editor.BasePopup    import BasePopup
from memops.editor.Util         import createDismissHelpButtonList

from HaddockBasic          import setPartnerChains

class EditSymmetryPopup(BasePopup):

    def __init__(self, parent, project):

        self.parent         = parent
        self.ccpnProject    = project
        self.hProject        = project.currentHaddockProject
        self.molPartner     = None
        self.molecules      = []
        self.symmetrySet    = None
        self.symmetryOp     = None
        self.symmetryCode    = None
        self.waiting        = False

        BasePopup.__init__(self, parent=parent, title='Symmetry Operations')

        self.font = 'Helvetica 12'
        self.setFont()

    def body(self, guiFrame):

        guiFrame.grid_columnconfigure(0, weight=1)
        guiFrame.grid_rowconfigure(1, weight=1)

        frame = LabelFrame(guiFrame, text='Options')
        frame.grid(row=0,column=0,sticky='ew')
        frame.grid_columnconfigure(5, weight=1)

        # Choose type of symmetry set to define or change (if allready present in the model)
        self.molLabel = Label(frame, text='Symmetry Operator:')
        self.molLabel.grid(row=0,column=2,sticky='w')
        self.symmCodePulldown = PulldownMenu(frame, callback=self.setSymmCode, entries=['NCS','C2','C3','C5'], do_initial_callback=False)
        self.symmCodePulldown.grid(row=0,column=3,sticky='w')

        frame = LabelFrame(guiFrame, text='Symmetry Operations')
        frame.grid(row=1,column=0,sticky='nsew')
        frame.grid_columnconfigure(0, weight=1)
        frame.grid_rowconfigure(0, weight=1)

        self.molSysPulldown   = PulldownMenu(self, callback=self.setMolSystem, do_initial_callback=False)
        self.chainSelect      = MultiWidget(self, CheckButton, callback=self.setChains, minRows=0, useImages=False)
        self.segStartEntry    = IntEntry(self, returnCallback=self.setSegStart, width=6)
        self.segLengthEntry   = IntEntry(self, returnCallback=self.setSegLength, width=6)

        headings = ['#','Symmetry\noperator','Mol System','Chains','Start\nresidue','Segment\nlength']  
        editWidgets      = [None, None, self.molSysPulldown, self.chainSelect, self.segStartEntry, self.segLengthEntry] 
        editGetCallbacks = [None, None, self.getMolSystem, self.getChains, self.getSegStart, self.getSegLength]
        editSetCallbacks = [None, self.setSymmCode, self.setMolSystem, self.setChains, self.setSegStart, self.setSegLength]
        
        self.symmetryMatrix = ScrolledMatrix(frame,headingList=headings,
                                             callback=self.selectSymmetry,
                                             editWidgets=editWidgets,
                                             editGetCallbacks=editGetCallbacks,
                                             editSetCallbacks=editSetCallbacks)
        self.symmetryMatrix.grid(row=0,column=0,sticky='nsew')

        texts = ['Add Symmetry Set','Remove Symmetry Set']
        commands = [self.addSymmetrySet,self.removeSymmetrySet]
        self.buttonList = createDismissHelpButtonList(guiFrame, texts=texts, commands=commands, expands=True)
        self.buttonList.grid(row=2,column=0,sticky='ew')
        
        self.updateMolPartners()
        self.notify(self.registerNotify)

        #Temporary report of parameters
        print self.molSystem
        print self.molecules
        print self.symmetrySet
        print self.symmetryOp
        print self.symmetryCode

    def getMolSystem(self, partner):

        """Select molecular system from list of molsystems stored in the project"""

        names = []; index = -1
        molSystem = partner.molSystem
        molSystems = self.ccpnProject.sortedMolSystems()
        if molSystems:
            names = [ms.code for ms in molSystems]
            if molSystem not in molSystems: 
                molSystem = molSystems[0]
                index = molSystems.index(molSystem)

        self.molSysPulldown.setup(names, index)

    def getChains(self, partner):
  
        names  = []
        values = []

        molSystem = partner.molSystem
        if molSystem:
            for chain in molSystem.sortedChains():
                names.append(chain.code)

                if not partner.chains:
                  values.append(True)

                elif partner.findFirstChain(chain=chain):
                  values.append(True)

                else:
                  values.append(False)   

                self.chainSelect.set(values=values,options=names)   

        else:
            showWarning('Warning','Set Mol System or ensemble first',parent=self)
            self.symmetryMatrix.keyPressEscape()        
    
    def getSegLength(self, symmetryOp):
        
        if symmetryOp and symmetryOp.segmentLength: self.segLengthEntry.set(symmetryOp.segmentLength)

    def getSegStart(self):
    
        pass

    def setMolSystem(self, partner, name=None):
  
        """Get all molsystems as stored in the project as list"""

        index =  self.molSysPulldown.getSelectedIndex()
        molSystems = self.ccpnProject.sortedMolSystems()
  
        if molSystems:
            molSystem = molSystems[index]
            self.molPartner.molSystem = molSystem
            chains = molSystem.sortedChains()
            if (not self.molPartner.chains) and chains: setPartnerChains(self.molPartner,chains)

        self.updateAllAfter()

    def setChains(self, obj):

        """Get the list of chains for the selected molsystem"""

        if self.molPartner and self.molPartner.molSystem:
            if obj is not None:
                chains = self.molPartner.molSystem.sortedChains()
                values = self.chainSelect.get()
                chains = [chains[i] for i in range(len(values)) if values[i]]
                setPartnerChains(self.molPartner,chains)

        self.symmetryMatrix.keyPressEscape()  
        self.updateAllAfter()
        
    def setSegLength(self, event):
        
        value = self.segLengthEntry.get() or 1
        self.symmetryOp.segmentLength = value

    def setSegStart(self):
        
        pass    

    def setSymmCode(self, index, name=None): 

        self.symmetryCode = self.symmCodePulldown.getSelected()    

    def selectSymmetry(self, obj, row, col):
        
        self.symmetryOp = obj
        self.updateMolSysPulldown()

    def notify(self, notifyFunc):
  
        for func in ('__init__', 'delete', 'setSymmetryCode','setSegmentLength'):
            notifyFunc(self.updateAllAfter, 'molsim.Symmetry.Symmetry', func)

        for func in ('__init__', 'delete','setfirstSeqId'):
            notifyFunc(self.updateAllAfter, 'molsim.Symmetry.Segment', func)

    def addSymmetrySet(self):
            
        if not self.ccpnProject.molSystems:
            showWarning('Warning','No molecular systems present in CCPN project',parent=self)
            return

            molSystem = self.ccpnProject.findFirstMolSystem()
            partner = self.hProject.newHaddockPartner(code=molSystem.code, molSystem=molSystem)
            setPartnerChains(partner, molSystem.chains)

            self.updateAllAfter()
            
    def removeSymmetrySet(self):
        
        pass            
    
    def updateAllAfter(self, obj=None):

        if self.waiting: return
        else:
            self.waiting = True
            self.after_idle(self.updateSymmetries,
                            self.updateMolPartners)
    
    def updateMolPartners(self):

        textMatrix = []
        objectList = []

        for partner in self.hProject.sortedHaddockPartners():
            datum = [partner.code,
                     self.symmetryCode,
                     partner.molSystem.code,
                     ','.join([c.chain.code for c in partner.chains]),
                     partner.autoHistidinePstate and 'Yes' or 'No',
                     partner.isDna and 'Yes' or 'No']
            
            objectList.append(partner)
            textMatrix.append(datum)

        self.symmetryMatrix.update(objectList=objectList, textMatrix=textMatrix)
        self.updateMolSysPulldown()
    
    def updateMolSysPulldown(self):

        names = []
        index = -1

        partners = self.hProject.sortedHaddockPartners()

        if partners:
            if self.molPartner not in partners:
                self.molPartner = partners[0]

            names = ['Partner %s' % p.code for p in partners]
            index = partners.index(self.molPartner)

        else: self.molPartner = None

        self.molSysPulldown.setup(names, index)
        
    def updateSymmetries(self):

        textMatrix  = []; objectList  = []
        if self.symmetrySet:
            for symmetryOp in self.symmetrySet.symmetries:
                chains = []; segments = [] 
                length = symmetryOp.segmentLength

                for segment in symmetryOp.sortedSegments():
                    code = segment.chainCode
                    chain = self.molSystem.findFirstChain(code=code)

                    if chain:
                        chains.append(code)
                        seqId = segment.firstSeqId
                        residue1 = chain.findFirstResidue(seqId=seqId)
                        residue2 = chain.findFirstResidue(seqId=seqId+length-1)
                        segments.append('%s:%d-%d' % (code,residue1.seqCode,residue2.seqCode))

                datum = [symmetryOp.serial, symmetryOp.symmetryCode, length, '\n'.join(chains), '\n'.join(segments)]
                objectList.append(symmetryOp)
                textMatrix.append(datum)

        self.symmetryMatrix.update(objectList=objectList, textMatrix=textMatrix)
        self.waiting = False

    def destroy(self):

        self.notify(self.unregisterNotify)
        BasePopup.destroy(self)



"""
  def getSegLength(self, symmetryOp):
  
    if symmetryOp and symmetryOp.segmentLength: self.segLengthEntry.set(symmetryOp.segmentLength)
  
  def getChains(self, symmetryOp): 
  
    chains = []
    for chain in self.molSystem.chains:
        if chain.residues:
            if chain.molecule in self.molecules: chains.append(chain.code)
    chains.sort()

    values = []
    for chain in chains:
        if symmetryOp.findFirstSegment(chainCode=chain): values.append(True)
        else: values.append(False)

    self.setChainMulti.set(values=values,options=chains)

  def getSegments(self, symmetryOp):
 
    values = []; names  = []

    if symmetryOp:
        for segment in symmetryOp.sortedSegments():
            names.append(segment.chainCode)
            values.append(segment.firstSeqId)

    n = len(values)
    self.setSegmentMulti.maxRows = n
    self.setSegmentMulti.minRows = n
    self.setSegmentMulti.set(values=values,options=names)

  def setSegLength(self, event):
  
    value = self.segLengthEntry.get() or 1
    self.symmetryOp.segmentLength = value

  def setChains(self, obj):
  
    if self.symmetryOp and obj:
        codes   = self.setChainMulti.options
        segment = self.symmetryOp.findFirstSegment()
        values  = self.setChainMulti.get()

        if segment: seqId0 = segment.firstSeqId
        else: seqId0 = 1

        for i in range(len(values)):
            segment = self.symmetryOp.findFirstSegment(chainCode=codes[i])

            if segment and not values[i]: segment.delete()
            elif values[i] and not segment:
                chain = self.molSystem.findFirstChain(code=codes[i])
                residue = chain.findFirstResidue(seqid=seqId0)

                if residue: seqId = seqId0
                else:
                    residue = chain.sortedResidues()[0]
                    seqId = residue.seqId
                    
                residue2 = chain.findFirstResidue(seqid=seqId+self.symmetryOp.segmentLength)
                if not residue2:
                    residue2 = chain.sortedResidues()[-1]
                    self.symmetryOp.segmentLength = (residue2.seqId - seqId) + 1

                segment = self.symmetryOp.newSegment(chainCode=codes[i],firstSeqId=seqId)
                
    self.symmetryMatrix.keyPressEscape()
    
  def setSegments(self, obj):
  
    if self.symmetryOp and obj:
        segments = self.symmetryOp.sortedSegments()
        values   = self.setSegmentMulti.get()

        for i in range(len(values)):
            seqCode = values[i]
            chain = self.molSystem.findFirstChain(code=segments[i].chainCode)
            residue = chain.findFirstResidue(seqCode=seqCode)
  
            if residue: seqId = residue.seqId
            if segments[i].firstSeqId != seqId:
                segments[i].delete()
                segments[i] = self.symmetryOp.newSegment(chainCode=chain.code,firstSeqId=seqId)

    self.symmetryMatrix.keyPressEscape()

  def addSymmOp(self):
  
    if self.molSystem:
        if not self.symmetrySet: self.symmetrySet = self.molSystem.findFirstMolSystemSymmetrySet()

        if not self.symmetrySet:
            objGen = self.project.newMolSystemSymmetrySet
            self.symmetrySet = objGen(symmetrySetId=1,molSystem=self.molSystem)

        segLen   = len(self.molSystem.findFirstChain().residues)
        symmetry = self.symmetrySet.newSymmetry(segmentLength=segLen)

  def removeSymmOp(self):
  
    if self.symmetryOp: self.symmetryOp.delete()

  def toggleSingleMolecule(self, boolean):
  
    self.singleMolecule = not boolean
    self.updateMolSystems()
    self.updateMolecules()

  def setMolecules(self, molecules):

    self.molecules = molecules
    if self.symmetrySet:
        for symmetryOp in self.symmetrySet.symmetries:
            for segment in symmetryOp.segments: chain = self.molSystem.findFirstChain(code=segment.chainCode)
              if chain and (chain.molecule not in molecules): segment.delete()

  def selectMolecule(self, index, name):
  
    self.setMolecules(self.getMolecules()[index])
    self.updateSymmetries()

  def getMolecules(self):
  
    moleculesList = []

    for chain in self.molSystem.chains: moleculesList.append(chain.molecule)

    return moleculesList

  def updateMolecules(self):
    
    names = []; index = -1
    moleculesList = self.getMolecules()
    if moleculesList:
        if self.molecules not in moleculesList:
            self.setMolecules(moleculesList[0])
            self.symmetrySet = self.molSystem.findFirstMolSystemSymmetrySet()
            
        index = moleculesList.index(self.molecules)  

        names = []
        for molecules in moleculesList: names.append(','.join([mol.name for mol in molecules]))
  
    else: self.molecules = []

    self.moleculePulldown.setup(names, index)

  def selectMolSystem(self, index, name):
  
    self.molSystem = self.getMolSystems()[index]
    self.symmetrySet = self.molSystem.findFirstMolSystemSymmetrySet()
    self.updateSymmetries()

  def getMolSystems(self):

    molSystems = []
    for molSystem in self.project.sortedMolSystems():
        n = len(molSystem.chains)

        if self.singleMolecule and (n > 1): molSystems.append(molSystem)
        elif n > 0: molSystems.append(molSystem)
        
    return molSystems
            
  def updateMolSystems(self):
  
    names = []; index = -1
    molSystems = self.getMolSystems()

    if molSystems:
        if self.molSystem not in molSystems: self.molSystem = molSystems[0]
    
        index = molSystems.index(self.molSystem)
        names = [ms.code for ms in molSystems]  
    else: self.molSystem = None

    self.molSystemPulldown.setup(names, index)

"""
