
"""
======================COPYRIGHT/LICENSE START==========================

PartitionedSelector.py: <write function here>

Copyright (C) 2005 Wayne Boucher, Rasmus Fogh, Tim Stevens and Wim Vranken (University of Cambridge and EBI/MSD)

=======================================================================

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
 
A copy of this license can be found in ../../../license/LGPL.license
 
This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
Lesser General Public License for more details.
 
You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA


======================COPYRIGHT/LICENSE END============================

for further information, please contact :

- CCPN website (http://www.ccpn.ac.uk/)
- PDBe website (http://www.ebi.ac.uk/pdbe/)

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
import Tkinter


from memops.gui.Frame  import Frame
from memops.gui.Button import Button
from memops.gui.Color  import invertColor, scaleColor
from memops.gui.Canvas import Canvas

class PartitionedSelector(Frame):

  def __init__(self, parent, callback=None, toggled=True, radio=False,
              objects=None, selected=None, labels=None, colors=None, fonts=None,
              maxRowObjects=18, toggledBg='grey75', toggledFg='grey30',
              font=None, sticky='ew', buttonBorderWidth=1, buttonRelief='raised',
              docKey=None, tipText=None, *args, **kw):
  
    Frame.__init__(self, parent, sticky=sticky, docKey=docKey, tipText=tipText, createToolTip=True, *args, **kw)
    
    # NOTE: toggledBg and toggledFg can be set to None for no change
    # callback should take object 
    
    self.radio     = radio
    self.callback  = callback
    self.toggled   = toggled
    self.toggledFg = toggledFg
    self.toggledBg = toggledBg
    self.state     = []
    self.objects   = objects
    self.labels    = labels
    self.colors    = colors
    self.fonts     = fonts
    self.maxRowObjects     = maxRowObjects
    self.buttonBorderWidth = buttonBorderWidth
    self.buttonRelief      = buttonRelief
    self.font = font or 'Helvetica 10'
    self.waiting = False
    
    self.buttons = []
    
    self.update(objects=objects, labels=labels, selected=selected)
   
    #self.bind('<Configure>',self.refreshSizeAfter)

  def update(self,objects=None,selected=None,labels=None,colors=None,fonts=None):
  
    if objects is None:
      if self.objects is None:
        self.objects = []
    else:
      self.objects = objects
      
    if labels is None:
      if self.labels is None:
        self.labels = []
    else:
      self.labels = labels
      
    if colors is None:
      if self.colors is None:
        self.colors = []
    else:
      self.colors = colors

    if fonts is None:
      if self.fonts is None:
        self.fonts = []
    else:
      self.fonts = fonts
    
    N = len(self.objects)
    for i in range(N):
      if i >= len(self.labels):
        self.labels.append( '%d' % i)
      if self.objects[i] is None:
        self.labels[i] = 'None'
      
      if i>= len(self.colors):
        self.colors.append( self.cget('background'))

      if i>= len(self.fonts):
        self.fonts.append(self.font)
	
      if i>= len(self.buttons):
        command = lambda n=i, s=self: s.selectButton(n)
        button  = Button(self,command=command,relief=self.buttonRelief,
                         font=self.fonts[i],borderwidth=self.buttonBorderWidth)
        self.buttons.append( button )

    if len(self.buttons) > N:
      for i in range(N,len(self.buttons)):
        self.buttons[i].grid_forget()
        self.grid_columnconfigure(i, weight = 0)
     
    row = 0
    col = 0
    for i in range(N):
      if i>= len(self.state):
        self.state.append( 1 )
 
      self.setButtonState(i, self.state[i])
      self.buttons[i].config(text=labels[i],font=self.fonts[i])
      self.buttons[i].grid(row=row,column=col, sticky=Tkinter.NSEW)
      self.grid_columnconfigure(col, weight = 1)
 
      col +=1
      if col >= self.maxRowObjects:
        row += 1
        col  = 0

    if selected is not None:
      self.setSelected(selected)
 
  def refreshSizeAfter(self, event):
    
    if self.waiting:
      return
  
    self.waiting = True
    self.after_idle(lambda: self.refreshSize(event))

  def refreshSize(self, event=None):
    
    canvas   = Canvas(self) 
    textItem = canvas.create_text(0,0,text='A', font=self.font)
    coords   = canvas.bbox(textItem)
    letterWidth = coords[2] - coords[0] + 2
    canvas.destroy()
    
    if not self.buttons:
      return
      
    if event:
      width  = event.width
      height = event.height
    else:
      width  = int(self.winfo_width())
      height = int(self.winfo_height())


    self.maxRowObjects = len(self.buttons)
    x = 0
    for i, button in enumerate(self.buttons):
      text = button.cget('text')
      
      if '\n' in text:
        texts = text.split('\n')
        text = texts[0]
        
        for t in texts[1:]:
          if len(t) > len(text):
            text = t
      
      n  = len(text)
      x1 = x + int(letterWidth*n)
      
      if x1 > width:
        self.maxRowObjects = max(1, i)
        break
    
      x = x1
       
    self.update(self.objects, None,
                self.labels, self.colors,
                self.fonts)  
                
    self.waiting = False

  def setButtonState(self, i, state):
  
    if i < len(self.buttons):
      self.state[i] = not state
      # below will invert the state back
      if self.radio:
        if state:
          self.selectButton(i, doCallback=False)
      else:
        self.selectButton(i, doCallback=False)
	
  def selectButton(self, i, doCallback=True):
  
    if i >= len(self.buttons):
      return
    
    button = self.buttons[i]
    
    if self.radio:
      self.state[i] = 0
   
    if self.toggled:

      if not self.state[i]:
        self.state[i] = 1
        if self.buttonRelief == 'sunken':
          relief = 'raised'
        else:
          relief = 'sunken'

        color = self.colors[i]
        fg = invertColor(self,color) or 'black'
        bg = color or self.cget('background')
        bg2 = scaleColor(self,bg,0.8)
        if bg2 == bg:
          bg2 = scaleColor(self,bg,1.2)
        
        fg2 = invertColor(self,bg2) or 'black'
        button.config(relief=relief, background=bg, fg=fg,
                      activebackground=bg2, activeforeground=fg2)

      elif not self.radio:
        self.state[i] = 0
        if self.buttonRelief == 'sunken':
          relief = 'sunken'
        else:
          relief = 'raised'

        color = self.colors[i]
        fg = self.toggledFg or 'black'
        bg = self.toggledBg or self.cget('background')
        fg2 = invertColor(self,color) or 'black'
        bg2 = color or self.cget('background')
        button.config(relief=relief, background=bg, fg=fg, 
                      activebackground=bg2, activeforeground=fg2)

      if self.radio:
        color = self.colors[i]
        fg = self.toggledFg or 'black'
        bg = self.toggledBg or self.cget('background')
        fg2 = invertColor(self,color) or 'black'
        bg2 = color or self.cget('background')
        self.states = []
        for j in range(len(self.buttons)):
          if j == i:
            self.states.append(1)
            continue
        
          self.states.append(0)
          if self.buttonRelief == 'sunken':
            relief = 'sunken'
          else:
            relief = 'raised'

          self.buttons[j].config(relief=relief, background=bg, fg=fg,
                                 activebackground=bg2, activeforeground=fg2)

    if doCallback and self.callback:
      self.callback(self.objects[i])

  def setSelected(self, objects):
  
    states = []
    for object in self.objects:
      states.append(0)
  
    for object in objects:
      if object in self.objects:
        i = self.objects.index(object)
        states[i] = 1
        
    for i in range(len(states)):
      self.setButtonState(i,states[i])
    
  def getSelected(self):
  
    objects = []
    for i in range(len(self.objects)):
      if self.state[i]:
        objects.append(self.objects[i])
  
    return objects
  

if __name__ == '__main__':

  import Tkinter
  
  root = Tkinter.Tk()
  
  obj = range(10)
  labels = ['%d' % (x+1) for x in obj]
  
  def clicked(*args):
    print "Clicked", args
  
  ps =  PartitionedSelector(root, callback=clicked, toggled=True, radio=True, grid=(0,0),
                            toggledBg='#000000', toggledFg='#80FF80',
                            objects=obj, selected=None, labels=labels, colors=None, fonts=None,
                            font=None, sticky='ew', buttonBorderWidth=1, buttonRelief='raised')
  
  colors = ['#FF0000','#FFFF00','#FF00FF','#00FFFF','#0000FF','#FF0000','#FFFF00','#FF00FF','#00FFFF','#0000FF',]
  ps2 =  PartitionedSelector(root, callback=clicked, toggled=True, radio=False, grid=(1,0),
                            objects=obj, selected=None, labels=labels, colors=colors, fonts=None,
                            font=None, sticky='ew', buttonBorderWidth=1, buttonRelief='sunken')
  
  print ps2.state

  root.mainloop()

