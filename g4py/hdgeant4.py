#
# hdgeant4.py : python interface to the hdgeant4 physics simulation context
#
# purpose:
#   This python executable provides a means for users to instantiate the
#   hdgeant4 simulation context (geometry, magnetic fields, event generators,
#   physics processes, tracking controls, hits collection and output classes,
#   etc.) within a flexible interactive framework. All of the primary geant4
#   classes and many of the user classes are exposed as python objects to 
#   enable runtime customization of many aspects of the simulation. This is
#   not intended to replace the standard C++ executable for a production tool
#   although it could be used that way, but rather as a tool for examining
#   the behavior of the simulation interactively without having to build a
#   custom executable to answer every question about the behavior of the
#   simulation.
#
# usage example:
#   $ python
#   >>> import hdgeant4
#   >>> hdgeant4.pick_point3D(0,0,65)
#   >>> hdgeant4.stdviews()
#   >>> hdgeant4.gUImanager.ApplyCommand("/vis/viewer/zoomTo 50")
#   >>> hdgeant4.gUImanager.ApplyCommand("/vis/ogl/export pic.eps")
#
# author: richard.t.jones at uconn.edu
# version: june 29, 2015
#

import sys

from G4fixes import *
from Geant4 import *
from HDGeant4 import *

def init():
  '''
  Initialize the jana framework and the Geant4 toolkit runtime,
  and load the geometry and magnetic field definitions for GlueX.
  '''
  # initialize the jana framework
  global dapp
  dapp = DApplication()
  dapp.Init()
  global opts
  opts = GlueXUserOptions()
  if not opts.ReadControl_in("control.in"):
    raise IOError("simulation needs control.in, init failed!")
    sys.exit(3)

  # define the detector geometry
  global geom
  geom = GlueXDetectorConstruction()
  for para in range(1, geom.GetParallelWorldCount() + 1):
    name = geom.GetParallelWorldName(para)
    topvol = geom.GetParallelWorldVolume(para)
    pworld = GlueXParallelWorld(name, topvol)
    geom.RegisterParallelWorld(pworld)
  gRunManager.SetUserInitialization(geom)

  # initialize physics processes
  global plist
  plist = GlueXPhysicsList(geom)
  gRunManager.SetUserInitialization(plist)

  # initialize event generators
  global gen
  gen = GlueXPrimaryGeneratorAction()
  gRunManager.SetUserAction(gen)

  # initialize run/event/step actions
  global runact
  global eventact
  global stepact
  runact = GlueXRunAction()
  eventact = GlueXEventAction()
  stepact = GlueXSteppingAction()
  gRunManager.SetUserAction(runact)
  gRunManager.SetUserAction(eventact)
  gRunManager.SetUserAction(stepact)

  # initialize G4 kernel
  gRunManager.Initialize()
      
def stdviews():
  '''
  Generate standard views of the simulation geometry and save
  them as graphics files.
  '''
  ui = gUImanager.GetUIpointer()
  ui.ApplyCommand("/control/execute vis1.mac")
  for view in "z-2340", "z-2300", "z-1750", "z-950", "z0", \
              "z15", "z20", "z65", "z100", "z160", "z170", \
              "z200", "z380", "z412", "z630", "z700":
    ui.ApplyCommand("/control/execute ../vis/stdviews/" + view + ".mac")
    ui.ApplyCommand("/vis/ogl/export " + view + ".eps")

def checkBfield(infile):
  '''
  Checks the magnetic field at a user-selected set of points 
  listed as 3D coordinates (cm) in a 3-column input text file.
  '''
  G4GeometryManager.GetInstance().CloseGeometry()
  for line in open(infile):
    coord = line.split()
    pickPoint3D(float(coord[0]), float(coord[1]), float(coord[2]))

# Automatic actions at startup
init()
