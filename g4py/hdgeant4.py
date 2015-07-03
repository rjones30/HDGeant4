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
#   >>> hdgeant4.init()
#
# author: richard.t.jones at uconn.edu
# version: june 29, 2015
#

from Geant4 import *
from HDGeant4 import *

def init():
  # initialize the jana framework
  global dapp
  dapp = DApplication()
  dapp.Init()
  global opts
  opts = GlueXUserOptions()
  opts.ReadControl_in("control.in")

  # define the detector geometry
  global geom
  geom = GlueXDetectorConstruction()
  for para in range(1, geom.GetParallelWorldCount() + 1):
    name = geom.GetParallelWorldName(para)
    topvol = geom.GetParallelWorldVolume(para)
    pworld = GlueXParallelWorld(name, topvol)
    geom.RegisterParallelWorld(pworld);
  gRunManager.SetUserInitialization(geom)

  # initialize physics processes
  global plist
  plist = GlueXPhysicsList()
  gRunManager.SetUserInitialization(plist)
   
  # initialize event generators
  global gen
  gen = GlueXPrimaryGeneratorAction(geom)
  gRunManager.SetUserAction(gen);

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
      
  # initialize graphics
  gVisManager.Initialize()

  return 0
