#
# cobrems.py : python class for computing coherent bremsstrahlung spectra
#              and rates, and for generating Monte Carlo samples of them.
#
# usage example:
#   $ python
#   >>> from ROOT import *
#   >>> import cobrems
#   >>> cobrems.init(E0=12, Epeak=8.7)
#   >>> cobrems.plotTotal(0)
#
# author: richard.t.jones at uconn.edu
# version: july 29, 2015
#

import sys
from array import array

from ROOT import *
from Cobrems import *

# set global variables to defaults
Epeak = 9
cur = 2.2
E0 = 12
Erms = 0.0006
emit = 2.5e-9
dist = 76
coldiam = 0.0034
radt = 20e-6
nbins = 200
Emin = 0
Emax = E0
pElim0 = 8.4
pElim1 = 9.0
bElim0 = 1.02
bElim1 = 2.04
eElim0 = 10.68
eElim1 = 11.7

# initialize the built-in generator to defaults
generator = CobremsGenerator(E0, Epeak)

def usage():
   print """
   Usage: cobrems.init(key=value, ...)
     where key should be one of the following, [defaults in backets]:
      * Epeak - photon energy (GeV) of primary coherent edge [9]
      * cur - current (uA) of electron beam [2.2]
      * E0 - mean energy (GeV) of electron beam [12]
      * Erms - rms energy width (GeV) of electron beam [0.0006]
      * emit - transverse emittance (m rad) of electron beam [2.5e-9]
      * dist - distance (m) from radiator to primary collimator [76]
      * coldiam - diameter (m) of primary collimator [0.0034]
      * radt - thickness (m) of radiator crystal [20e-6]
      * nbins - number of bins to use for plotting spectra [200]
      * Emin - low energy limit (GeV) for spectrum historgrams [0]
      * Emax - high energy limit (GeV) for spectrum historgrams [12]
      * pElim0 - lower limit (GeV) of primary peak integration window [8.4]
      * pElim1 - upper limit (GeV) of primary peak integration window [9.0]
      * bElim0 - lower limit (GeV) of beam background integration window [1.02]
      * bElim1 - upper limit (GeV) of beam background integration window [2.04]
      * eElim0 - lower limit (GeV) of beam endpoint integration window [10.68]
      * eElim1 - upper limit (GeV) of beam endpoint integration window [11.70]
   """

def init(**kwargs):
   """
   Loads a new set of beamline parameters and generates the acceptance
   function and collimated beam intensity spectra as ROOT TF1/TH1 objects.
   To see a list of supported options as keyword arguments, call init as
    >>> cobrems.init(o='help')
   """
   for arg, value in kwargs.iteritems():
      if arg == "E0":
         global E0
         E0 = value
         generator.setBeamEnergy(E0)
      elif arg == "Epeak":
         global Epeak
         Epeak = value
         generator.setCoherentEdge(Epeak)
      elif arg == "cur":
         global cur
         cur = value
      elif arg == "Erms":
         global Erms
         Erms = value
         generator.setBeamErms(Erms)
      elif arg == "emit":
         global emit
         emit = value
         generator.setBeamEmittance(emit)
      elif arg == "dist":
         global dist
         dist = value
         generator.setCollimatorDistance(dist)
      elif arg == "coldiam":
         global coldiam
         coldiam = value
         generator.setCollimatorDiameter(coldiam)
      elif arg == "radt":
         global radt
         radt = value
         generator.setTargetThickness(radt)
      elif arg == "nbins":
         global nbins
         nbins = value
      elif arg == "Emin":
         global Emin
         Emin = value
      elif arg == "Emax":
         global Emax
         Emax = value
      elif arg == "pElim0":
         global pElim0
         pElim0 = value
      elif arg == "pElim1":
         global pElim1
         pElim1 = value
      elif arg == "bElim0":
         global bElim0
         bElim0 = value
      elif arg == "bElim1":
         global bElim1
         bElim1 = value
      elif arg == "eElim0":
         global eElim0
         eElim0 = value
      elif arg == "eElim1":
         global eElim1
         eElim1 = value
      else:
         usage()
         return

   generator.printBeamlineInfo()
   generator.setCollimatedFlag(true)
   generator.setPolarizedFlag(false)

   # plot the acceptance curve
   global acceptF1
   acceptF1 = TF1("acceptF1", acceptance, 0, 2, 0)
   acceptF1.SetTitle("\\mbox{acceptance vs }\\theta/\\theta_C")

   # plot the beam energy spectrum
   x0 = Emin / E0
   x1 = Emax / E0
   global dRtdxF1
   dRtdxF1 = TF1("dRtdxF1", dRtdx, x0, x1, 0)
   dRtdxF1.SetTitle("\\mbox{photon beam spectrum vs }E_\\gamma\\,/E_0" +
                    "\\mbox{ (/GeV/s)}")
   dRtdxF1.SetNpx(nbins)

   # apply the beam-crystal convolution
   xvals = array('d', nbins * [0])
   yvals = array('d', nbins * [0])
   for i in range(0, nbins):
      xvals[i] = x0 + (i + 0.5) * (x1 - x0) / nbins;
      yvals[i] = dRtdx([xvals[i]])
   generator.applyBeamCrystalConvolution(nbins, xvals, yvals)

   # convert the result to a spectrum vs photon energy
   global dRtdkH1
   dRtdkH1 = 0
   title = "\\mbox{photon beam spectrum vs }E_\\gamma" + \
           "\\mbox{ (/GeV/s)}"
   dRtdkH1 = TH1D("dRtdkH1", title, nbins, Emin, Emax)
   dRtdkH1.GetXaxis().SetRangeUser(Emin + (Emax - Emin)/10., Emax)
   for i in range(0, nbins):
      dRtdkH1.Fill(xvals[i] * E0, yvals[i] / E0)
   dRtdkH1.SetStats(0)
   dRtdkH1.Draw("c")

   # print the rates in the three integration windows
   persec = (Emax - Emin) * 1./nbins
   prate = dRtdkH1.Integral(dRtdkH1.FindBin(pElim0),
                            dRtdkH1.FindBin(pElim1) - 1) * persec
   brate = dRtdkH1.Integral(dRtdkH1.FindBin(bElim0),
                            dRtdkH1.FindBin(bElim1) - 1) * persec
   erate = dRtdkH1.Integral(dRtdkH1.FindBin(eElim0),
                            dRtdkH1.FindBin(eElim1) - 1) * persec
   print "beam rate in the peak window [" + str(pElim0) + "," + \
         str(pElim1) + "] GeV = ", prate, "/s"
   print "beam rate in the background window [" + str(bElim0) + "," + \
         str(bElim1) + "] GeV =", brate, "/s"
   print "beam rate in the endpoint window [" + str(eElim0) + "," + \
         str(eElim1) + "] GeV =", erate, "/s"

def plotCoherent(collimated=1):
   """
   Plot the coherent part of the spectrum, either pre-collimator
   (collimated=0) or post-collimator (collimated=1) as a TH1D object.
   """

   # apply the beam-crystal convolution
   x0 = Emin / E0
   x1 = Emax / E0
   xvals = array('d', nbins * [0])
   yvals = array('d', nbins * [0])
   colFlag = generator.getCollimatedFlag()
   generator.setCollimatedFlag(collimated)
   generator.setPolarizedFlag(false)
   for i in range(0, nbins):
      xvals[i] = x0 + (i + 0.5) * (x1 - x0) / nbins;
      yvals[i] = dRcdx([xvals[i]])
   generator.setCollimatedFlag(colFlag)
   generator.applyBeamCrystalConvolution(nbins, xvals, yvals)

   # convert the result to a spectrum vs photon energy
   global dRcdkH1
   dRcdkH1 = 0
   title = "\\mbox{coherent photon beam spectrum vs }E_\\gamma" + \
           "\\mbox{ (/GeV/s)}"
   dRcdkH1 = TH1D("dRcdkH1", title, nbins, Emin, Emax)
   for i in range(0, nbins):
      dRcdkH1.Fill(xvals[i] * E0, yvals[i] / E0)
   dRcdkH1.SetStats(0)
   dRcdkH1.Draw("c")
   return dRcdkH1

def plotIncoherent(collimated=1):
   """
   Plot the incoherent part of the spectrum, either pre-collimator
   (collimated=0) or post-collimator (collimated=1) as a TH1D object.
   """

   # apply the beam-crystal convolution
   x0 = Emin / E0
   x1 = Emax / E0
   xvals = array('d', nbins * [0])
   yvals = array('d', nbins * [0])
   colFlag = generator.getCollimatedFlag()
   generator.setCollimatedFlag(collimated)
   generator.setPolarizedFlag(false)
   for i in range(0, nbins):
      xvals[i] = x0 + (i + 0.5) * (x1 - x0) / nbins;
      yvals[i] = dRidx([xvals[i]])
   generator.setCollimatedFlag(colFlag)
   generator.applyBeamCrystalConvolution(nbins, xvals, yvals)

   # convert the result to a spectrum vs photon energy
   global dRidkH1
   dRidkH1 = 0
   title = "\\mbox{incoherent photon beam spectrum vs }E_\\gamma" + \
           "\\mbox{ (/GeV/s)}"
   dRidkH1 = TH1D("dRidkH1", title, nbins, Emin, Emax)
   dRidkH1.GetXaxis().SetRangeUser(Emin + (Emax - Emin)/10., Emax)
   for i in range(0, nbins):
      dRidkH1.Fill(xvals[i] * E0, yvals[i] / E0)
   dRidkH1.SetStats(0)
   dRidkH1.Draw("c")
   return dRidkH1

def plotTotal(collimated=1):
   """
   Plot the total coherent+incoherent spectrum, either pre-collimator
   (collimated=0) or post-collimator (collimated=1) as a TH1D object.
   """

   # apply the beam-crystal convolution
   x0 = Emin / E0
   x1 = Emax / E0
   xvals = array('d', nbins * [0])
   yvals = array('d', nbins * [0])
   colFlag = generator.getCollimatedFlag()
   generator.setCollimatedFlag(collimated)
   generator.setPolarizedFlag(false)
   for i in range(0, nbins):
      xvals[i] = x0 + (i + 0.5) * (x1 - x0) / nbins;
      yvals[i] = dRtdx([xvals[i]])
   generator.setCollimatedFlag(colFlag)
   generator.applyBeamCrystalConvolution(nbins, xvals, yvals)

   # convert the result to a spectrum vs photon energy
   global dRtdkH1
   dRtdkH1 = 0
   title = "\\mbox{total photon beam spectrum vs }E_\\gamma" + \
           "\\mbox{ (/GeV/s)}"
   dRtdkH1 = TH1D("dRtdkH1", title, nbins, Emin, Emax)
   dRtdkH1.GetXaxis().SetRangeUser(Emin + (Emax - Emin)/10., Emax)
   for i in range(0, nbins):
      dRtdkH1.Fill(xvals[i] * E0, yvals[i] / E0)
   dRtdkH1.SetStats(0)
   dRtdkH1.Draw("c")
   return dRtdkH1

def plotPolarization(collimated=1):
   """
   Plot the linear polarization spectrum, either pre-collimator
   (collimated=0) or post-collimator (collimated=1) as a TH1D object.
   """

   # apply the beam-crystal convolution
   x0 = Emin / E0
   x1 = Emax / E0
   xvals = array('d', nbins * [0])
   yvals = array('d', nbins * [0])
   ypols = array('d', nbins * [0])
   colFlag = generator.getCollimatedFlag()
   generator.setCollimatedFlag(collimated)
   for i in range(0, nbins):
      xvals[i] = x0 + (i + 0.5) * (x1 - x0) / nbins;
      generator.setPolarizedFlag(true)
      ypols[i] = dRcdx([xvals[i]])
      generator.setPolarizedFlag(false)
      yvals[i] = dRtdx([xvals[i]])
   generator.setCollimatedFlag(colFlag)
   generator.applyBeamCrystalConvolution(nbins, xvals, yvals)
   generator.applyBeamCrystalConvolution(nbins, xvals, ypols)

   # convert the result to a spectrum vs photon energy
   global polarH1
   polarH1 = 0
   title = "\\mbox{photon beam polarization vs }E_\\gamma" + \
           "\\mbox{ (/GeV/s)}"
   polarH1 = TH1D("polarH1", title, nbins, Emin, Emax)
   for i in range(0, nbins):
      polarH1.Fill(xvals[i] * E0, ypols[i] / yvals[i])
   polarH1.SetStats(0)
   polarH1.Draw("c")
   return polarH1

def acceptance(vars):
   return generator.Acceptance(vars[0] ** 2)

def polarization(vars):
   return generator.Polarization(vars[0])

def dRtdx(vars):
   return generator.Rate_dNtdx(vars[0]) * cur / 1.6e-13

def dRcdx(vars):
   return generator.Rate_dNcdx(vars[0]) * cur / 1.6e-13

def dRidx(vars):
   return generator.Rate_dNidx(vars[0]) * cur / 1.6e-13

def plotTotal_rc(rchist, collimated=1):
   """
   Plot the total coherent+incoherent spectrum, either pre-collimator
   (collimated=0) or post-collimator (collimated=1) as a TH1D object,
   using a measured rocking curve in the place of an assumed Gaussian.
   """

   saved_thetax = generator.getTargetThetax()
   global dRtdkRC
   dRtdkRC = 0
   nsteps = rchist.GetNbinsX()
   threshold = 0.02 * rchist.GetMaximum()
   mean_tilt = rchist.GetMean()
   mean_thetax = generator.getTargetThetax()
   sum_intens = 0
   for step in range(1, nsteps + 1):
      intens = rchist.GetBinContent(step)
      if intens < threshold:
         continue
      tilt = rchist.GetXaxis().GetBinCenter(step) - mean_tilt
      generator.setTargetThetax(mean_thetax + tilt * 1e-6)
      dRtdkH1 = plotTotal(collimated)      
      if not dRtdkRC:
         dRtdkRC = dRtdkH1.Clone("dRtdkRC")
         dRtdkRC.Scale(0)
      dRtdkRC.Add(dRtdkH1, intens)
      sum_intens += intens
      dRtdkRC.Draw("c")
      dRtdkH1.Delete()
      c1.Update()
   generator.setTargetThetax(saved_thetax)
   dRtdkRC.Scale(1 / sum_intens)
   dRtdkRC.SetStats(0)
   dRtdkRC.Draw("c")
   return dRtdkRC

def plotPolarization_rc(rchist, collimated=1):
   """
   Plot the linear polarization spectrum, either pre-collimator
   (collimated=0) or post-collimator (collimated=1) as a TH1D object,
   using a measured rocking curve in the place of an assumed Gaussian.
   """

   saved_thetax = generator.getTargetThetax()
   global polarRC
   polarRC = 0
   nsteps = rchist.GetNbinsX()
   threshold = 0.02 * rchist.GetMaximum()
   mean_tilt = rchist.GetMean()
   mean_thetax = generator.getTargetThetax()
   sum_intens = 0
   for step in range(1, nsteps + 1):
      intens = rchist.GetBinContent(step)
      if intens < threshold:
         continue
      tilt = rchist.GetXaxis().GetBinCenter(step) - mean_tilt
      generator.setTargetThetax(mean_thetax + tilt * 1e-6)
      polarH1 = plotPolarization(collimated)      
      if not polarRC:
         polarRC = polarH1.Clone("polarRC")
         polarRC.Scale(0)
      polarRC.Add(polarH1, intens)
      sum_intens += intens
      polarRC.Draw("c")
      polarH1.Delete()
      c1.Update()
   generator.setTargetThetax(saved_thetax)
   polarRC.Scale(1 / sum_intens)
   polarRC.SetStats(0)
   polarRC.Draw("c")
   return polarRC
