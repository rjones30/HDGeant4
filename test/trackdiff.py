#!/usr/bin/python
#
# trackdiff.py - script to compare the track stepping log output
#                from hdgeant and hdgeant4 (/tracking/verbose 2)
#
# author: richard.t.jones at uconn.edu
# version: august 10, 2015

import sys
import re

def usage():
   print
   print "Usage: trackdiff.py <g3.log> <g4.log>"
   sys.exit(1)

if len(sys.argv) != 3:
   usage()

g3log = open(sys.argv[1])
g4log = open(sys.argv[2])
if not (g3log and g4log):
   usage()

def read_g3log():
   """
   Read the next event from the geant3 tracking log and return the step
   information at each geometry transition line in a tuple object with
   the contents (step, name, num, x, y, z, ds, s), using generator semantics.
   """
   count = 0
   name = ""
   num = 0
   global g3log
   for line in g3log:
      m = re.match(r" +([-.0-9]+) +([-.0-9]+) +([-.0-9]+) +([.0-9]+)" +
                   r" +([A-Z0-9]+) +([0-9]+) +([.0-9]+) +([.0-9]+)" +
                   r" +([.0-9]+) +([A-Za-z]+) +([.0-9]+) +([A-Za-z]+) +(.*)",
                   line)
      if m:
         if float(m.group(7)) == 0:
            count = 0
            name = m.group(5)
            num = int(m.group(6))
            mlist = [m]
         elif m.group(5) != name or int(m.group(6)) != num:
            count += 1
            x = float(mlist[0].group(1))
            y = float(mlist[0].group(2))
            z = float(mlist[0].group(3))
            r = float(mlist[0].group(4))
            s = float(mlist[0].group(7))
            ds = sum(float(mi.group(8)) for mi in mlist)
            ds += float(m.group(8)) - float(mlist[0].group(8))
            dE = sum(energy_in_GeV(mi.group(9), mi.group(10)) for mi in mlist)
            Ek = energy_in_GeV(mlist[0].group(11), mlist[0].group(12))
            yield (count, name, num, x, y, z, ds, s)
            name = m.group(5)
            num = int(m.group(6))
            mlist = [m]
         else:
            mlist.append(m)
      elif re.match(r" *X *Y *Z *R *NAME *NUMBER *SLENG" +
                    r" *STEP *DESTEP *GEKIN *MECHANISMS",
                    line):
         return
   g3log = 0

def read_g4log():
   """
   Read the next event from the geant4 tracking log and return the step
   information at each geometry transition line in a tuple object with
   the contents (step, name, num, x, y, z, ds, s), using generator semantics.
   """
   global g4log
   name = ""
   num = 0
   for line in g4log:
      m = re.match(r" +([0-9]+) +([-.0-9]+) +([a-z]+) +([-.0-9]+) +([a-z]+)" +
                   r" +([-.0-9]+) +([a-z]+) +([.0-9]+) +([A-Za-z]+)" +
                   r" +([.0-9]+) +([A-Za-z]+) +([.0-9]+) +([a-z]+)" +
                   r" +([.0-9]+) +([a-z]+) +([:A-Za-z0-9]+):([0-9]+)" +
                   r" +([^ ].*)",
                   line)
      if m:
         if length_in_cm(m.group(14), m.group(15)) == 0:
            name = m.group(16)
            num = int(m.group(17))
            mlist = [m]
         elif m.group(16) != name or int(m.group(17)) != num:
            n = int(mlist[0].group(1))
            x = length_in_cm(mlist[0].group(2), mlist[0].group(3))
            y = length_in_cm(mlist[0].group(4), mlist[0].group(5))
            z = length_in_cm(mlist[0].group(6), mlist[0].group(7))
            Ek = energy_in_GeV(mlist[0].group(8), mlist[0].group(9))
            dE = sum(energy_in_GeV(mi.group(10), mi.group(11)) for mi in mlist)
            ds = sum(length_in_cm(mi.group(12), mi.group(13)) for mi in mlist)
            ds -= length_in_cm(mlist[0].group(12), mlist[0].group(13))
            ds += length_in_cm(m.group(12), m.group(13))
            s = length_in_cm(mlist[0].group(14), mlist[0].group(15))
            if ds > 1e-12:
               yield (n, name, num, x, y, z, ds, s)
            name = m.group(16)
            num = int(m.group(17))
            mlist = [m]
         else:
            mlist.append(m)
      elif re.match(r"Step# *X *Y *Z *KineE *dEStep *" +
                    r"StepLeng *TrakLeng *Volume *Process",
                    line):
         return
   g4log = 0

def length_in_cm(value, unit):
   """
   Returns length specified in some arbitrary length unit in cm.
   """
   if unit == "fm":
      return float(value) * 1e-13
   elif unit == "nm":
      return float(value) * 1e-7
   elif unit == "um":
      return float(value) * 1e-4
   elif unit == "mm":
      return float(value) * 1e-1
   elif unit == "cm":
      return float(value)
   elif unit == "m":
      return float(value) * 1e2
   elif unit == "km":
      return float(value) * 1e5
   else:
      print "unknown length unit", unit, "quitting!"
      sys.exit(1)

def energy_in_GeV(value, unit):
   """
   Returns energy specified in some arbitrary length unit in GeV.
   """
   if unit == "meV":
      return float(value) * 1e-12
   elif unit == "eV":
      return float(value) * 1e-9
   elif unit == "keV":
      return float(value) * 1e-6
   elif unit == "MeV":
      return float(value) * 1e-3
   elif unit == "GeV":
      return float(value)
   elif unit == "TeV":
      return float(value) * 1e3
   else:
      print "unknown energy unit", unit, "quitting!"
      sys.exit(1)

def compare(g3, g4):
   """
   Runs a sequential comparison against the two step sequences
   represented by iterables g3, g4 which return tuples of the
   form (step, name, num, x, y, z, ds, s) with generator semantics.
   """
   ig3 = iter(g3)
   ig4 = iter(g4)
   t3 = ig3.next()
   t4 = ig4.next()
   l3 = []
   l4 = []
   while True:
      g3name = t3[1] + ":" + str(t3[2])
      # Automatic division names in the g4 listing are in lowercase;
      # their indices need to be incremented by one to match the
      # division numbers in the g3 listing.
      if t4[1] != t4[1].upper():
         nc = t4[2] + 1
      else:
         nc = t4[2]
      # Double-colon layer suffixes on volume names in the g4 listing
      # must be stripped off to match the names in the g3 listing.
      m = re.match(r"([A-Za-z0-9]+)::.*", t4[1])
      if m:
         g4name = m.group(1) + ":" + str(nc)
      else:
         g4name = t4[1] + ":" + str(nc)
      # Convert the World volume in the g4 listing to SITE to match the
      # name for the root volume in the g3 geometry.
      if g4name == "World:1":
         g4name = "SITE:1"
      if abs(t3[7] - t4[7]) < 0.01 + (t3[7] + t4[7]) * 1e-3 and \
                              g3name.upper() == g4name.upper():
         if len(l3) + len(l4) > 0:
            for tl3 in l3:
               xyzs = [tl3[3], tl3[4], tl3[5], tl3[7], tl3[6]]
               print "  ", tl3[0], tl3[1], " ".join(["%9.4f" % x for x in xyzs])
            for tl4 in l4:
               xyzs = [tl4[3], tl4[4], tl4[5], tl4[7], tl4[6]]
               print "G4", tl4[0], tl4[1], " ".join(["%9.4f" % x for x in xyzs])
            xyzs = [t3[3], t3[4], t3[5], t3[7], t3[6]]
            print "  ", t3[0], g3name, " ".join(["%9.4f" % x for x in xyzs])
            xyzs = [t4[3], t4[4], t4[5], t4[7], t4[6]]
            print "G4", t4[0], g4name, " ".join(["%9.4f" % x for x in xyzs])
            print
            l3 = []
            l4 = []
         try:
            t3 = ig3.next()
            t4 = ig4.next()
            continue
         except:
            break
      else:
         try:
            if t3[7] < t4[7]:
               l3.append((t3[0], g3name, t3[2], t3[3], t3[4],
                                         t3[5], t3[6], t3[7]))
               t3 = ig3.next()
            else:
               l4.append((t4[0], g4name, t4[2], t4[3], t4[4],
                                         t4[5], t4[6], t4[7]))
               t4 = ig4.next()
         except:
            break;

   while True:
      try:
         ig3.next()
      except:
         break
   while True:
      try:
         ig4.next()
      except:
         break

# advance the input logs to the start of the first event
sum(0 for m in read_g3log())
sum(0 for m in read_g4log())

# loop over events, should be the same number in both logs
if False:
   while g3log and g4log:
      lines3 = sum(1 for m in read_g3log())
      lines4 = sum(1 for m in read_g4log())
      print sys.argv[1], "has", lines3, "lines,", sys.argv[2], "has", lines4, "lines"

if True:
   count = 0
   while g3log and g4log:
      count += 1
      print "**************** event", count, "********************"
      compare(read_g3log(), read_g4log())
