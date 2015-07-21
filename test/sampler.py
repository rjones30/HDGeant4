#!/usr/bin/python
#
# sampler.py - utility functions for comparing hdds simulation geometries
#              in hdgeant (geant3) and hdgeant4 (geant4).
#
# author: richard.t.jones at uconn.edu
# version: july 21, 2015
#
# usage pattern:
#

import random
import sys
import re

def usage():
   print "Usage: sampler.py <command>"
   print "  where <command> is one of:"
   print "  *) sample - generate a set of random coordinates for sampling"
   print "  *) scan3 <logfile> - scan logfile for picking output from hdgeant"
   print "  *) scan4 <logfile> - scan logfile for picking output from hdgeant4"
   print "  *) compare <file3> <file4> - compare output from file3 and file4"
   print "     which should be output from prior scan3 and scan4 runs."

def sample():
   random.seed("It's a wonderful day in the neighborhood.")
   for n in range(0,100000):
      s = "{0:7.2f} ".format(random.uniform(-100,100)) \
        + "{0:7.2f} ".format(random.uniform(-100,100)) \
        + "{0:7.2f}".format(random.uniform(-2400,800))
      print s

def scan4(infile):
   fin = open(infile)
   fld = 1
   volname = "NONE"
   for line in fin:
      m1 = re.match(r"\(([-0-9.eE]+),([-0-9.eE]+),([-0-9.eE]+)\) +found in +([A-Za-z0-9]+).*", line)
      if m1:
         coords = 1
         if not fld:
           print "{0:7.2f} ".format(x), \
                 "{0:7.2f} ".format(y), \
                 "{0:7.2f}    ".format(z), \
                 volname, "0 0 0"
         fld = None
         x = float(m1.group(1))
         y = float(m1.group(2))
         z = float(m1.group(3))
         volname = m1.group(4)
      m2 = re.match(r" +magnetic field \(Tesla\): ([-0-9.eE]+),([-0-9.eE]+),([-0-9.eE]+)", line)
      if m2:
         fld = 1
         if not coords:
           print "missing coordinate line"
         coords = None
         Bx = float(m2.group(1))
         By = float(m2.group(2))
         Bz = float(m2.group(3))
         print "{0:7.2f} ".format(x), \
               "{0:7.2f} ".format(y), \
               "{0:7.2f}    ".format(z), \
               volname, \
               "{0:12.7f} ".format(Bx), \
               "{0:12.7f} ".format(By), \
               "{0:12.7f} ".format(Bz)

def scan3(infile):
   fin = open(infile)
   fld = 1
   volname = "NONE"
   for line in fin:
      m1 = re.match(r" *point *[*0-9]+: *([-0-9.eE]+) +([-0-9.eE]+) +([-0-9.eE]+).*/([ A-Z0-9]{4}) *[0-9]+/", line)
      if m1:
         coords = 1
         if not fld:
           print "{0:7.2f} ".format(x), \
                 "{0:7.2f} ".format(y), \
                 "{0:7.2f}    ".format(z), \
                 volname, "0 0 0"
         fld = None
         x = float(m1.group(1))
         y = float(m1.group(2))
         z = float(m1.group(3))
         volname = m1.group(4)
      m2 = re.match(r".* inhomogeneous field *\( *([-0-9.eE]+) *, *([-0-9.eE]+) *, *([-0-9.eE]+) *\) *kG.*", line)
      if m2:
         fld = 1
         if not coords:
           print "missing coords"
         coords = None
         Bx = float(m2.group(1))
         By = float(m2.group(2))
         Bz = float(m2.group(3))
         print "{0:7.2f} ".format(x), \
               "{0:7.2f} ".format(y), \
               "{0:7.2f}    ".format(z), \
               volname, \
               "{0:12.7f} ".format(Bx/10), \
               "{0:12.7f} ".format(By/10), \
               "{0:12.7f} ".format(Bz/10)

def scan34(infile3, infile4):
   ifx = open(infile3)
   ify = open(infile4)
   xline = []
   yline = []
   for line in ifx:
      xline.append(line)
   for line in ify:
      yline.append(line)
   for line in range(0,len(xline)):
      xwords = xline[line].split()
      ywords = yline[line].split()
      if len(xwords) < 3 or len(ywords) < 3:
         continue
      if xwords[0] == ywords[0] and xwords[1] == ywords[1] and \
                                    xwords[2] == ywords[2]:
         dx = float(xwords[4]) - float(ywords[4])
         dy = float(xwords[5]) - float(ywords[5])
         dz = float(xwords[6]) - float(ywords[6])
         dr2 = dx ** 2 + dy ** 2 + dz ** 2
         if xwords[3].upper() != ywords[3].upper():
           print "{0:7.2f} ".format(float(xwords[0])), \
                 "{0:7.2f} ".format(float(xwords[1])), \
                 "{0:7.2f} ".format(float(xwords[2])), \
                 xwords[3], "!=", ywords[3]
         elif dr2 > 1e-4:
           print "{0:7.2f} ".format(float(xwords[0])), \
                 "{0:7.2f} ".format(float(xwords[1])), \
                 "{0:7.2f} ".format(float(xwords[2])), \
                 "{0:9.7f} ".format(float(ywords[4])), \
                 "{0:9.7f} ".format(float(ywords[5])), \
                 "{0:9.7f} ".format(float(ywords[6])), \
                 "{0:9.7f} ".format(float(xwords[4])), \
                 "{0:9.7f} ".format(float(xwords[5])), \
                 "{0:9.7f} ".format(float(xwords[6]))

# parse commands from command line

if len(sys.argv) < 2:
   usage()
elif sys.argv[1] == "sample":
   sample()
elif sys.argv[1] == "scan3" and len(sys.argv) > 2:
   scan3(sys.argv[2])
elif sys.argv[1] == "scan4" and len(sys.argv) > 2:
   scan4(sys.argv[2])
elif sys.argv[1] == "compare" and len(sys.argv) > 3:
   scan34(sys.argv[2], sys.argv[3])
else:
   usage()
