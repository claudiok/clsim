# FIXME: The functionality of this script is now integrated into env-shell.sh,
#        therefore this script should be removed at some point to avoid confusion.
#
# Copyright (c) 2011, 2012
# Claudio Kopper <claudio.kopper@icecube.wisc.edu>
# and the IceCube Collaboration <http://www.icecube.wisc.edu>
# 
# Permission to use, copy, modify, and/or distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
# 
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
# SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION
# OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
# CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
# 
# 
# $Id: AutoSetGeant4Environment.py 137653 2015-09-16 22:30:32Z hdembinski $
# 
# @file AutoSetGeant4Environment.py
# @version $Revision: 137653 $
# @date $Date: 2015-09-16 18:30:32 -0400 (Wed, 16 Sep 2015) $
# @author Claudio Kopper
#

from __future__ import print_function
from distutils import spawn

import os
import subprocess
import pickle

def AutoSetGeant4Environment(force=True):
    # older versions used these hard-coded locations:
    hardCodedForGeant4_9_3and4 = \
        {"G4LEVELGAMMADATA" : os.path.expandvars("$SROOT/share/geant4/data/PhotonEvaporation2.0"),
         "G4RADIOACTIVEDATA" : os.path.expandvars("$SROOT/share/geant4/data/RadioactiveDecay3.2"),
         "G4LEDATA" : os.path.expandvars("$SROOT/share/geant4/data/G4EMLOW6.9"),
         "G4NEUTRONHPDATA" : os.path.expandvars("$SROOT/share/geant4/data/G4NDL3.13"),
         "G4ABLADATA" : os.path.expandvars("$SROOT/share/geant4/data/G4ABLA3.0"),
         "G4NEUTRONXSDATA" : "(none)", # geant 4.9.3/4 does not have these
         "G4PIIDATA" : "(none)",       # geant 4.9.3/4 does not have these
         "G4REALSURFACEDATA" : "(none)",# geant 4.9.3/4 does not have these
         "G4SAIDXSDATA": "(none)", # geant 4.9.3/4 does not have these
         "G4ENSDFSTATEDATA": "(none)", # geant 4.9.3/4 does not have these
        }
    hasOldGeant4 = (os.path.isdir(os.path.expandvars("$I3_PORTS/lib/geant4_4.9.4")) or os.path.isdir(os.path.expandvars("$I3_PORTS/lib/geant4_4.9.3"))) and (not os.path.isfile(os.path.expandvars("$I3_PORTS/bin/geant4.sh")))
    
    if not hasOldGeant4:
        
        if os.path.isfile(os.path.expandvars("$I3_PORTS/bin/geant4.sh")):
            geant4_script = os.path.expandvars("$I3_PORTS/bin/geant4.sh")
        elif os.path.isfile(os.path.expandvars("$SROOT/bin/geant4.sh")):
            geant4_script = os.path.expandvars("$SROOT/bin/geant4.sh")
        elif spawn.find_executable("geant4-config"):
            pfx = os.path.normpath(subprocess.check_output(["geant4-config", "--prefix"]).strip())
            geant4_script = os.path.normpath("{}/bin/geant4.sh".format(pfx))
        elif force:
            print("cannot forcibly overwrite the Geant4 environment variables because geant4.sh is not available. Using what's available.")
            force = False
            geant4_script = ''
    
    Geant4Variables = set(hardCodedForGeant4_9_3and4.keys())
        
    Geant4Variables_unset = set()
    for var in Geant4Variables:
        if var not in os.environ:
            Geant4Variables_unset.add(var)
        
    Geant4Variables_set = Geant4Variables.difference(Geant4Variables_unset)
    
    # if the user chose to force overwrite everything, reset the "unset" set to the full set
    if force:
        Geant4Variables_setnew = Geant4Variables
    else:
        Geant4Variables_setnew = Geant4Variables_unset
    
    if len(Geant4Variables_unset)>0:
        if not hasOldGeant4:
            print("Not all Geant4 environment variables are set. Trying to get some of them from geant4.sh..")
        else:
            print("Not all Geant4 environment variables are set. Trying to use defaults for geant4.9.3/4.9.4..")
            
        if force:
            print("already set: (will be overwritten)")
        else:
            print("already set:")
        for var in Geant4Variables_set:
            print("  *", var, "->", os.environ[var])

        if not force:
            print("missing:")
            for var in Geant4Variables_unset:
                print("  *", var)
                
        if hasOldGeant4:
            geant4env = hardCodedForGeant4_9_3and4
        else:
            if not os.path.isfile(geant4_script):
                raise RuntimeError("Cannot automatically set missing environment variables. ($I3_PORTS/bin/geant4.sh is missing.) Please set them yourself.")

            # get the environment after loading geant4.sh
            #source = os.path.expandvars(". $I3_PORTS/bin/geant4.sh")
            g4scriptloc = geant4_script
            #dump = '/usr/bin/env python -c "import os,pickle;print pickle.dumps(os.environ)"'
            penv = subprocess.Popen(['/bin/bash', '-c', '. ' + g4scriptloc + ' && /usr/bin/env python -c "import os,pickle;print pickle.dumps(os.environ)"'], stdout=subprocess.PIPE).communicate()[0]
            #penv = subprocess.Popen('%s && %s' %(source,dump), stdout=subprocess.PIPE, shell=True).communicate()[0]
            #penv = subprocess.check_output('%s && %s' %(source,dump), shell=True)
            geant4env = pickle.loads(penv)
            #penv = os.popen('%s && %s' %(source,dump))
            #geant4env = pickle.loads(penv.read())

        print("setting from geant4.sh:")
        for var in Geant4Variables_setnew:
            if var not in geant4env:
                if hardCodedForGeant4_9_3and4[var] != '(none)':
                    raise RuntimeError("Cannot find the %s environment variable in the geant4.sh script." % var)
            else:
                os.environ[var] = geant4env[var]
                print("  *", var, "->", os.environ[var])

    
