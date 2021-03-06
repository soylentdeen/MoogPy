#!/usr/bin/python
import MoogTools
import sys

moogPyConfigFile = sys.argv[1]
flavor = sys.argv[2]
Moog = MoogTools.MoogStokes(moogPyConfigFile, fileBase=flavor, moogInstance='Alpha',  progressBar = True)
Moog.run(saveRaw=True)
