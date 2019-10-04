#this script is used to run geneLevel.R
#note it relies on us generating a list of commands first, we do that in iPython

import subprocess
import sys

geneLevelCommandsPath = '/ifs/work/taylorlab/friedman/myAdjustedDataFiles/geneLevelCommands.txt'
with open(geneLevelCommandsPath) as f:
  lines = f.readlines()
  for line in lines:
    proc = subprocess.Popen(line, shell=True)
    proc.wait()
