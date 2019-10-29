#!/usr/bin/python

import random
import os
import subprocess
import math
import sys
import time
import copy
from numpy.random import choice as choices

def switch_nt(nt):
	switcher = {
		"A": list(choices(["C","T","G"],1,p=[0.14, 0.04, 0.82])),
		"T": list(choices(["C","A","G"],1,p=[0.84, 0.03, 0.13])),
		"G": list(choices(["C","A","T"],1,p=[0.19, 0.7, 0.11])),
		"C": list(choices(["G","A","T"],1,p=[0.17, 0.12, 0.71])),
		"N": ["N"]
	}
	return switcher.get(nt, "N")
