#!/usr/bin/env python3.8

########################################################################

import sys
import quick_wavelet

infile1 = sys.argv[1]
infile2 = sys.argv[2]

quick_wavelet.run_double_wavelet_analysis(infile1, infile2, dt=10000., mirror=True, cut1=None, cut2=None, write_output=True, wf='morlet', dj=0.1, om0=6, normmean=True, mirrormethod=1)
