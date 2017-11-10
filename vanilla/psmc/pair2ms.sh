#!/bin/bash


PYTHON=/usr/local/var/homebrew/linked/python3/bin/python3.6


${PYTHON} ./pair2ms.py pairs_result_10_NN.ccf.txt ../vanilla.hdf5
${PYTHON} ./pair2ms.py pairs_result_10_RD.ccf.txt ../vanilla.hdf5


