#!/bin/bash



./rvage ibd -o ibd_plot_tru -i truH.bin --kList 3 --region 18875000 18885000 -m hmm --hmm HMM.initial.prob.txt HMM.emission.prob.txt -b 1000 --writeIBD --Ne 7300
./rvage ibd -o ibd_plot_tru -i truH.bin --kList 3 --region 18875000 18885000 -m fgt -b 1000 --writeIBD
./rvage ibd -o ibd_plot_tru -i truH.bin --kList 3 --region 18875000 18885000 -m dhg -b 1000 --writeIBD

./rvage ibd -o ibd_plot_err -i errH.bin --kList 3 --region 18875000 18885000 -m hmm --hmm HMM.initial.prob.txt HMM.emission.prob.txt -b 1000 --writeIBD --Ne 7300
./rvage ibd -o ibd_plot_err -i errH.bin --kList 3 --region 18875000 18885000 -m fgt -b 1000 --writeIBD
./rvage ibd -o ibd_plot_err -i errH.bin --kList 3 --region 18875000 18885000 -m dhg -b 1000 --writeIBD


