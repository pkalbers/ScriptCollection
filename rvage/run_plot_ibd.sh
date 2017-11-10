#!/bin/bash



./rvage ibd -o ibd_plot_tru -i truH.bin --kList 4 --region 28875000 28895000 -m hmm --hmm HMM.initial.prob.txt HMM.emission.prob.txt -b 1000 --writeIBD --Ne 7300
./rvage ibd -o ibd_plot_tru -i truH.bin --kList 4 --region 28875000 28895000 -m fgt -b 1000 --writeIBD 
./rvage ibd -o ibd_plot_tru -i truH.bin --kList 4 --region 28875000 28895000 -m dhg -b 1000 --writeIBD 

./rvage ibd -o ibd_plot_err -i errH.bin --kList 4 --region 28875000 28895000 -m hmm --hmm HMM.initial.prob.txt HMM.emission.prob.txt -b 1000 --writeIBD --Ne 7300
./rvage ibd -o ibd_plot_err -i errH.bin --kList 4 --region 28875000 28895000 -m fgt -b 1000 --writeIBD 
./rvage ibd -o ibd_plot_err -i errH.bin --kList 4 --region 28875000 28895000 -m dhg -b 1000 --writeIBD 


