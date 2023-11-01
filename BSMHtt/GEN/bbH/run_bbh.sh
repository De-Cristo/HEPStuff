#!/bin/bash

python ./run_pwg_condor.py -p 0 -i m125/bbH_m125.input -m bbH -f bbH_m125_output -d 1

# step 1-1
python ./run_pwg_condor.py -p 1 -x 1 -i m125/bbH_m125.input -m bbH -f bbH_m125_output -q longlunch -j 10
# step 1-2
python ./run_pwg_condor.py -p 1 -x 2 -i m125/bbH_m125.input -m bbH -f bbH_m125_output -q longlunch -j 10
# step 1-3
python ./run_pwg_condor.py -p 1 -x 3 -i m125/bbH_m125.input -m bbH -f bbH_m125_output -q longlunch -j 10
# step 1-4
python ./run_pwg_condor.py -p 1 -x 4 -i m125/bbH_m125.input -m bbH -f bbH_m125_output -q longlunch -j 10
# step 1-5
python ./run_pwg_condor.py -p 1 -x 5 -i m125/bbH_m125.input -m bbH -f bbH_m125_output -q longlunch -j 10

# ...
# step 2
python ./run_pwg_condor.py -p 2 -i m125/bbH_m125.input -m bbH -f bbH_m125_output -q longlunch -j 10
# step 3
python ./run_pwg_condor.py -p 3 -i m125/bbH_m125.input -m bbH -f bbH_m125_output -q longlunch -j 10

python ./run_pwg_condor.py -p 9 -i m125/bbH_m125.input -m bbH -f bbH_m125_output -k 1



# python ./run_pwg_condor.py -p 0 -i m100/bbH_m100.input -m bbH -f bbH_m100_output -d 1

# # step 1-1
# python ./run_pwg_condor.py -p 1 -x 1 -i m100/bbH_m100.input -m bbH -f bbH_m100_output -q longlunch -j 10
# # step 1-2
# python ./run_pwg_condor.py -p 1 -x 2 -i m100/bbH_m100.input -m bbH -f bbH_m100_output -q longlunch -j 10
# # step 1-3
# python ./run_pwg_condor.py -p 1 -x 3 -i m100/bbH_m100.input -m bbH -f bbH_m100_output -q longlunch -j 10
# # step 1-4
# python ./run_pwg_condor.py -p 1 -x 4 -i m100/bbH_m100.input -m bbH -f bbH_m100_output -q longlunch -j 10
# # step 1-5
# python ./run_pwg_condor.py -p 1 -x 5 -i m100/bbH_m100.input -m bbH -f bbH_m100_output -q longlunch -j 10

# # ...
# # step 2
# python ./run_pwg_condor.py -p 2 -i m100/bbH_m100.input -m bbH -f bbH_m100_output -q longlunch -j 10
# # step 3
# python ./run_pwg_condor.py -p 3 -i m100/bbH_m100.input -m bbH -f bbH_m100_output -q longlunch -j 10

# python ./run_pwg_condor.py -p 9 -i m100/bbH_m100.input -m bbH -f bbH_m100_output -k 1

mkdir bbH_slc7_amd64_gcc700_CMSSW_10_6_29_bbH_m125_output
cd bbH_slc7_amd64_gcc700_CMSSW_10_6_29_bbH_m125_output
tar xvzf ../bbH_slc7_amd64_gcc700_CMSSW_10_6_29_bbH_m125_output.tgz
./runcmsgrid.sh 10 1 4
