#!/bin/sh

cd /u/user/sako/ModHEEP/CMSSW_14_1_9/work/combined/verCWR_v2/Z0p8/550/

ls /cvmfs/cms.cern.ch
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`

input=${1} # proc starts from 0
ntoy=$((${input}+1)) # ntoy starts from 1

hmass=("250"" 275" "300" "325" "350" "375" "400" "425" "450" "500" "550" "650" "750" "850" "1000" "1250" "1500" "1750" "2000")
#hmass=("250")

cd filesLEE

for mass in ${hmass[@]}
do

  rmax=100

  if [[ ${mass} -gt 1000 ]]; then
    rmax=10
  fi

  combine -M HybridNew ../../${mass}/combined.root --LHCmode LHC-significance --saveToys --fullBToys --saveHybridResult --rMax=${rmax} -m ${mass} -n .LEEmX${mass}toy${ntoy} -T 2000 -i 1 -D higgsCombineTest.GenerateOnly.mH550.123456.root:toys/toy_${ntoy} --pvalue
done

