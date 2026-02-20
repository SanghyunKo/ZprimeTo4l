#!/bin/sh

pwddir=$PWD

j="${1}"
isZA="${2}"

hmass=("250" "275" "300" "325" "350" "375" "400" "425" "450" "500" "550" "650" "750" "850" "1000" "1250" "1500" "1750" "2000")
i="${hmass[${3}]}" # ${hmass[${3}]}

cd /u/user/sako/ModHEEP/CMSSW_14_1_9/work/combined/verCWR_v2/${isZA}${j}_mu/${i}

ls /cvmfs/cms.cern.ch
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`

for quant in 0.025 0.16 0.5 0.84 0.975
do
  combine -d combined.root -M HybridNew --LHCmode LHC-limits --grid=higgsCombine.X${i}${isZA}${j}.mH${i}.root --clsAcc 0 -T 8000 -n .X${i}${isZA}${j} -m ${i} --readHybridResults --expectedFromGrid ${quant}
done

combine -d combined.root -M HybridNew --LHCmode LHC-limits --grid=higgsCombine.X${i}${isZA}${j}.mH${i}.root --clsAcc 0 -T 8000 -n .X${i}${isZA}${j} -m ${i} --readHybridResults

cd $pwddir

