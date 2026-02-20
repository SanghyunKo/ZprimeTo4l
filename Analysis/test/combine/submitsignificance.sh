#!/bin/sh

pwddir=$PWD

myarr=()
isZA="${2}"

if [ "${1}" = "1st" ]; then
  myarr=("0p4" "0p6" "0p8" "1" "1p5") # "0p4" "0p6" "0p8" "1" "1p5"
elif [ "${1}" = "2nd" ]; then
  myarr=("2" "5" "10" "50" "100")
elif [ "${1}" = "3rd" ]; then
  myarr=("250" "500" "750" "1000" "1500")

  for j in ${myarr[@]}
  do
    myarr2=()

    if [ "${j}" = "250" ]; then
      if [ "${isZA}" = "A" ]; then
        myarr2=("750" "1000" "1500")
      elif [ "${isZA}" = "Z" ]; then
        myarr2=("500" "750")
      fi
    elif [ "${j}" = "500" ]; then
      if [ "${isZA}" = "A" ]; then
        myarr2=("1500")
      elif [ "${isZA}" = "Z" ]; then
        myarr2=("750" "1000" "1500" "2000")
      fi
    elif [ "${j}" = "750" ]; then
      if [ "${isZA}" = "A" ]; then
        myarr2=("2000")
      elif [ "${isZA}" = "Z" ]; then
        myarr2=("1000")
      fi
    elif [ "${j}" = "1000" ]; then
      if [ "${isZA}" = "Z" ]; then
        myarr2=("1500")
      fi
    elif [ "${j}" = "1500" ]; then
      if [ "${isZA}" = "Z" ]; then
        myarr2=("2000")
      fi
    fi

    cd ${isZA}${j}

    for i in ${myarr2[@]}
    do
      cd ${i}

      rmax=100

      if [[ ${i} -gt 1000 ]]; then
        rmax=10
      fi

      combineTool.py --LHCmode LHC-significance --saveToys --fullBToys --saveHybridResult --rMax=${rmax} -T 4000 -i 1 -M HybridNew -m ${i} -s -1 -d combined.root -n .significanceCWR.X${i}${isZA}${j}.sig1 --job-mode condor --task-name X${i}${isZA}${j}sig1 --sub-opts='RequestMemory = 4 GB'

      combineTool.py --LHCmode LHC-significance --saveToys --fullBToys --saveHybridResult --rMax=${rmax} -T 4000 -i 1 -M HybridNew -m ${i} -s -1 -d combined.root -n .significanceCWR.X${i}${isZA}${j}.sig2 --job-mode condor --task-name X${i}${isZA}${j}sig2 --sub-opts='RequestMemory = 4 GB'

      combineTool.py --LHCmode LHC-significance --saveToys --fullBToys --saveHybridResult --rMax=${rmax} -T 4000 -i 1 -M HybridNew -m ${i} -s -1 -d combined.root -n .significanceCWR.X${i}${isZA}${j}.sig3 --job-mode condor --task-name X${i}${isZA}${j}sig3 --sub-opts='RequestMemory = 4 GB'

      combineTool.py --LHCmode LHC-significance --saveToys --fullBToys --saveHybridResult --rMax=${rmax} -T 4000 -i 1 -M HybridNew -m ${i} -s -1 -d combined.root -n .significanceCWR.X${i}${isZA}${j}.sig4 --job-mode condor --task-name X${i}${isZA}${j}sig4 --sub-opts='RequestMemory = 4 GB'

      combineTool.py --LHCmode LHC-significance --saveToys --fullBToys --saveHybridResult --rMax=${rmax} -T 4000 -i 1 -M HybridNew -m ${i} -s -1 -d combined.root -n .significanceCWR.X${i}${isZA}${j}.sig5 --job-mode condor --task-name X${i}${isZA}${j}sig5 --sub-opts='RequestMemory = 4 GB'

      cd ..
    done

    cd ..
  done

  cd ${pwddir}
  exit 0
fi

hmass=("250" "275" "300" "325" "350" "375" "400" "425" "450" "500" "550" "650" "750" "850" "1000" "1250" "1500" "1750" "2000")

for j in ${myarr[@]}
do
  cd ${isZA}${j}

  for i in ${hmass[@]}
  do
    cd ${i}

    rmax=100

    if [[ ${i} -gt 1000 ]]; then
      rmax=10
    fi

    ntoy=4000

    if [[ ${i} -eq 550 ]]; then
      if [[ "${j}" == "0p8" ]]; then
        ntoy=8000
      fi
    fi

    combineTool.py --LHCmode LHC-significance --saveToys --fullBToys --saveHybridResult --rMax=${rmax} -T ${ntoy} -i 1 -M HybridNew -m ${i} -s -1 -d combined.root -n .significanceCWR.X${i}${isZA}${j}.sig1 --job-mode condor --task-name X${i}${isZA}${j}sig1 --sub-opts='RequestMemory = 4 GB'

    combineTool.py --LHCmode LHC-significance --saveToys --fullBToys --saveHybridResult --rMax=${rmax} -T ${ntoy} -i 1 -M HybridNew -m ${i} -s -1 -d combined.root -n .significanceCWR.X${i}${isZA}${j}.sig2 --job-mode condor --task-name X${i}${isZA}${j}sig2 --sub-opts='RequestMemory = 4 GB'

    combineTool.py --LHCmode LHC-significance --saveToys --fullBToys --saveHybridResult --rMax=${rmax} -T ${ntoy} -i 1 -M HybridNew -m ${i} -s -1 -d combined.root -n .significanceCWR.X${i}${isZA}${j}.sig3 --job-mode condor --task-name X${i}${isZA}${j}sig3 --sub-opts='RequestMemory = 4 GB'

    combineTool.py --LHCmode LHC-significance --saveToys --fullBToys --saveHybridResult --rMax=${rmax} -T ${ntoy} -i 1 -M HybridNew -m ${i} -s -1 -d combined.root -n .significanceCWR.X${i}${isZA}${j}.sig4 --job-mode condor --task-name X${i}${isZA}${j}sig4 --sub-opts='RequestMemory = 4 GB'

    combineTool.py --LHCmode LHC-significance --saveToys --fullBToys --saveHybridResult --rMax=${rmax} -T ${ntoy} -i 1 -M HybridNew -m ${i} -s -1 -d combined.root -n .significanceCWR.X${i}${isZA}${j}.sig5 --job-mode condor --task-name X${i}${isZA}${j}sig5 --sub-opts='RequestMemory = 4 GB'

    cd ..
  done

  cd ..
done

cd ${pwddir}

