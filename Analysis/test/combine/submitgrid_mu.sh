#!/bin/sh

pwddir=$PWD

myarr=()
isZA="${2}"

if [ "${1}" = "1st" ]; then
  myarr=("0p4" "0p6" "0p8" "1" "1p5") # "0p4" "0p6" "0p8" "1" "1p5"
elif [ "${1}" = "2nd" ]; then
  myarr=("2" "5" "10" "50" "100") # "2" "5" "10" "50" "100"
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

    cd ${isZA}${j}_mu

    for i in ${myarr2[@]}
    do
      cd ${i}

      if [[ ${i} -le 1000 ]]; then
        combineTool.py -M HybridNew combined.root --LHCmode LHC-limits --clsAcc 0 -T 8000 -s -1 -n .grid1 --singlePoint 10:100:2 --merge 2 --saveToys --saveHybridResult -m ${i} --job-mode condor --task-name X${i}${isZA}${j}grid1 --sub-opts='RequestMemory = 4 GB' --rMin=-50 --rMax=100
      fi

      combineTool.py -M HybridNew combined.root --LHCmode LHC-limits --clsAcc 0 -T 8000 -s -1 -n .grid2 --singlePoint 1:9.8:0.2 --merge 2 --saveToys --saveHybridResult -m ${i} --job-mode condor --task-name X${i}${isZA}${j}grid2 --sub-opts='RequestMemory = 4 GB' --rMin=-5 --rMax=10

      if [[ ${i} -ge 500 ]]; then
        combineTool.py -M HybridNew combined.root --LHCmode LHC-limits --clsAcc 0 -T 8000 -s -1 -n .grid3 --singlePoint 0.1:0.98:0.02 --merge 2 --saveToys --saveHybridResult -m ${i} --job-mode condor --task-name X${i}${isZA}${j}grid3 --sub-opts='RequestMemory = 8 GB' --rMin=-0.5 --rMax=1
      fi

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
  cd ${isZA}${j}_mu

  for i in ${hmass[@]}
  do
    cd ${i}

    jtemp=`echo ${j} | tr p .`

    boolstr=`echo "${jtemp} < 1" | bc -l`

    if [ "${boolstr}" = "1" ]; then

      if [[ ${i} -le 1000 ]]; then
        combineTool.py -M HybridNew combined.root --LHCmode LHC-limits --clsAcc 0 -T 8000 -s -1 -n .gridExt --singlePoint 1200:10000:200 --merge 2 --saveToys --saveHybridResult -m ${i} --job-mode condor --task-name X${i}${isZA}${j}gridExt --sub-opts='RequestMemory = 4 GB' --rMin=-100 --rMax=10000

	combineTool.py -M HybridNew combined.root --LHCmode LHC-limits --clsAcc 0 -T 8000 -s -1 -n .grid0 --singlePoint 120:1000:20 --merge 2 --saveToys --saveHybridResult -m ${i} --job-mode condor --task-name X${i}${isZA}${j}grid0 --sub-opts='RequestMemory = 4 GB' --rMin=-100 --rMax=1000
      fi
    fi

    combineTool.py -M HybridNew combined.root --LHCmode LHC-limits --clsAcc 0 -T 8000 -s -1 -n .grid1 --singlePoint 10:100:2 --merge 2 --saveToys --saveHybridResult -m ${i} --job-mode condor --task-name X${i}${isZA}${j}grid1 --sub-opts='RequestMemory = 4 GB' --rMin=-50 --rMax=100

    combineTool.py -M HybridNew combined.root --LHCmode LHC-limits --clsAcc 0 -T 8000 -s -1 -n .grid2 --singlePoint 1:9.8:0.2 --merge 2 --saveToys --saveHybridResult -m ${i} --job-mode condor --task-name X${i}${isZA}${j}grid2 --sub-opts='RequestMemory = 4 GB' --rMin=-5 --rMax=10

    if [[ ${i} -ge 500 ]]; then
      combineTool.py -M HybridNew combined.root --LHCmode LHC-limits --clsAcc 0 -T 8000 -s -1 -n .grid3 --singlePoint 0.1:0.98:0.02 --merge 2 --saveToys --saveHybridResult -m ${i} --job-mode condor --task-name X${i}${isZA}${j}grid3 --sub-opts='RequestMemory = 8 GB' --rMin=-0.5 --rMax=1
    fi


    cd ..
  done

  cd ..
done

cd ${pwddir}

