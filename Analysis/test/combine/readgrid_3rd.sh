#!/bin/sh

pwddir=$PWD

isZA="${1}"

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

  cd ${isZA}${j}_el

  for i in ${myarr2[@]}
  do
    cd ${i}

    for quant in 0.025 0.16 0.5 0.84 0.975
    do
      combine -d combined.root -M HybridNew --LHCmode LHC-limits --grid=higgsCombine.X${i}${isZA}${j}.mH${i}.root --clsAcc 0 -T 8000 -n .X${i}${isZA}${j} -m ${i} --readHybridResults --expectedFromGrid ${quant}
    done

    combine -d combined.root -M HybridNew --LHCmode LHC-limits --grid=higgsCombine.X${i}${isZA}${j}.mH${i}.root --clsAcc 0 -T 8000 -n .X${i}${isZA}${j} -m ${i} --readHybridResults

    cd ..
  done

  cd ..
done

cd ${pwddir}
exit 0



