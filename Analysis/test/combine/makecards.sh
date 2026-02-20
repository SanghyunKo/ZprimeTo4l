#!/bin/sh

pwddir=$PWD

for j in 0p4 0p6 0p8 1 1p5 2 5 10 50 100
do
  mkdir A${j}
  commandstr="combineCards.py"
  jtemp=`echo ${j} | tr p .`

  boolstr=`echo "${jtemp} <= 10" | bc -l`

  if [ "${boolstr}" = "1" ]; then
    commandstr="${commandstr} mergedEle2E=datacard_ME2E_Y1.txt"
  fi

  boolstr=`echo "${jtemp} <= 2" | bc -l`

  if [ "${boolstr}" = "1" ]; then
    commandstr="${commandstr} mergedEMu1M=datacard_MEMu1M_Y1.txt"
  fi

  boolstr=`echo "${jtemp} <= 10" | bc -l`

  if [ "${boolstr}" = "1" ]; then
    commandstr="${commandstr} mergedEMu2M=datacard_MEMu2M_Y1.txt"

    tmpstr=`echo "${jtemp} >= 5" | bc -l`

    if [ "${tmpstr}" = "1" ]; then
      commandstr="${commandstr} mergedEle3E=datacard_ME3E_Y1.txt"
    fi

    tmpstr=`echo "${jtemp} < 5" | bc -l`

    if [ "${tmpstr}" = "1" ]; then
      commandstr="${commandstr} mergedEle3E=datacard_ME3E_CR_Y1.txt"
    fi
  fi

  boolstr=`echo "${jtemp} <= 2" | bc -l`

  if [ "${boolstr}" = "1" ]; then
    commandstr="${commandstr} mergedMu3M=datacard_MMFF_Y1.txt"
    commandstr="${commandstr} mergedMu2E=datacard_MM2E_CR_Y1.txt"
  fi

  boolstr=`echo "${jtemp} >= 5" | bc -l`

  if [ "${boolstr}" = "1" ]; then
    commandstr="${commandstr} resolvedEle=datacard_REFF_Y1.txt"
  fi

  boolstr=`echo "${jtemp} >= 2" | bc -l`

  if [ "${boolstr}" = "1" ]; then
    commandstr="${commandstr} resolvedEMu=datacard_REMuFF_Y1.txt"
  fi

  boolstr=`echo "${jtemp} >= 1" | bc -l`

  if [ "${boolstr}" = "1" ]; then
    commandstr="${commandstr} resolvedMu=datacard_RMFF_Y1.txt"

    tmpstr=`echo "${jtemp} < 5" | bc -l`

    if [ "${tmpstr}" = "1" ]; then
      commandstr="${commandstr} resolvedEle=datacard_REFF_CR_Y1.txt"
    fi

    tmpstr=`echo "${jtemp} < 2" | bc -l`

    if [ "${tmpstr}" = "1" ]; then
      commandstr="${commandstr} resolvedEMu=datacard_REMuFF_CR_Y1.txt"
    fi
  fi

  commandstr="${commandstr} &> A${j}/combined.txt"
  echo ${commandstr}
  eval "${commandstr}"
  echo "s/A1/A${j}/g"
  sed -i "s/A1/A${j}/g" A${j}/combined.txt

  mkdir Z${j}
  commandstr="combineCards.py"

  boolstr=`echo "${jtemp} <= 10" | bc -l`

  if [ "${boolstr}" = "1" ]; then
    commandstr="${commandstr} mergedEle3E=datacard_ME3E_Y1.txt"
    commandstr="${commandstr} mergedEMu2M=datacard_MEMu2M_Y1.txt"
    commandstr="${commandstr} mergedEle2E=datacard_ME2E_CR_Y1.txt"
  fi

  boolstr=`echo "${jtemp} <= 2" | bc -l`

  if [ "${boolstr}" = "1" ]; then
    commandstr="${commandstr} mergedMu2E=datacard_MM2E_Y1.txt"
    commandstr="${commandstr} mergedMu3M=datacard_MMFF_Y1.txt"
    commandstr="${commandstr} mergedEMu1M=datacard_MEMu1M_CR_Y1.txt"
  fi

  boolstr=`echo "${jtemp} >= 5" | bc -l`

  if [ "${boolstr}" = "1" ]; then
    commandstr="${commandstr} resolvedEle=datacard_REFF_Y1.txt"
  fi

  boolstr=`echo "${jtemp} >= 1" | bc -l`

  if [ "${boolstr}" = "1" ]; then
    commandstr="${commandstr} resolvedEMu=datacard_REMuFF_Y1.txt"
  fi

  boolstr=`echo "${jtemp} >= 1" | bc -l`

  if [ "${boolstr}" = "1" ]; then
    commandstr="${commandstr} resolvedMu=datacard_RMFF_Y1.txt"

    tmpstr=`echo "${jtemp} < 5" | bc -l`

    if [ "${tmpstr}" = "1" ]; then
      commandstr="${commandstr} resolvedEle=datacard_REFF_CR_Y1.txt"
    fi
  fi

  commandstr="${commandstr} &> Z${j}/combined.txt"
  echo ${commandstr}
  eval "${commandstr}"
  echo "s/A1/Z${j}/g"
  sed -i "s/A1/Z${j}/g" Z${j}/combined.txt

done

for j in 250 500 750 1000 1500
do
  if [ ${j} -le 750 ]; then
    mkdir A${j}
    commandstr="combineCards.py resolvedEle=datacard_REFF_Y1.txt resolvedEMu=datacard_REMuFF_Y1.txt resolvedMu=datacard_RMFF_Y1.txt &> A${j}/combined.txt"
    echo ${commandstr}
    eval "${commandstr}"
    echo "s/A1/A${j}/g"
    sed -i "s/A1/A${j}/g" A${j}/combined.txt
  fi

  mkdir Z${j}
  commandstr="combineCards.py resolvedEle=datacard_REFF_Y1.txt resolvedEMu=datacard_REMuFF_Y1.txt resolvedMu=datacard_RMFF_Y1.txt &> Z${j}/combined.txt"
  echo ${commandstr}
  eval "${commandstr}"
  echo "s/A1/Z${j}/g"
  sed -i "s/A1/Z${j}/g" Z${j}/combined.txt
done 

