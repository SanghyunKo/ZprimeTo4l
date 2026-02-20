#!/bin/sh
set -e

pwddir=$PWD

for j in 0p4 0p6 0p8 1 1p5 2 5 10 50 100
do

  mkdir A${j}_el
  commandstr="combineCards.py"
  jtemp=`echo ${j} | tr p .`

  boolstr=`echo "${jtemp} <= 10" | bc -l`

  if [ "${boolstr}" = "1" ]; then
    commandstr="${commandstr} mergedEle2E=supplemental_datacards/datacard_ME2E_YY_Y1.txt"
    commandstr="${commandstr} mergedEMu2M=datacard_MEMu2M_CR_Y1.txt"
  fi

  boolstr=`echo "${jtemp} <= 10" | bc -l`

  if [ "${boolstr}" = "1" ]; then
    tmpstr=`echo "${jtemp} >= 5" | bc -l`

    if [ "${tmpstr}" = "1" ]; then
      commandstr="${commandstr} mergedEle3E=supplemental_datacards/datacard_ME3E_YY_Y1.txt"
    fi

    tmpstr=`echo "${jtemp} < 5" | bc -l`

    if [ "${tmpstr}" = "1" ]; then
      commandstr="${commandstr} mergedEle3E=datacard_ME3E_CR_Y1.txt"
    fi
  fi

  boolstr=`echo "${jtemp} >= 5" | bc -l`

  if [ "${boolstr}" = "1" ]; then
    commandstr="${commandstr} resolvedEle=supplemental_datacards/datacard_REFF_YY_Y1.txt"
    commandstr="${commandstr} resolvedEMu=datacard_REMuFF_CR_Y1.txt"
    commandstr="${commandstr} resolvedMu=datacard_RMFF_CR_Y1.txt"
  fi

  commandstr="${commandstr} &> A${j}_el/combined.txt"
  echo ${commandstr}
  eval "${commandstr}"
  echo "s/A1/A${j}/g"
  sed -i "s/A1/A${j}/g" A${j}_el/combined.txt

  mkdir Z${j}_el
  commandstr="combineCards.py"
  jtemp=`echo ${j} | tr p .`

  boolstr=`echo "${jtemp} <= 10" | bc -l`

  if [ "${boolstr}" = "1" ]; then
    commandstr="${commandstr} mergedEle3E=supplemental_datacards/datacard_ME3E_ZY_Y1.txt"
    commandstr="${commandstr} mergedEMu2M=supplemental_datacards/datacard_MEMu2M_ZY_Y1.txt"
    commandstr="${commandstr} mergedEle2E=datacard_ME2E_CR_Y1.txt"
  fi

  boolstr=`echo "${jtemp} >= 5" | bc -l`

  if [ "${boolstr}" = "1" ]; then
    commandstr="${commandstr} resolvedEle=supplemental_datacards/datacard_REFF_ZY_Y1.txt"
  fi

  boolstr=`echo "${jtemp} >= 1" | bc -l`

  if [ "${boolstr}" = "1" ]; then
    commandstr="${commandstr} resolvedEMu=supplemental_datacards/datacard_REMuFF_ZY_el_Y1.txt"
    commandstr="${commandstr} resolvedMu=datacard_RMFF_CR_Y1.txt"

    tmpstr=`echo "${jtemp} < 5" | bc -l`

    if [ "${tmpstr}" = "1" ]; then
      commandstr="${commandstr} resolvedEle=datacard_REFF_CR_Y1.txt"
    fi
  fi

  commandstr="${commandstr} &> Z${j}_el/combined.txt"
  echo ${commandstr}
  eval "${commandstr}"
  echo "s/A1/Z${j}/g"
  sed -i "s/A1/Z${j}/g" Z${j}_el/combined.txt

done

for j in 250 500 750 1000 1500
do

  if [ ${j} -le 750 ]; then
    mkdir A${j}_el
    commandstr="combineCards.py resolvedEle=supplemental_datacards/datacard_REFF_YY_Y1.txt resolvedEMu=datacard_REMuFF_CR_Y1.txt resolvedMu=datacard_RMFF_CR_Y1.txt &> A${j}_el/combined.txt"
    echo ${commandstr}
    eval "${commandstr}"
    echo "s/A1/A${j}/g"
    sed -i "s/A1/A${j}/g" A${j}_el/combined.txt
  fi

  mkdir Z${j}_el
  commandstr="combineCards.py resolvedEle=supplemental_datacards/datacard_REFF_ZY_Y1.txt resolvedEMu=supplemental_datacards/datacard_REMuFF_ZY_el_Y1.txt resolvedMu=datacard_RMFF_CR_Y1.txt &> Z${j}_el/combined.txt"
  echo ${commandstr}
  eval "${commandstr}"
  echo "s/A1/Z${j}/g"
  sed -i "s/A1/Z${j}/g" Z${j}_el/combined.txt

done 

