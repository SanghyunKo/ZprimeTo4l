#!/bin/sh
set -e

pwddir=$PWD

for j in 0p4 0p6 0p8 1 1p5 2 5 10 50 100
do

  mkdir A${j}_mu
  commandstr="combineCards.py"
  jtemp=`echo ${j} | tr p .`

  boolstr=`echo "${jtemp} <= 2" | bc -l`

  if [ "${boolstr}" = "1" ]; then
    commandstr="${commandstr} mergedEMu1M=datacard_MEMu1M_CR_Y1.txt"
    commandstr="${commandstr} mergedMu2E=datacard_MM2E_CR_Y1.txt"
  fi

  boolstr=`echo "${jtemp} <= 2" | bc -l`

  if [ "${boolstr}" = "1" ]; then
    commandstr="${commandstr} mergedMu3M=supplemental_datacards/datacard_MMFF_YY_Y1.txt"
  fi

  boolstr=`echo "${jtemp} >= 1" | bc -l`

  if [ "${boolstr}" = "1" ]; then
    commandstr="${commandstr} resolvedMu=supplemental_datacards/datacard_RMFF_YY_Y1.txt"
    commandstr="${commandstr} resolvedEMu=datacard_REMuFF_CR_Y1.txt"
    commandstr="${commandstr} resolvedEle=datacard_REFF_CR_Y1.txt"
  fi

  commandstr="${commandstr} &> A${j}_mu/combined.txt"
  echo ${commandstr}
  eval "${commandstr}"
  echo "s/A1/A${j}/g"
  sed -i "s/A1/A${j}/g" A${j}_mu/combined.txt

  mkdir Z${j}_mu
  commandstr="combineCards.py"
  jtemp=`echo ${j} | tr p .`

  boolstr=`echo "${jtemp} <= 2" | bc -l`

  if [ "${boolstr}" = "1" ]; then
    commandstr="${commandstr} mergedMu2E=supplemental_datacards/datacard_MM2E_ZY_Y1.txt"
    commandstr="${commandstr} mergedMu3M=supplemental_datacards/datacard_MMFF_ZY_Y1.txt"
    commandstr="${commandstr} mergedEMu1M=datacard_MEMu1M_CR_Y1.txt"
  fi

  boolstr=`echo "${jtemp} >= 1" | bc -l`

  if [ "${boolstr}" = "1" ]; then
    commandstr="${commandstr} resolvedEle=datacard_REFF_CR_Y1.txt"
    commandstr="${commandstr} resolvedMu=supplemental_datacards/datacard_RMFF_ZY_Y1.txt"
    commandstr="${commandstr} resolvedEMu=supplemental_datacards/datacard_REMuFF_ZY_mu_Y1.txt"
  fi

  commandstr="${commandstr} &> Z${j}_mu/combined.txt"
  echo ${commandstr}
  eval "${commandstr}"
  echo "s/A1/Z${j}/g"
  sed -i "s/A1/Z${j}/g" Z${j}_mu/combined.txt

done

for j in 250 500 750 1000 1500
do

  if [ ${j} -le 750 ]; then
    mkdir A${j}_mu
    commandstr="combineCards.py resolvedEle=datacard_REFF_CR_Y1.txt resolvedEMu=datacard_REMuFF_CR_Y1.txt resolvedMu=supplemental_datacards/datacard_RMFF_YY_Y1.txt &> A${j}_mu/combined.txt"
    echo ${commandstr}
    eval "${commandstr}"
    echo "s/A1/A${j}/g"
    sed -i "s/A1/A${j}/g" A${j}_mu/combined.txt
  fi

  mkdir Z${j}_mu
  commandstr="combineCards.py resolvedEle=datacard_REFF_CR_Y1.txt resolvedEMu=supplemental_datacards/datacard_REMuFF_ZY_mu_Y1.txt resolvedMu=supplemental_datacards/datacard_RMFF_ZY_Y1.txt &> Z${j}_mu/combined.txt"
  echo ${commandstr}
  eval "${commandstr}"
  echo "s/A1/Z${j}/g"
  sed -i "s/A1/Z${j}/g" Z${j}_mu/combined.txt

done 

