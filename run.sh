#!/bin/bash
source /opt/intel/oneapi/setvars.sh > /dev/null 2>&1
/bin/echo "##" $(whoami) is compiling Cylinder Flow Simulation -- Shq Project - 1 of 1 cylinder-flow-simulation
export OverrideDefaultFP64Settings=1 
export IGC_EnableDPEmulation=1
rm ./a.out
icpx -fsycl src/main.cpp
if [ $? -eq 0 ]; then ./a.out > output ; fi
