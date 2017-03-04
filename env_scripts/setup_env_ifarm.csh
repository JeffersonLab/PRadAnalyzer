#!/bin/csh
# set environment on JLab ifarm based on the version
# Linux ifarm1402.jlab.org 3.10.0-327.el7.x86_64
# 1 SMP Thu Nov 19 22:10:57 UTC 2015 x86_64 x86_64 x86_64 GNU/Linux

setenv PRAD_PATH /work/hallb/prad/PRadAnalyzer

# setup the coda library and caen library
# these are needed for online monitoring and high voltage control
if (`uname -m` == 'x86_64') then
    setenv ET_LIB /u/site/coda/3.06/Linux/lib64
    setenv ET_INC /u/site/coda/3.06/common/include
    setenv THIRD_LIB  ${PRAD_PATH}/thirdparty/lib64
else
    setenv ET_LIB /u/site/coda/2.6.2/Linux/lib
    setenv ET_INC /u/site/coda/2.6.2/common/include
    setenv THIRD_LIB  ${PRAD_PATH}/thirdparty/lib
endif

# comment this if you have your own root installed
# setup root environment
source /u/apps/root/6.08.00/setroot_CUE.csh

# setup the libraries path
setenv PRAD_LIB ${PRAD_PATH}/libs

# setup LD_LIBRARY_PATH
if ! $?LD_LIBRARY_PATH then
	setenv LD_LIBRARY_PATH ${THIRD_LIB}:${ET_LIB}:${PRAD_LIB}
else
	setenv LD_LIBRARY_PATH ${THIRD_LIB}:${ET_LIB}:${PRAD_LIB}:${LD_LIBRARY_PATH}
endif
