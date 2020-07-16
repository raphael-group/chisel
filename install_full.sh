#!/bin/bash

set -e
set -x

# Finding whether os is Linux or MacOSX
OS=$(uname -s)
case ${OS} in
    Linux*)     OS=Linux;;
    Darwin*)    OS=MacOSX;;
    *)          echo "Unknown OS ${OS}; please use manual installation." && exit 1;;
esac

# Finding whether machine is 32bit or 64bit
case ${OS} in
    Linux)
	VER=$(arch)
	case ${VER} in
	    x86_64)    VER="x86_64";;
	    *)         VER="x86";;
	esac;;
    MacOSX)
	VER=$(arch)
	case ${VER} in
	    x86_64)    VER="x86_64";;
	    *)         echo "Mac OSX 32 bit is not supported by Miniconda, please use custom installation." && exit 1;;
	esac;;
    *)
        echo "Unknown OS ${OS}; please use manual installation." && exit 1;;
esac

# Installing Miniconda
CHISEL_HOME=$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )
cd ${CHISEL_HOME}
curl -L https://repo.anaconda.com/miniconda/Miniconda2-latest-${OS}-${VER}.sh > miniconda.sh
rm -rf ./conda/
bash miniconda.sh -b -f -p ./conda/
export CONDA_HOME=${CHISEL_HOME}/conda/bin/

# Installing chisel
${CONDA_HOME}/conda config --add channels defaults
${CONDA_HOME}/conda config --add channels bioconda
${CONDA_HOME}/conda config --add channels conda-forge
${CONDA_HOME}/conda create -n chisel chisel -y

# Activating CHISEL
source ${CONDA_HOME}/activate chisel
echo -e "\nInstallation was succesfull and CHISEL is ready!\nPlease remember to run the following command now and during every new session before using CHISEL:\n\t\tsource ${CONDA_HOME}/activate chisel\n\n"
