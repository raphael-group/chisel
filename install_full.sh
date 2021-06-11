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
	VER=$(uname -i)
	case ${VER} in
	    x86_64)    MINICONDA="https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh";;
	    *)         MINICONDA="https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86.sh";;
	esac;;
    MacOSX)
	VER=$(uname -m)
	case ${VER} in
	    *)         MINICONDA="https://repo.anaconda.com/miniconda/Miniconda2-latest-MacOSX-x86_64.sh";;
	esac;;
    *)
        echo "Unknown OS ${OS}; please use manual installation." && exit 1;;
esac

# Installing Miniconda
CHISEL_HOME=$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )
cd ${CHISEL_HOME}
curl -L ${MINICONDA} > miniconda.sh
rm -rf ./conda/
bash miniconda.sh -b -f -p ./conda/
export CONDA_HOME=${CHISEL_HOME}/conda/bin

# Installing chisel
${CONDA_HOME}/conda config --add channels defaults
${CONDA_HOME}/conda config --add channels bioconda
${CONDA_HOME}/conda config --add channels conda-forge
${CONDA_HOME}/conda create -n chisel chisel -y

# Activating CHISEL
source ${CONDA_HOME}/activate chisel
echo -e "\nInstallation was succesfull and CHISEL is ready!\nPlease remember to run the following command now and during every new session before using CHISEL:\n\n\t\tsource ${CONDA_HOME}/activate chisel\n\n"
