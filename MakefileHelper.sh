#!/bin/bash

function export_yml_env_sh {
    yml_file_name=$1
    conda env export | grep -v "^prefix: " > $yml_file_name
}

function create_conda_env_sh {
    condaEvnName=$1
    wcl=$(conda env list | grep $condaEvnName | wc -l)
    if [ $wcl = 0 ]; then
	conda create --name $condaEvnName -y numpy pyserial notebook 
	conda env list | grep -i $condaEvnName
	echo " "
	echo "> conda activate $condaEvnName"
    else
        conda env list | grep -i $condaEvnName
        echo "> conda activate $condaEvnName"
    fi
}

function rm_conda_env_sh {
    condaEvnName=$1
    condaEvnCurrentName=$(conda env list | grep '*' | awk {'print $1'})
    if [ $condaEvnCurrentName = $condaEvnName ]; then
	echo " "
	echo "> conda deactivate"
    else
	conda env remove -n $condaEvnName
    fi
}

function pipinstall_sh {
    pkgName=$1
    pkgVersion=$2
    wcl=$(echo $pkgVersion | grep -i x | wc -l)
    if [ $wcl = 0 ]; then
	echo "pkgName = $pkgName"
	echo "pkgVersion = $pkgVersion"
	pip install .
    else
	echo "pkgName = $pkgName"
	echo "pkgVersion = $pkgVersion"
	pip install -e .
    fi
    conda list | grep -i $pkgName 
}

# Libraries to be added
function install_packages_sh {
    condaEvnName=$1
    wcl=$(conda env list | grep -i $condaEvnName | wc -l)
    if [ $wcl = 0 ]; then
	echo "install_packages_sh"
	# pandas-profiling
	#conda search -c conda-forge pandas-profiling=2.9.0
	#echo "installing --> pandas-profiling=2.9.0"
	#conda install -y -c conda-forge pandas-profiling=2.9.0
	# scikit-learn
	#conda search scikit-learn=0.23.2
	#echo "installing --> scikit-learn==0.23.2"
	#conda install -y scikit-learn=0.23.2
	# plotly
	#conda search plotly=4.13.0
	#echo "installing --> plotly=4.13.0"
	#conda install -y plotly=4.13.0    
    else
        conda env list | grep -i $condaEvnName
        echo "> conda activate $condaEvnName"
    fi
}

function nm_sh {
    liblist=obj/wflib.so 
    for filename in $liblist ; do
        echo " "
        echo "----> $filename"
        nm -CD --defined-only $filename | awk '{if($2=="T"){{print }}}' | grep -v _fini | grep -v _init | nl
        echo " "
    done
}

function printHelp {
    echo " --> ERROR in input arguments"
    echo " -h           : print help"
    echo " --nm         : list shared library functions and classes"
    echo " --condaenv   : create conda environment"
    echo "        [1]   : environment name"
    echo " --condaenvrm : remove conda environment"
    echo "        [1]   : environment name"
    echo " --exportenv  : export env into yml file"
    echo "        [1]   : yml file name"
    echo " --pipinstall : pip install python package"
    echo "        [1]   : package name"
    echo "        [2]   : package version"
    echo " --install    : install conda packages"
    echo "        [1]   : environment name"
}

if [ $# -eq 0 ]; then
    printHelp
else
    if [ "$1" = "-h" ]; then
	printHelp
    elif [ "$1" = "--condaenv" ]; then
	if [ $# -eq 2 ]; then
	    create_conda_env_sh $2
	else
	    printHelp
	fi
    elif [ "$1" = "--condaenvrm" ]; then
	if [ $# -eq 2 ]; then
	    rm_conda_env_sh $2
	else
	    printHelp
	fi
    elif [ "$1" = "--exportenv" ]; then
	if [ $# -eq 2 ]; then
	    export_yml_env_sh $2
	else
	    printHelp
	fi
    elif [ "$1" = "--pipinstall" ]; then
	if [ $# -eq 3 ]; then
	    pipinstall_sh $2 $3
	else
	    printHelp
	fi
    elif [ "$1" = "--install" ]; then
	if [ $# -eq 2 ]; then
	    install_packages_sh $2
	else
	    printHelp
	fi
    elif [ "$1" = "--nm" ]; then
	nm_sh
    else
        printHelp
    fi
fi
