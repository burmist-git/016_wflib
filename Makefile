PYTHONSETUPFILE=setup.py
PKGNAME := $(shell grep 'name=' ${PYTHONSETUPFILE} | sed 's/,//'| sed 's/name=//' | sed 's/"//g' | sed 's/ //g')
PKGVERSION := $(shell grep 'version=' ${PYTHONSETUPFILE} | sed 's/,//'| sed 's/version=//' | sed 's/"//g' | sed 's/ //g')
CONDAENVNAME=comporttest
CONDAENVYML=${CONDAENVNAME}_env.yml

.PHONY: environmentsetup environmentsetupyml exportpkglist rmenv installpkg uninstallpkg printinfo clean printhelp

##############################################################
ROOTCFLAGS  = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS    = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS   = $(shell $(ROOTSYS)/bin/root-config --glibs)
ROOTLDFLAGS = $(shell $(ROOTSYS)/bin/root-config --ldflags)
# add xml files
ROOTLIBS += -lMLP -lXMLIO

OUTLIB = ./obj/
SRC = ./src/

PATHTOSHAREDLIB=$(OUTLIB)

CXX  = g++
CXX += -I./     

CXXFLAGS  = -g -Wall -fPIC -Wno-deprecated
CXXFLAGS += $(ROOTCFLAGS)
CXXFLAGS += $(ROOTLIBS)
CXXFLAGS += $(ROOTGLIBS)
CXXFLAGS += $(ROOTLDFLAGS)
CXXFLAGSDICT = -g -Wall -fPIC -Wno-deprecated

#----------------------------------------------------#

#all: makedir convertFebTemp_main convertMergerTemp_main runanaGL840_main

## Makes obj and other support directory 
makedir:
	mkdir -p $(OUTLIB);
	mkdir -p tmp;

## prints cxx info
printCXXinfo:
	$(info CXX          = "$(CXX)")
	$(info CXXFLAGS     = "$(CXXFLAGS)")
	$(info CXXFLAGSDICT = "$(CXXFLAGSDICT)")

## make makefile helper/reminder
printmakehelp_and_reminder: setup.py wflib.h
	$(info  /**********************************************************************/)
	$(info  * task --> printmakehelp_and_reminder: setup.py wflib.h *)
	$(info  * $$@ ----> $@                                 *)
	$(info  * $$< --------------------------------> $<                      *)
	$(info  * $$^ --------------------------------> $^ *)
	$(info  /**********************************************************************/)

## build waveform object file
waveform.o: $(SRC)/waveform.cpp $(SRC)/waveform.hh
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)/$@ $<

## build waveform object file 
wfSim.o: $(SRC)/wfSim.cpp  $(SRC)/wfSim.hh $(SRC)/wfSimConst.hh
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)/$@ $<

## creating wflib - shared library
wflib.so: obj/waveform.o obj/wfSim.o
	$(CXX) -shared -o $(PATHTOSHAREDLIB)$@ $^ $(ROOTLIBS) $(ROOTGLIBS)

#libarichstandaloneDict.cxx: libarichstandalone.h LinkDef.h
#	rootcling -f $@ -c $(CXXFLAGSDICT) -p $^

#libarichstandaloneDictLib.so: libarichstandaloneDict.cxx ARICHChannelHist.cc testARICHChannelHist.cc ARICHSectorHist.cc testARICHSectorHist.cc ARICHmappingCopperMergerFeb.cc convertFebTemp.cc testARICHmapping.cc
#	g++ -shared -o $@ $(CXXFLAGS) -I. -I$(ROOTSYS)/include $^

## Set up python interpreter environment via conda
environmentsetup:
	./MakefileHelper.sh --condaenv ${CONDAENVNAME}

## Set up python interpreter environment
environmentsetupyml:
	conda env create -f ${CONDAENVYML}

## Export packages list to yml file
exportpkglist:
	./MakefileHelper.sh --exportenv ${CONDAENVYML}

## remove conda environment
rmenv:
	./MakefileHelper.sh --condaenvrm ${CONDAENVNAME}

## pip install package
installpkg:
	./MakefileHelper.sh --pipinstall ${PKGNAME} ${PKGVERSION}

## pip uninstall package
uninstallpkg:
	pip uninstall -y ${PKGNAME}
	@echo conda list | grep -i ${PKGNAME} | wc -l



## Print info test
printinfo:
	@echo MAKEFILE_LIST = ${MAKEFILE_LIST}
	@echo PYTHONSETUPFILE = ${PYTHONSETUPFILE}
	@echo PKGNAME = ${PKGNAME}
	@echo PKGVERSION = ${PKGVERSION}
	@echo CONDAENVNAME = ${CONDAENVNAME}
	@echo CONDAENVYML = ${CONDAENVYML}

## Clean the folders from unused files and others (created by python) 
clean:
	rm -f *~
	rm -f .*~
	rm -f */*~
	rm -f */.*~
	rm -rf ${PKGNAME}.egg-info
	rm -rf .ipynb_checkpoints/
	rm -rf ./${PKGNAME}/.ipynb_checkpoints
	rm -rf ./${PKGNAME}/__pycache__
	rm -f *~
	rm -f .*~
	rm -rf $(OUTLIB)/*

########################################################################
# Self Documenting Commands                                            #
# Copied with modification from :                                      #
# https://github.com/hgrif/example-project/blob/master/Makefile        #
########################################################################

.DEFAULT_GOAL := printhelp

## Self-Documented Makefile : print descriptions of the described commands
printhelp:
	@echo "$$(tput bold)Available rules:$$(tput sgr0)"
	@echo
	@sed -n -e "/^## / { \
		h; \
		s/.*//; \
		:doc" \
		-e "H; \
		n; \
		s/^## //; \
		t doc" \
		-e "s/:.*//; \
		G; \
		s/\\n## /---/; \
		s/\\n/ /g; \
		p; \
	}" ${MAKEFILE_LIST} \
	| LC_ALL='C' sort --ignore-case \
	| awk -F '---' \
		-v ncol=$$(tput cols) \
		-v indent=19 \
		-v col_on="$$(tput setaf 6)" \
		-v col_off="$$(tput sgr0)" \
	'{ \
		printf "%s%*s%s ", col_on, -indent, $$1, col_off; \
		n = split($$2, words, " "); \
		line_length = ncol - indent; \
		for (i = 1; i <= n; i++) { \
			line_length -= length(words[i]) + 1; \
			if (line_length <= 0) { \
				line_length = ncol - indent - length(words[i]) - 1; \
				printf "\n%*s ", -indent, " "; \
			} \
			printf "%s ", words[i]; \
		} \
		printf "\n"; \
	}'
