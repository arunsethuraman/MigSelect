

###
###  minimal makefile for  meanpartition   
###  
###   Makefile will build all .c and .cpp files in all subdirs and
###    the current directory into a single executable binary file.
###   All .o files are created in the $(BLD) directory.
###
###   Files are built in following order:
###     $(SUBDIRS)/*.c	---> $(BLD)/*.o
###     *.c		---> $(BLD)/*.o
###     $(SUBDIRS)/*.cpp	---> $(BLD)/*.o
###     *.cpp		---> $(BLD)/*.o
###
###  Files listed in $(EXTRA_DIST) are excluded from the build.
###
###  Including another .c or .cpp file to the resulting executable binary 
###  does not require any changes to this file.  To add a another 
###  subdirectory, edit $(SUBDIRS).
###  
###


PROGRAM=meanpartition

SUBDIRS=

EXTRA_DIST = msdirent.c msdirent.h 

SUFFIX=

BLDDIR =BLD

HEADERS := $(foreach DIR,$(SUBDIRS), -I $(DIR)/.)

ALL= mkbld_dir $(PROGRAM) 

RM=rm -f
#########################################################
#
# setup any paths to various tools on the build machine.
#
#########################################################

#SHAREDIR=../share

#########################################################
#  Compile 
#########################################################

CC = gcc
#CC=g++

#
#   debugging flags   --  turn on debugging here 
#
#DEBUG=-g3 -save-temps -D "DEBUG" -print-search-dirs
#DEBUG=-g3 -save-temps -D "DEBUG" -D "DEBUGL"
#DEBUG=-g3 -D "DEBUG" -D "DEBUGL" -D "DEBUG1" -D "DEBUG_ALL" 
#DEBUG=-g3 -D "DEBUG"

#
#   C flags
#
#   note: compile flags set by existing automake creation process 
# 	CFLAGS = -g -O2	 
# 	CPP = gcc -E      (set but not used)
#
#   The default, if no C language dialect options are given, is `-std=gnu89'
CFLAGS = -g -O2 


#
#  any extra flags needed set here
#
#EXTRAFLAG= -Wno-deprecated 
#EXTRAFLAG= -Wall -ansi -pedantic -Wno-deprecated 
#EXTRAFLAG= -Wall -Wextra 
EXTRAFLAG=-D"NOCONFIG"

#EXTRAFLAG=-D"SANITY_TEST" -D"NOCONFIG"


#
#  options for code coverage tools  needs -g also
#
#PROFILEOPTS=-static -fprofile-arcs -ftest-coverage 

#
#   C++ flags
#
# The default, if no C++ language dialect options are given, is `-std=gnu++98'
#CXXFLAGS= -g -O2 -std=c++0x
CXXFLAGS= -g -O2 


#
#
#

ALL_CFLAGS=	\
	$(DEBUG) \
	$(CFLAGS) \
	$(EXTRAFLAG) \
	$(PROFILEOPTS)


ALL_CXXFLAGS=	\
	$(DEBUG) \
	$(CXXFLAGS) \
	$(EXTRAFLAG)	\
	$(PROFILEOPTS)


#########################################################
#  Link Flags
#########################################################

LINKOPTS=

CCLD=gcc
#CCLD=g++

LINK_ALL= \
	$(LINKOPTS) \
	$(PROFILEOPTS)

LIBS= \
	-lm 

# create full paths for all extra files that get excluded from binary
SUBDIR_EXTRA_DIST := $(foreach DIR, $(SUBDIRS), $(EXTRA_DIST:%=$(DIR)/%)) 

# create list of all C objects in sub-directories, filter out the extra files
SUBDIR_C_OBJS := $(filter-out $(SUBDIR_EXTRA_DIST), \
		$(foreach DIR, $(SUBDIRS), $(wildcard $(DIR)/*.c)))

# create list of all C files in current directory, filter out the extra files
LOCAL_C_OBJS := $(filter-out $(EXTRA_DIST), $(wildcard *.c))

# combine the lists of C files and change suffix to .o
C_OBJS := $(patsubst %.c,%.o, $(SUBDIR_C_OBJS) $(LOCAL_C_OBJS))

# create list of all C++ objects in sub-directories, filter out the extra files
SUBDIR_CXX_OBJS := $(filter-out $(SUBDIR_EXTRA_DIST), \
		$(foreach DIR, $(SUBDIRS), $(wildcard $(DIR)/*.cpp)))

# create list of all C++ files in current directory, filter out the extra files
LOCAL_CXX_OBJS := $(filter-out $(EXTRA_DIST), $(wildcard *.cpp))

# combine the lists C++ files and change suffix to .o
CXX_OBJS := $(patsubst %.cpp,%.o, $(SUBDIR_CXX_OBJS) $(LOCAL_CXX_OBJS))

# combine the lists of C and C++ files 
ALLOBJS := $(C_OBJS) $(CXX_OBJS)


all:  $(ALL)
	@echo -e "--------------------  target: < $@ > complete -----\n"

.c.o: 
	@echo -e "--------------------  .c to .o target: < $@ > -----"
	$(CC) $(ALL_CFLAGS) -I. $(HEADERS) -c  $^ -o $(BLDDIR)/$(notdir $@)
	@echo -e "\nC target: < $@ > complete -----\n"

.cpp.o: 
	@echo -e "--------------------  .cpp to .o target: < $@ > -----"
	$(CXX) $(ALL_CXXFLAGS) -I. $(HEADERS) -c  $^ -o $(BLDDIR)/$(notdir $@)
	@echo -e "\nC++ target: < $@ > complete -----\n"

$(PROGRAM): $(ALLOBJS) 
	@echo -e "--------------------------------------"
	@echo -e "\nobject build done ...."    
	@echo -e "--------------------------------------"
	@echo -e "linking   < $@ >  ....  \n"
	
#        libtool --tag=CC --mode=link 
	cd $(BLDDIR); $(CCLD) $(CFLAGS) $(LINK_ALL) -o ../$@$(SUFFIX) $(notdir $(ALLOBJS)) $(LIBS) 
	@echo -e "\n--------------------  target: < $@ > complete -----\n"

mkbld_dir:
	@if test ! -d $(BLDDIR); then mkdir $(BLDDIR); fi
	@echo 

install: 
	@echo
	@echo current target: --$@-- 
	@echo

clean:
	@echo
	@echo Cleaning up.
	${RM} $(BLDDIR)/*
	${RM} ${PROGRAM} 
	@echo

cleanCoverage:
	@echo
	@echo Cleaning up.
	${RM} *.gcov
