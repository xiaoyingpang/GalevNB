#!make

###############################################
#                                             #
#  Makefile for GalevNB                       #
#                                             #
#  by C. Olczak                               #
#                                             #
#  last modification: 2013-03-06              #
#                                             #
###############################################


#------------------------
# Set global parameters.
#------------------------

.KEEP_STATE:

# Set exectuable.
RESULT = GalevNBsevs

# Set path of source files.
PATH_GALEV=$(PWD)


#----------------------------------------------------------------------------
# Standard compilation.
#----------------------------------------------------------------------------
standard:
	$(MAKE) $(RESULT) -j4 \
	"FC    = gfortran" \
	"CPP   = g++" \
	"DEBUG = no"
#----------------------------------------------------------------------------
# Compile for maximum debugging.
#----------------------------------------------------------------------------
debug:
	$(MAKE) $(RESULT) -j4 \
	"FC    = gfortran" \
	"CPP   = g++" \
	"DEBUG = yes"
#----------------------------------------------------------------------------
# Clean directory: keep only source files.
#----------------------------------------------------------------------------
clean:
	rm -f \
	*.o \
	*~ \
	"#"*"#" \
	*.lst \
	*.mod \
	*__genmod.* \
	core \
	loadmap \
	$(RESULT)
#----------------------------------------------------------------------------


#------------------------------
# Set standard compiler flags.
#------------------------------
#ifneq ($(DEBUG), yes)
  OPTIMIZE  = -O3
#else
#  OPTIMIZE  = -O0
#endif
CPP_FLAGS = $(OPTIMIZE) -DPATH_GALEV="$(PATH_GALEV)" $(DEBUG_FLAGS)
FC_FLAGS  = $(OPTIMIZE) --free-line-length-0 $(DEBUG_FLAGS)
LD_FLAGS  = -lm -lgsl -lgslcblas -lstdc++



#--------------------------
# Check optional settings.
#--------------------------
ifeq ($(DEBUG), yes)
  DEBUG_FLAGS += -Wall -Wextra -ggdb3 -fbounds-check -ftrapv -fexceptions -fstack-check -fbacktrace
endif



#-------------------
# Set object lists.
#-------------------
OBJECTS = \
	filter.o \
	filterlib.o \
	spec.o \
	speclib.o \
	spectra.o \
	spectra_asciiprecalib.o \
	spectra_blackbody.o \
	spectra_blank.o \
	spectra_lejeune.o \
	GalevNB.o \

INCLUDES = $(wildcard *.h)



#---------------------
# Set compiler rules.
#---------------------

# Set suffixes.
.SUFFIXES: .f .f90 .F .F90 .o .mod .cpp

# Set compiler rules.
.f.o:
	$(FC)  $(FC_FLAGS)  -c $<

.f90.o:
	$(FC)  $(FC_FLAGS)  -c $<

.f90.mod:
	$(FC)  $(FC_FLAGS)  -c $<

.F.o:
	$(FC)  $(FC_FLAGS)  -c $<

.F90.o:
	$(FC)  $(FC_FLAGS)  -c $<

.F90.mod:
	$(FC)  $(FC_FLAGS)  -c $<

.cpp.o:
	$(CPP) $(CPP_FLAGS) -c $<



# Set object dependencies.
$(OBJECTS): $(INCLUDES)



# Set dependencies.
$(RESULT): $(OBJECTS)
	$(FC) $(FC_FLAGS) $(OBJECTS) -o $@ $(LD_FLAGS)
