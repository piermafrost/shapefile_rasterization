# This makefile includes one optional argument variable: MODE, defining the compilation model
# By default, MODE = standard, but two other modes can be specified while calling make: "debug" and "optim"
#     make MODE=standard    (standard options)
#     make MODE=optim       (for fast use. DEFAULT MODE)
#     make MODE=debug       (extra debugging options)


#################################################
# Main customization variables (default values) #
#################################################

# Fortran compiler
FC ?= gfortran

# Compilation mode (default):
MODE = optim

# Path of netCDF library (must contains subrepertories lib/ and include/)
NCPATH ?= /usr


#############################
# Compiler-dependent flags: #
#############################


ifeq ($(findstring pgfortran, $(FC)), pgfortran)

#=========================#
# PGI compiler: pgfortran #
#=========================#

# Fortran language option
lang_flags = -Mfreeform

# Error and warning options:
ifeq ($(MODE), optim)
warn_flags = -Minform=severe
else ifeq ($(MODE), debug)
warn_flags = -Minform=inform
else # assume standard mode
warn_flags = -Minform=warn
endif

# Debugging options:
ifeq ($(MODE), debug)
debug_flags = -g
endif

# Code generation options:
ifneq ($(MODE), optim)
code_flags = -Mbounds
endif
ifeq ($(MODE), debug)
code_flags += -traceback
endif

# Optimization options:
ifeq ($(MODE), debug)
optim_flags = -O0
else ifeq ($(MODE), optim)
optim_flags = -O4 -fast
else
# if none of previous: assume standard mode:
optim_flags = -O2
endif



else ifeq ($(findstring gfortran, $(FC)), gfortran)

#========================#
# GNU compiler: gfortran #
#========================#

# Fortran language options:
lang_flags = -ffree-form -std=gnu

# Error and warning options:
ifeq ($(MODE), debug)
warn_flags = -Wall -Wconversion-extra -Wimplicit-interface -Wunderflow -Wextra -Wunreachable-code
else
warn_flags = -Waliasing -Wampersand -Wline-truncation -Wcharacter-truncation -Wconversion -Wunderflow -Wunreachable-code
ifeq ($(MODE), standard)
warn_flags += -Wpedantic -Wunused-dummy-argument
endif
endif

# Debugging options:
ifneq ($(MODE), optim)
debug_flags = -ffpe-trap=invalid,zero,overflow
endif

# Code generation options:
ifneq ($(MODE), optim)
code_flags = -fbounds-check
endif
ifeq ($(MODE), debug)
code_flags += -fbacktrace -g3
endif

# Optimization options:
ifeq ($(MODE), debug)
optim_flags = -O0 -fstack-protector-all
else ifeq ($(MODE), optim)
optim_flags = -O3
else
# if none of previous: assume standard mode:
optim_flags = -O1
endif



else ifeq ($(findstring ifort, $(FC)), ifort)

#========================
# Intel compiler: ifort #
#========================

# Fortran langauge options:
lang_flags = -free -132

# Error and warning options:
ifeq ($(MODE), debug)
warn_flags = -warn all
else
warn_flags = -warn general,usage,declaration,truncated_source,ignore_loc
endif

# Debugging options:
ifeq ($(MODE), debug)
debug_flags = -debug full
endif

# Options at run time
ifeq ($(MODE), debug)
code_flags = -check all -fp-stack-check -traceback
else ifeq ($(MODE), standard)
code_flags = -check bounds
endif

# Optimization options:
ifeq ($(MODE), debug)
optim_flags = -O0 -fstack-protector-all -fstack-security-check
else ifeq ($(MODE), optim)
optim_flags = -Ofast
else
# if none of previous: assume standard mode:
optim_flags = -O1
endif



endif



########################
# LIBRARY and INCLUDE: #
########################


# NetCDF library:
inc_flags += -I$(NCPATH)/include
lib_flags += -L$(NCPATH)/lib
NETCDF_FLAGS ?= -lnetcdf -lnetcdff


##############################################
# Complete list of fortran compilation flags #
##############################################

main_flags = $(lang_flags) $(warn_flags) $(debug_flags) $(code_flags) $(optim_flags)
FFLAGS = $(main_flags) $(inc_flags) $(lib_flags) $(NETCDF_FLAGS)


#############
# Root path #
#############

# Check if path assignment uses correct path
# If not, add 'updatepath' to the executable dependencies

# get root path (remove 'source' from current directory path)
pathcommand = character(len=*), parameter:: root_path = "$(CURDIR:source=)"
# get path assignment command (last line of 'path.inc' file)
currpathcommand = $(shell tail -1 path.inc)

ifneq ($(pathcommand), $(currpathcommand))
pathsource = updatepath path.inc
else
pathsource = path.inc
endif


#############################################################################


# Specific executable making rule:
shp2grid: shp2grid.o shp2grid_functions.o shapefile_io_module.o netcdf_output.o geography.o miscellaneous_functions.o scan_arguments.o
	$(FC) $^ -o $@ $(FFLAGS)

test_read_shapefile: test_read_shapefile.o shapefile_io_module.o miscellaneous_functions.o
	$(FC) $^ -o $@ $(FFLAGS)


# Redefine implicit rule to for making objects:
%.o: %.f90
	$(FC) $< -o $@ -c $(FFLAGS)


# Additional dependencies
test_read_shapefile.o: shapefile_io_module.o miscellaneous_functions.o
netcdf_output.o: geography.o miscellaneous_functions.o
shp2grid_functions.o: netcdf_output.o geography.o
shp2grid.o: shp2grid_functions.o shapefile_io_module.o netcdf_output.o geography.o miscellaneous_functions.o scan_arguments.o $(pathsource)


############################
# Update root path command #
############################

# create path.inc file if it doesn't exist
path.inc:
	touch $@

.PHONY: updatepath # make sure to always execute that command when called
updatepath:
	echo '! File automatically generated by `make`' > path.inc
	echo '$(pathcommand)' >> path.inc


#############################################################################


check:
	@echo "$(FC) $(FFLAGS)"

.PHONY: clean clc

clean:
	rm -f *.o *.mod
clc:
	rm -f *.o *.mod shp2grid test_read_shapefile

