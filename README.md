# shapefile\_rasterization
An efficient tool to rasterize shapefile polygons as  "polygon area on grid cells"

## Description of repository

### Fortran source files
There are two Fortran main programs in current repository:
* shp2grid.f90:             the key program, to rasterize shapefiles
* test\_read\_shapefile.f90:  a simple script to test the shapefile reading functions

The main programs use the following modules:
* shapefile\_io\_module.f90:     central module, to read shapefile data (\*.shp and \*.shx files)
* netcdf\_output.f90:           contains subroutine to generate and write netCDF outputs
* scan\_arguments.f90:          subroutine for analyzing command-line arguments passed to the main program
* geography.f90:               geographic tools: Earth model, radius, and authalic latitude
* miscellaneous\_functions.f90: various tools: modify file names, time string, etc.

### Other subrepertories
* backups/:      backup of old or incomplete code
* data\_example/: shapefiles for testing

## Required libraries
The main program (shp2grid) need the netcdf-fortran library (https://docs.unidata.ucar.edu/netcdf-fortran/current/index.html).  
On Linux operating systems: `apt install libnetcdff`.

Because I couldn't find a Fortran shapefile library, I wrote my own module to read shapefile data from binary files \*.shx
and \*.shp (with Fortan option "form='unformatted'"): shapefile\_io\_module.f90

## Compilation
Before compiling, make sure that the netCDF-Fortran library is installed and configured for your Fortran compiler.
Then, update the environement variables `export FC=_your_fortran_compiler_` and `export NCPATH=_current_netcdf_path_`.
Tip: use `nc-config` to determine those elements.  
Alternatively, open the Makefile with any text editor and directly edit the Makefile variables "FC", "NCPATH",
"inc\_flags", "lib\_flags" or "NETCDF\_FLAGS".

Compile the executable with `make` (shortcut for `make shp2grid`), or `make test_read_shapefile`.
Note that this second program does not need the netCDF library.

To compile with debugging options, do `make clc` and `make shp2grid MODE=debug`.

## General description of the rasterization method
### Purpose
...

### Computation steps
...

### KNOWN ISSUES
...

## User guide of `shp2grid` command
A Comprehensive "help" text is written in 'help\_message.txt' (and displayed when typing `./shp2grid --help`)

## References
Shapefiles example in subrepertory data\_example/ are taken from database:  
Global Self-consistent Hierarchical High-resolution Geography, GSHHG  
Version 2.3.7 June 15, 2017  
Distributed under the Lesser GNU Public License  
Credit: Paul Wessel (primary contact: pwessel@hawaii.edu), Walter H. F. Smith
Url: http://www.soest.hawaii.edu/pwessel/gshhg/

