
# NCCMP

## 2.1.0 - 2015-03-18
Added function to recursively compare subgroups within groups
Also automatically set -m when -g opt specified on command line.

### Source Changed
    * nccmp.cpp
    * nccmp.hpp
    * opt.c
    
## 2.0.1 - 2015-03-12
### Added
  * Changelog.md

## 2.0.0 - 2015-03-16

Converted program to C++ in order to convert macros to templates for easier debugging.
Added nanequal to bring code closer to public version
Added options:
    
    -C <n> Output only the first n error messages
    -n     Do not use default tolerance of 0.00001% (-T) for floats and doubles
    -N     Nans are equal 

### Source Changed
  deleted code
    * nccmp.c
    * nccmp.h
  
  added code
    * nccmp.cpp
    * nccmp.hpp
    
  modified code
    * opt.c
    * opt.h