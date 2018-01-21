# val_extract, sstval_extract, argpar, phash, shash, pqueue, vincenty

## ToDo
### Fix
  * Any stack arrays over ~4KB, switch to malloc (done?)
  * Clean up or refactor free()ing arrays before returning (remove copy-pasta)
  * val_extract.c: Change convert_lat_lons to do things in smaller chunks

### Add
  * CMake paradigm for data files
  * Unit tests for val_extract and sstval_extract (and *-mains?)
  * argpar-help: start using flags for which sections to print
  * argpar: handle quoted keys and values?
  * argpar: option to prefix options with --? (separate library, maybe?)
  * argpar-help: run-time modification of columns
  * vincenty and pqueue versioning

## val_extract 2.5.9 - 2016-09-12
### Fixed
  * val_extract.c: fixed utime calculation for times occuring in local DST

## sstval_extract - 2016-09-12
### Removed
  * Removed all sstval_extract content

## val_extract 2.5.8 - 2016-09-12
### Added
  * val_extract: added platform and instrument to the base output file, units
    to each variable

## val_extract 2.5.7 - 2016-07-29
### Added
  * val_extract: extra documentation for clustering algorithms

## val_extract 2.5.6 - 2016-05-24
### Fixed
  * val_extract: bug involving gap detection and bad geometry values

## val_extract 2.5.5, sstval_extract 1.0.2 - 2016-04-22
### Fixed
  * val_extract: bug preventing products from being selected

## val_extract 2.5.4, sstval_extract 1.0.2 - 2015-10-28
### Fixed
  * val_extract: changed error codes, 1-99 non-fatal, >= 100 fatal
  * sstval_extract: changed error codes, 1-109 non-fatal, >= 110 fatal
  
## val_extract 2.5.3 - 2015-10-26
### Fixed
  * val_extract.c: fixed bug in reading flag_masks
  
## val_extract-main 1.2.1 - 2015-10-19
### Added
  * versions printed to output file, added version=1 argument

## val_extract 2.5.2 - 2015-10-14
### Fixed
  * val_extract.c: fixed older style flag readings
  * val_extract.c, val_extract.h: changed flags to uint32_t

## val_extract 2.5.1 - 2015-10-05
### Fixed
  * val_extract.c: IQR no longer done with less than eight pixels
  * val_extract.c: fixed a bug finding the median in even-length arrays
  
## val_extract-main 1.2.0 - 2015-10-01
### Added
  * Updated file output to reflect various val_extract changes
  
## val_extract 2.5.0 - 2015-10-01
### Added
  * val_extract.c, val_extract.h: count_flags argument
  
## vincenty 1.0.1 - 2015-10-01
### Fixed
  * vincenty.h: Documentation fixes
  
## val_extract 2.4.0 - 2015-10-01
### Fixed
  * val_extract.c: moved some stack arrays to heap 
  * all: renamed some structure variables, chiefly related to pixel counts 
  
### Added
  * val_extract.c: filtered and IQR stats

## vincenty 1.0.0 - 2015-09-23
### Added
  * vincenty.c: library to calculate distance between two coordinate pairs

## val_extract 2.3.0 - 2015-09-23
### Fixed
  * val_extract.c: added support for older L2 files

## val_extract 2.2.2 - 2015-08-13
### Fixed
  * val_extract.c: for lat/lon box extracts: fixed logic for checking whether a point is in the
   granule, including dateline detection; added data gap and failed geometry detection

## pqueue 1.0.0 - 2015-08-12
### Added
  * pqueue.h: priority queue interface
  * pqueue-ll.c: linked list based priority queue
  * tests/pqueue/: unit tests
 
## val_extract 2.2.1 - 2015-06-14
### Fixed
  * val_extract.c: fixed center-of-region detections

## argpar 2.0.2 - 2015-06-08
### Fixed
  * argpar.c, argpar-help.c, argpar.h: silenced -Wextra warnings and changed some types

## val_extract 2.2.0 - 2015-06-08
### Added
  * val_extract.c: changed L2QC style and added some filtering features

## argpar 2.0.1 - 2015-06-03
### Fixed
  * argpar.c: fixed compiler warning by adding best-practice return check

## val_extract-main 1.1.1 - 2015-06-02
### Fixed
  * val_extract-main.c: added a newline at the end of the base output files

## val_extract-main 1.1.0 - 2015-05-28
### Added
  * val_extract-main.c: added RMS to variable output files

## val_extract 2.1.0 - 2015-05-28
### Fixed
  * val_extract.h: buffer overflow error when C optimization flags are used

### Added
  * val_extract.h, val_extract.c: functions to retrieve version numbers

## 1.0.0 - 2015-03-12
### Added
  * Changelog.md

