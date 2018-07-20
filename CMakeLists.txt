# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /apps\\RNAMotif/CMakeLists.txt
#
# CMakeLists.txt file for RNAMotif.
# ===========================================================================

cmake_minimum_required (VERSION 3.0.0)
project (seqan_apps_RNAMotif CXX)
message (STATUS "Configuring apps/RNAMotif")

ENABLE_LANGUAGE(C)

set(CMAKE_MODULE_PATH $ENV{HOME}/Entwicklung/seqan/util/cmake)
set(SEQAN_INCLUDE_PATH $ENV{HOME}/Entwicklung/seqan/include)

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
set (SEQAN_FIND_DEPENDENCIES ZLIB BZip2 OpenMP)
find_package (SeqAn REQUIRED)
# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS_MAIN})
include_directories ("/usr/include/ViennaRNA")
include_directories (SYSTEM ${SEQAN_INCLUDE_DIRS_DEPS})
include_directories ("folding_utils/ipknot")

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})
add_definitions (-DWITH_GLPK -DHAVE_CONFIG_H)

set (IPPATH folding_utils/ipknot)
set (COPATH ${IPPATH}/contrafold)
set (NUPATH ${IPPATH}/nupack)

set (CONTRA_SOURCES ${COPATH}/Config.hpp ${COPATH}/SStruct.hpp ${COPATH}/SStruct.cpp ${COPATH}/InferenceEngine.hpp ${COPATH}/Utilities.hpp ${COPATH}/Utilities.cpp)
set (NUPACK_SOURCES ${NUPATH}/nupack.h ${NUPATH}/nupack.cpp ${NUPATH}/def_param.h ${NUPATH}/dptable.h)
set (IPKNOT_SOURCES ${IPPATH}/config.h ${IPPATH}/aln.h ${IPPATH}/aln.cpp ${IPPATH}/fold.h ${IPPATH}/fold.cpp ${IPPATH}/ip.h ${IPPATH}/ip.cpp ${CONTRA_SOURCES} ${NUPACK_SOURCES})

# Update the list of file names below if you add source files to your application.
add_executable (RNAMotif RNAMotif.cpp motif.h motif_structures.h motif_search.h stockholm_file.h stockholm_io.h folding_utils/RNAlib_utils.h folding_utils/IPknot_utils.h ${IPKNOT_SOURCES} stored_interval_tree.h)

# Add dependencies found by find_package (SeqAn).
#target_link_libraries (RNAMotif ${SEQAN_LIBRARIES} "/usr/lib/x86_64-linux-gnu/libRNA.a" glpk gmp)
target_link_libraries (RNAMotif ${SEQAN_LIBRARIES} "/usr/local/lib/libRNA.a" glpk gmp)

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# App Test
# ----------------------------------------------------------------------------

#seqan_add_app_test (RNAMotif)

# ----------------------------------------------------------------------------
# CPack Install
# ----------------------------------------------------------------------------

if (SEQAN_BUILD_SYSTEM STREQUAL "APP:RNAMotif")
  set (CPACK_PACKAGE_NAME "RNAMotif")
  set (CPACK_PACKAGE_DESCRIPTION_SUMMARY "RNAMotif")
  set (CPACK_DEBIAN_PACKAGE_MAINTAINER "Benjamin Strauch")
  set (CPACK_PACKAGE_VENDOR "Benjamin Strauch")

  seqan_configure_cpack_app (RNAMotif "RNAMotif")
endif (SEQAN_BUILD_SYSTEM STREQUAL "APP:RNAMotif")

