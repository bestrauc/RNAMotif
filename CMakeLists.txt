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

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search SeqAn and select dependencies.
find_package (SeqAn REQUIRED)
find_package (OpenMP REQUIRED)
# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS_MAIN})
include_directories ("/usr/include/ViennaRNA")
include_directories (SYSTEM ${SEQAN_INCLUDE_DIRS_DEPS})
include_directories ("ipknot")

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

#set(IPKNOT_SOURCES ipknot/aln.h ipknot/aln.cpp ipknot/boltzmann_param.h ipknot/boltzmann_param.h ipknot/bpseq.h ipknot_bpseq.cpp ipknot/config.h ipknot/fa.h ipknot/fa.cpp ipknot/fold.h ipknot/fold.cpp ipknot.i
file (GLOB IPKNOT_SOURCES ipknot/*.h ipknot/*.cpp ipknot/*.ipp ipknot/*.c)
file (GLOB CONTRA_SOURCES ipknot/contrafold/*.hpp ipknot/contrafold/*.cpp ipknot/contrafold/*.ipp)
file (GLOB NUPACK_SOURCES ipknot/nupack/*.h ipknot/nupack/*.cpp)

# Update the list of file names below if you add source files to your application.
add_executable (RNAMotif RNAMotif.cpp motif.h motif_structures.h RNAlib_utils.h IPknot_utils.h ${IPKNOT_SOURCES} ${CONTRA_SOURCES} ${NUPACK_SOURCES})

# Add dependencies found by find_package (SeqAn).
target_link_libraries (RNAMotif ${SEQAN_LIBRARIES} "/usr/lib/x86_64-linux-gnu/libRNA.a")

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Installation
# ----------------------------------------------------------------------------

# Set variables for installing, depending on the selected build type.
if (NOT SEQAN_PREFIX_SHARE_DOC)
  seqan_setup_install_vars (RNAMotif)
endif (NOT SEQAN_PREFIX_SHARE_DOC)

# Install RNAMotif in ${PREFIX}/bin directory
install (TARGETS RNAMotif
         DESTINATION bin)

# Install non-binary files for the package to "." for app builds and
# ${PREFIX}/share/doc/RNAMotif for SeqAn release builds.
install (FILES LICENSE
               README
         DESTINATION ${SEQAN_PREFIX_SHARE_DOC})
#install (FILES example/fasta1.fa
#               example/fasta2.fa
#         DESTINATION ${SEQAN_PREFIX_SHARE_DOC}/example)

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

