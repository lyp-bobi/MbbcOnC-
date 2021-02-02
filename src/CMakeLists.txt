###############################################################################
#
# src/CMakeLists.txt controls building of SIDX library
#
# Copyright (c) 2009 Mateusz Loskot <mateusz@loskot.net>
#
###############################################################################

###############################################################################
# Source files specification
#
# Naming format:
#   SIDX_CPP - all the .cpp files
#   SIDX_HPP - all the .hpp files
#   SIDX_<dirname>_CPP - all the .cpp files for the given subdir/namespace
#   ...

set(SIDX_HEADERS_DIR "${PROJECT_SOURCE_DIR}/include/spatialindex/")
set(SIDX_HEADERS_CAPI_DIR "${PROJECT_SOURCE_DIR}/include/spatialindex/capi")
set(SIDX_HPP "")
set(SIDX_CPP "")

set(SIDX_SRC_DIR "${PROJECT_SOURCE_DIR}/src")
set(SIDX_CAPI_DIR "${PROJECT_SOURCE_DIR}/src/capi")

#
# base
#
set(SIDX_BASE_HPP
  "${SIDX_HEADERS_DIR}/LineSegment.h"
  "${SIDX_HEADERS_DIR}/MovingPoint.h"
  "${SIDX_HEADERS_DIR}/Point.h"
  "${SIDX_HEADERS_DIR}/Region.h"
  "${SIDX_HEADERS_DIR}/MBCRTree.h"
  "${SIDX_HEADERS_DIR}/SpatialIndex.h"
  "${SIDX_HEADERS_DIR}/TimePoint.h"
  "${SIDX_HEADERS_DIR}/STPoint.h"
  "${SIDX_HEADERS_DIR}/TimeRegion.h"
  "${SIDX_HEADERS_DIR}/Version.h"
  "${SIDX_HEADERS_DIR}/../storagemanager/TrajStore.h"
  "${SIDX_HEADERS_DIR}/ShapeList.h"
)
list (APPEND SIDX_HPP ${SIDX_BASE_HPP} )

set(SIDX_CAPI_HPP
  "${SIDX_HEADERS_CAPI_DIR}/BoundsQuery.h"
  "${SIDX_HEADERS_CAPI_DIR}/CountVisitor.h"
  "${SIDX_HEADERS_CAPI_DIR}/CustomStorage.h"
  "${SIDX_HEADERS_CAPI_DIR}/Error.h"
  "${SIDX_HEADERS_CAPI_DIR}/IdVisitor.h"
  "${SIDX_HEADERS_CAPI_DIR}/LeafQuery.h"
  "${SIDX_HEADERS_CAPI_DIR}/ObjVisitor.h"
#  "${SIDX_HEADERS_CAPI_DIR}/sidx_api.h"
  "${SIDX_HEADERS_CAPI_DIR}/sidx_config.h"
  "${SIDX_HEADERS_CAPI_DIR}/sidx_impl.h"
  "${SIDX_HEADERS_CAPI_DIR}/Utility.h"
)
list (APPEND SIDX_HPP ${SIDX_CAPI_HPP} )

set(SIDX_CAPI_CPP
  "${SIDX_CAPI_DIR}/BoundsQuery.cc"
  "${SIDX_CAPI_DIR}/CountVisitor.cc"
  "${SIDX_CAPI_DIR}/CustomStorage.cc"
  "${SIDX_CAPI_DIR}/Error.cc"
  "${SIDX_CAPI_DIR}/IdVisitor.cc"
  "${SIDX_CAPI_DIR}/LeafQuery.cc"
  "${SIDX_CAPI_DIR}/ObjVisitor.cc"
#  "${SIDX_CAPI_DIR}/sidx_api.cc"
  "${SIDX_CAPI_DIR}/Utility.cc")


set(SIDX_SPATIALINDEX_CPP
  "${SIDX_SRC_DIR}/spatialindex/LineSegment.cc"
  "${SIDX_SRC_DIR}/spatialindex/MovingPoint.cc"
  "${SIDX_SRC_DIR}/spatialindex/MovingRegion.cc"
  "${SIDX_SRC_DIR}/spatialindex/MBC.cc"
  "${SIDX_SRC_DIR}/spatialindex/Cylinder.cc"
  "${SIDX_SRC_DIR}/spatialindex/Trajectory.cc"
  "${SIDX_SRC_DIR}/spatialindex/Point.cc"
  "${SIDX_SRC_DIR}/spatialindex/Region.cc"
  "${SIDX_SRC_DIR}/spatialindex/SpatialIndexImpl.cc"
  "${SIDX_SRC_DIR}/spatialindex/TimePoint.cc"
  "${SIDX_SRC_DIR}/spatialindex/STPoint.cc"
  "${SIDX_SRC_DIR}/spatialindex/TimeRegion.cc"
  "${SIDX_SRC_DIR}/spatialindex/ShapeList.cc"
    "${SIDX_SRC_DIR}/spatialindex/xMBC.cc"
    "${SIDX_SRC_DIR}/spatialindex/xCylinder.cc"
    "${SIDX_SRC_DIR}/spatialindex/xTrajectory.cc"
    "${SIDX_SRC_DIR}/spatialindex/xPoint.cc"
        "${SIDX_SRC_DIR}/spatialindex/xMBR.cc"
        "${SIDX_SRC_DIR}/spatialindex/xLine.cc"
        "${SIDX_SRC_DIR}/spatialindex/xSBB.cc"
)
list (APPEND SIDX_CPP ${SIDX_SPATIALINDEX_CPP} )



set(SIDX_MBCRTREE_CPP
        "${SIDX_SRC_DIR}/mbcrtree/BulkLoader.cc"
        "${SIDX_SRC_DIR}/mbcrtree/BulkLoader.h"
        "${SIDX_SRC_DIR}/mbcrtree/Index.cc"
        "${SIDX_SRC_DIR}/mbcrtree/Index.h"
        "${SIDX_SRC_DIR}/mbcrtree/Leaf.cc"
        "${SIDX_SRC_DIR}/mbcrtree/Leaf.h"
        "${SIDX_SRC_DIR}/mbcrtree/Node.cc"
        "${SIDX_SRC_DIR}/mbcrtree/Node.h"
        "${SIDX_SRC_DIR}/mbcrtree/PointerPoolNode.h"
        "${SIDX_SRC_DIR}/mbcrtree/MBCRTree.cc"
        "${SIDX_SRC_DIR}/mbcrtree/MBCRTree.h"
        "${SIDX_SRC_DIR}/mbcrtree/Statistics.cc"
        "${SIDX_SRC_DIR}/mbcrtree/Statistics.h"
        "${SIDX_SRC_DIR}/mbcrtree/TrajStoreRelated.cc"
        )
list (APPEND SIDX_CPP ${SIDX_MBCRTREE_CPP})


set(SIDX_XRTREE_CPP
        "${SIDX_SRC_DIR}/xrtree/BulkLoader.cc"
        "${SIDX_SRC_DIR}/xrtree/BulkLoader.h"
        "${SIDX_SRC_DIR}/xrtree/Index.cc"
        "${SIDX_SRC_DIR}/xrtree/Index.h"
        "${SIDX_SRC_DIR}/xrtree/Leaf.cc"
        "${SIDX_SRC_DIR}/xrtree/Leaf.h"
        "${SIDX_SRC_DIR}/xrtree/Node.cc"
        "${SIDX_SRC_DIR}/xrtree/Node.h"
        "${SIDX_SRC_DIR}/xrtree/PointerPoolNode.h"
        "${SIDX_SRC_DIR}/xrtree/xRTree.cc"
        "${SIDX_SRC_DIR}/xrtree/xRTree.h"
        "${SIDX_SRC_DIR}/xrtree/Statistics.cc"
        "${SIDX_SRC_DIR}/xrtree/Statistics.h"
        "${SIDX_SRC_DIR}/xrtree/IndexBuildHelper.cc"
        "${SIDX_SRC_DIR}/xrtree/SBBForest.cc"
        )
list (APPEND SIDX_CPP ${SIDX_XRTREE_CPP})


set(SIDX_STORAGEMANAGER_CPP
  "${SIDX_SRC_DIR}/storagemanager/Buffer.cc"
  "${SIDX_SRC_DIR}/storagemanager/DiskStorageManager.cc"
  "${SIDX_SRC_DIR}/storagemanager/DiskStorageManager.h"
  "${SIDX_SRC_DIR}/storagemanager/MemoryStorageManager.cc"
  "${SIDX_SRC_DIR}/storagemanager/MemoryStorageManager.h"
  "${SIDX_SRC_DIR}/storagemanager/RandomEvictionsBuffer.cc"
  "${SIDX_SRC_DIR}/storagemanager/RandomEvictionsBuffer.h"
  "${SIDX_SRC_DIR}/storagemanager/TrajStore.cc"
  "${SIDX_SRC_DIR}/storagemanager/xStore.cc"
)
list (APPEND SIDX_CPP ${SIDX_STORAGEMANAGER_CPP})

set(SIDX_RAND48 )

if (NOT HAVE_SRAND48)
    set(SIDX_RAND48 "${SIDX_SRC_DIR}/tools/rand48.cc")
endif()

set(SIDX_TOOLS_CPP
  ${SIDX_RAND48}
  "${SIDX_SRC_DIR}/tools/Tools.cc"
)
list (APPEND SIDX_CPP ${SIDX_TOOLS_CPP})

#
# Group source files for IDE source explorers (e.g. Visual Studio)
#

source_group("CMake Files" FILES CMakeLists.txt)

source_group("Header Files" FILES ${SIDX_BASE_HPP})
source_group("C API Header Files" FILES ${SIDX_CAPI_HPP})

source_group("Source Files" FILES ${SIDX_BASE_CPP})
source_group("C API Source Files" FILES ${SIDX_CAPI_CPP})


# Standard include directory of SIDX library
include_directories(../include)

set (APPS_CPP_DEPENDENCIES "${SIDX_LIB_NAME}" CACHE INTERNAL "libraries to link")

###############################################################################
# Targets settings

set(SIDX_SOURCES
  ${SIDX_HPP}
  ${SIDX_CPP})

# NOTE:
# This hack is required to correctly link static into shared library.
# Such practice is not recommended as not portable, instead each library,
# static and shared should be built from sources separately.
#if(UNIX)
#  add_definitions("-fPIC")
#endif()

if(WIN32)
    add_definitions("-DSIDX_DLL_EXPORT=1")
if (NOT WITH_STATIC_LASZIP)
    add_definitions("-DSIDX_DLL_IMPORT=1")
endif()
endif()

add_library(${SIDX_LIB_NAME} SHARED ${SIDX_SOURCES})

add_library(${SIDX_C_LIB_NAME} SHARED ${SIDX_CAPI_CPP})

target_link_libraries(${SIDX_C_LIB_NAME}
  ${SIDX_LIB_NAME}
)



set_target_properties(${SIDX_LIB_NAME}
    PROPERTIES VERSION "${SIDX_LIB_VERSION}"
               SOVERSION "${SIDX_LIB_SOVERSION}"  )

set_target_properties(${SIDX_C_LIB_NAME}
    PROPERTIES VERSION "${SIDX_LIB_VERSION}"
               SOVERSION "${SIDX_LIB_SOVERSION}" )



if (APPLE)
  set_target_properties(
    ${SIDX_LIB_NAME}
    PROPERTIES
    INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib" BUILD_WITH_INSTALL_RPATH OFF)

  set_target_properties(
    ${SIDX_C_LIB_NAME}
    PROPERTIES
    INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib" BUILD_WITH_INSTALL_RPATH OFF)
endif()

###############################################################################
# Targets installation

install(TARGETS ${SIDX_LIB_NAME} ${SIDX_C_LIB_NAME}
  RUNTIME DESTINATION ${SIDX_BIN_DIR}
  LIBRARY DESTINATION ${SIDX_LIB_DIR}
  ARCHIVE DESTINATION ${SIDX_LIB_DIR})

install(DIRECTORY ${SIDX_HEADERS_DIR}
  DESTINATION include/spatialindex
  FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")
