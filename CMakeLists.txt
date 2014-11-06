# =====================================================================================
#
#       Filename:  CMakeLists.txt
#
#    Description:  Extract Topologically distinct surfaces form data set
#
#        Version:  1.0
#        Created:  12/14/2013 04:35:23 PM
#       Revision:  none
#       Compiler:  gcc
#
#         Author:  Siavash Ameli
#   Organization:  University of California, Berkeley
#
# =====================================================================================

cmake_minimum_required(VERSION 2.8)
project(FeatureSurface CXX)

# ======================
# VTK/ParaView Inclusion
# ======================

# ParaView plugin option
set(BUILD_PARAVIEW_PLUGIN ON CACHE BOOL "Build ParaView plugin?")

# Avoid confliction of Paraview and VTK inclusion
if(BUILD_PARAVIEW_PLUGIN)
    # include ParaView instead of VTK
    find_package(ParaView REQUIRED)
    if(ParaView_FOUND)
        include(${PARAVIEW_USE_FILE})
    else(ParaView_FOUND)
        message(FATAL_ERROR "ParaView not found.")
    endif(ParaView_FOUND)
else(BUILD_PARAVIEW_PLUGIN)
    # include VTK
    find_package(VTK REQUIRED)
    if(VTK_FOUND)
        include(${VTK_USE_FILE})
    else(VTK_FOUND)
        message(FATAL_ERROR "VTK not found.")
    endif(VTK_FOUND)
endif(BUILD_PARAVIEW_PLUGIN)

# ==============
# Compiler Flags
# ==============

# Initialize CXX FLAGS
set(CMAKE_CXX_FLAGS                "-Wall -std=c++0x ${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG          "-O0 -g")
set(CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE        "-O4 -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")

# Compiler-specific C++11 activation
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
    if (NOT (GCC_VERSION VERSION_GREATER 4.4 OR GCC_VERSION VERSION_EQUAL 4.4))
        message(FATAL_ERROR "${PROJECT_NAME} Requires g++ 4.4 or greater.")
    endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libstdc++")
# else()
    # message(FATAL_ERROR "Your C++ compiler does not support C++11.")
endif()

# ================
# Source Inclusion
# ================

# Headers
set(PROJECT_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
include_directories(${PROJECT_INCLUDE_DIR})

# Main Source
set(PROJECT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

set(MAIN_EXECUTABLE_FILENAME "FeatureSurface")
add_executable(${MAIN_EXECUTABLE_FILENAME} 
    ${PROJECT_SOURCE_DIR}/main.cxx
    ${PROJECT_SOURCE_DIR}/FeatureSurface.cxx)
target_link_libraries(${MAIN_EXECUTABLE_FILENAME} ${VTK_LIBRARIES})

# Test Source
# set(TEST_EXECUTABLE_FILENAME "TestFeatureSurface")
# add_executable(${TEST_EXECUTABLE_FILENAME}
#     ${PROJECT_SOURCE_DIR}/TestFeatureSurface.cxx
#     ${PROJECT_SOURCE_DIR}/FeatureSurface.cxx)
# target_link_libraries(${TEST_EXECUTABLE_FILENAME} ${VTK_LIBRARIES})

# ==================
# Output Directories
# ==================

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin CACHE PATH "Directory for all executables.")
set(LIBRARY_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin CACHE PATH "Directory for all libraries.")

# ===============
# ParaView Plugin
# ===============

if(BUILD_PARAVIEW_PLUGIN)
    ADD_PARAVIEW_PLUGIN(
        FeatureSurfacePlugin
        "1.0"
        DOCUMENTATION_DIR       ${CMAKE_CURRENT_SOURCE_DIR}/doc
        SERVER_MANAGER_XML      ${PROJECT_SOURCE_DIR}/FeatureSurface.xml
        GUI_RESOURCE_FILES      ${PROJECT_SOURCE_DIR}/FeatureSurfaceGUI.xml
        SERVER_MANAGER_SOURCES  ${PROJECT_SOURCE_DIR}/FeatureSurface.cxx)
endif(BUILD_PARAVIEW_PLUGIN)

# ============
# Installation
# ============

option(INSTALL_PROJECT "Install project on host system" ON)

if(${INSTALL_PROJECT})
    # set(PROJECT_INSTALL_PATH ${CMAKE_CURRENT_SOURCE_DIR}/install)
    set(PROJECT_INSTALL_PATH /usr/local/bin)
    install(TARGETS ${MAIN_EXECUTABLE_FILENAME} DESTINATION ${PROJECT_INSTALL_PATH})
endif(${INSTALL_PROJECT})

# ================
# Generate Package
# ================

option(GENERATE_PACKAGE "Generate system dependent package" ON)
if(${GENERATE_PACKAGE})
    include(InstallRequiredSystemLibraries)

    # Custom variables
    set(AUTHOR "Siavash Ameli")
    set(PACKAGE_NAME "FeatureSurface")
    set(PACKAGE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bin)

    # cpack variables
    set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/License.txt")
    set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/ReadMe.md")
    set(CPACK_PACKAGE_NAME ${PACKAGE_NAME})
    set(CPACK_DEBIAN_PACKAGE_MAINTAINER ${AUTHOR})
    set(CPACK_PACKAGE_CONTACT "sameli@berkeley.edu")
    set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "i386")

    # Package directory
    set(CPACK_PACKAGE_FILE_NAME "${PACKAGE_DIR}/${PACKAGE_NAME}")

    # Package path
    # set(CPACK_TOPLEVEL_TAG ${PACKAGE_DIR})
    # set(CPACK_PACKAGING_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/bin)
    # set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/bin)
    # set(CPACK_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/bin)
    # set(CPACK_PACKAGE_INSTALL_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)
    # set(CPACK_PACKAGING_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/bin)
    
    # Package install path
    if(NOT DEFINED PROJECT_INSTALL_PATH)
        set(CPACK_SET_DESTDIR true)
        set(CPACK_INSTALL_PREFIX ${CMAKE_CURRENT_PROJECT}/bin)
    endif(NOT DEFINED PROJECT_INSTALL_PATH)

    # Binary option
    option(BINARY_PACKAGE "Generate binary package" ON)

    # OS
    if(CMAKE_HOST_SYSTEM_NAME MATCHES "Linux")
        if(${BINARY_PACKAGE})
            set(CPACK_GENERATOR "DEB")
        else(${BINARY_PACKAGE})
            set(CPACK_GENERATOR "SOURCE_TGZ")
        endif(${BINARY_PACKAGE})
    elseif()
        # TODO: other OS
    endif()

    include(CPack)
endif(${GENERATE_PACKAGE})
