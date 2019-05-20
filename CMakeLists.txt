# Created by the script cgal_create_cmake_script_with_options
# This is the CMake script for compiling a set of CGAL applications.

project( offset_surface_mesh )


cmake_minimum_required(VERSION 2.8.11)



# CGAL and its components
find_package( CGAL QUIET COMPONENTS Core )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()

# include helper file
include( ${CGAL_USE_FILE} )


# Boost and its components
find_package( Boost REQUIRED)

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()
  
add_definitions( -DUSE_FORWARD_DECL)

# include for local package
include_directories( BEFORE ../../include )
include_directories (BEFORE ./include)

# Creating entries for all .cpp/.C files with "main" routine
# ##########################################################

include( CGAL_CreateSingleSourceCGALProgram )

create_single_source_cgal_program( "main.cpp" "kdprocessor.cpp" "kdtree.c")

set_target_properties("main" PROPERTIES OUTPUT_NAME "off_set_mesh")
