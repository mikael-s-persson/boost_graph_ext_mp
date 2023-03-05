
cmake_minimum_required(VERSION 3.5)

project(bgl_extras)

set(CMAKE_CXX_STANDARD 17)

enable_testing()

include(CheckCXXCompilerFlag)

# Needed to set the namespace for both the export targets and the
# alias libraries
set(cmake_package_name bgl_extras CACHE INTERNAL "")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -stdlib=libc++ -fPIC -O3 -Wall -Woverloaded-virtual -Wold-style-cast -Wnon-virtual-dtor")

# Process googletest's CMake build system
add_subdirectory(third_party/googletest)

# Look for Boost:

set(Boost_DEBUG TRUE)

set(Boost_ADDITIONAL_VERSIONS "1.70" "1.70.0" "1.71" "1.71.0" "1.72" "1.72.0" "1.73" "1.73.0" "1.74" "1.74.0")
find_package(Boost 1.70 COMPONENTS graph REQUIRED)
if(Boost_FOUND)
  include_directories(SYSTEM ${Boost_INCLUDE_DIR})
  link_directories(${Boost_LIBRARY_DIRS})
  message(STATUS "Boost library version ${Boost_LIB_VERSION} found, with headers at '${Boost_INCLUDE_DIR}' and libraries at '${Boost_LIBRARY_DIRS}' for libraries: \n${Boost_LIBRARIES}")
endif()

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")

add_subdirectory(test)
