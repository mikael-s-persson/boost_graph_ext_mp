cmake_minimum_required(VERSION 2.8)

project(boost_graph_contrib)

enable_testing()

if (WIN32)
  if (MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Ox /W3 /bigobj -D_SCL_SECURE_NO_WARNINGS")
    message(STATUS "Configured compiler options and output directories for MSVC toolset.")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -O3 -Wall -Woverloaded-virtual -Wold-style-cast -Wnon-virtual-dtor")
    set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} --enable-stdcall-fixup")
    message(STATUS "Configured compiler options and output directories for MinGW toolset.")
  endif()
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -std=c++0x -O3 -Wall -Woverloaded-virtual -Wold-style-cast -Wnon-virtual-dtor")
  message(STATUS "Configured compiler options and output directories for *nix GCC toolset.")
endif()



set(Boost_ADDITIONAL_VERSIONS "1.45" "1.45.0" "1.46" "1.46.0" "1.46.1" "1.47" "1.47.0" "1.48" "1.48.0" "1.49" "1.49.0" "1.50" "1.50.0" "1.51" "1.51.0" "1.52" "1.52.0" "1.53" "1.53.0")
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)

if (NOT WIN32)
  # make sure that the *nix suffixes and prefixes are correct (some cmake installs of findBoost.cmake are wrong with this).
  set(_ORIGINAL_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})
  set(_ORIGINAL_CMAKE_FIND_LIBRARY_PREFIXES ${CMAKE_FIND_LIBRARY_PREFIXES})
  if( Boost_USE_STATIC_LIBS )
    set(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
  else()
    set(CMAKE_FIND_LIBRARY_SUFFIXES .so ${CMAKE_FIND_LIBRARY_SUFFIXES})
  endif()
  set(CMAKE_FIND_LIBRARY_PREFIXES lib ${CMAKE_FIND_LIBRARY_PREFIXES})
endif()

find_package(Boost 1.48 COMPONENTS unit_test_framework REQUIRED)
if(Boost_FOUND)
  include_directories(SYSTEM ${Boost_INCLUDE_DIR})
  link_directories(${Boost_LIBRARY_DIRS})
  message(STATUS "Boost library version ${Boost_LIB_VERSION} found, with headers at '${Boost_INCLUDE_DIR}' and libraries at '${Boost_LIBRARY_DIRS}' for libraries: \n${Boost_LIBRARIES}")
  if(MSVC)
    #disable the libraries, since it uses automatic linking:
        set(Boost_LIBRARIES "")
  endif()
endif()

if( NOT WIN32 )
  set(CMAKE_FIND_LIBRARY_SUFFIXES ${_ORIGINAL_CMAKE_FIND_LIBRARY_SUFFIXES})
  set(CMAKE_FIND_LIBRARY_PREFIXES ${_ORIGINAL_CMAKE_FIND_LIBRARY_PREFIXES})
endif()




set(BOOST_GRAPH_CONTRIB_INCLUDE_DIRS "../../../")

include_directories(BEFORE SYSTEM ${BOOST_INCLUDE_DIRS})
include_directories(BEFORE ${BOOST_GRAPH_CONTRIB_INCLUDE_DIRS})
include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})


set(EDGE_STORAGE_TYPES 
  vecBC poolBC listBC 
  setBC unordered_setBC multisetBC unordered_multisetBC)

foreach(EDGE_STORAGE_TYPE ${EDGE_STORAGE_TYPES})
  set(TEST_TREES_TNAME "test_trees_t${EDGE_STORAGE_TYPE}")
  configure_file("test_trees_t.cpp.in" "${TEST_TREES_TNAME}.cpp")
  add_executable(${TEST_TREES_TNAME} "${TEST_TREES_TNAME}.cpp")
  add_test(NAME ${TEST_TREES_TNAME} COMMAND "./${TEST_TREES_TNAME}")
  target_link_libraries(${TEST_TREES_TNAME} ${Boost_LIBRARIES})
  
  set(TEST_TREES_ANAME "test_trees_a${EDGE_STORAGE_TYPE}")
  configure_file("test_trees_a.cpp.in" "${TEST_TREES_ANAME}.cpp")
  add_executable(${TEST_TREES_ANAME} "${TEST_TREES_ANAME}.cpp")
  add_test(NAME ${TEST_TREES_ANAME} COMMAND "./${TEST_TREES_ANAME}")
  target_link_libraries(${TEST_TREES_ANAME} ${Boost_LIBRARIES})
  
  set(TEST_GRAPHS_ANAME "test_graphs_a${EDGE_STORAGE_TYPE}")
  configure_file("test_graphs_a.cpp.in" "${TEST_GRAPHS_ANAME}.cpp")
  add_executable(${TEST_GRAPHS_ANAME} "${TEST_GRAPHS_ANAME}.cpp")
  add_test(NAME ${TEST_GRAPHS_ANAME} COMMAND "./${TEST_GRAPHS_ANAME}")
  target_link_libraries(${TEST_GRAPHS_ANAME} ${Boost_LIBRARIES})
  
  set(TEST_PROPMAPS_ANAME "test_prop_maps_a${EDGE_STORAGE_TYPE}")
  configure_file("test_prop_maps_a.cpp.in" "${TEST_PROPMAPS_ANAME}.cpp")
  add_executable(${TEST_PROPMAPS_ANAME} "${TEST_PROPMAPS_ANAME}.cpp")
  add_test(NAME ${TEST_PROPMAPS_ANAME} COMMAND "./${TEST_PROPMAPS_ANAME}")
  target_link_libraries(${TEST_PROPMAPS_ANAME} ${Boost_LIBRARIES})
  
  set(TEST_PROPMAPS_TNAME "test_prop_maps_t${EDGE_STORAGE_TYPE}")
  configure_file("test_prop_maps_t.cpp.in" "${TEST_PROPMAPS_TNAME}.cpp")
  add_executable(${TEST_PROPMAPS_TNAME} "${TEST_PROPMAPS_TNAME}.cpp")
  add_test(NAME ${TEST_PROPMAPS_TNAME} COMMAND "./${TEST_PROPMAPS_TNAME}")
  target_link_libraries(${TEST_PROPMAPS_TNAME} ${Boost_LIBRARIES})
  
endforeach(EDGE_STORAGE_TYPE)

add_executable(test_trees_bfl "test_trees_bfl.cpp")
add_test(NAME test_trees_bfl COMMAND "./test_trees_bfl")
target_link_libraries(test_trees_bfl ${Boost_LIBRARIES})

add_executable(test_prop_maps_bfl "test_prop_maps_bfl.cpp")
add_test(NAME test_prop_maps_bfl COMMAND "./test_prop_maps_bfl")
target_link_libraries(test_prop_maps_bfl ${Boost_LIBRARIES})

add_executable(test_trees_vebl "test_trees_vebl.cpp")
add_test(NAME test_trees_vebl COMMAND "./test_trees_vebl")
target_link_libraries(test_trees_vebl ${Boost_LIBRARIES})

add_executable(test_prop_maps_vebl "test_prop_maps_vebl.cpp")
add_test(NAME test_prop_maps_vebl COMMAND "./test_prop_maps_vebl")
target_link_libraries(test_prop_maps_vebl ${Boost_LIBRARIES})


add_executable(basic_vebl_tests "basic_vebl_tests.cpp")
add_test(NAME basic_vebl_tests COMMAND "./basic_vebl_tests")
target_link_libraries(basic_vebl_tests ${Boost_LIBRARIES})

add_executable(dummy_tests "dummy_tests.cpp")
add_test(NAME dummy_tests COMMAND "./dummy_tests")
target_link_libraries(dummy_tests ${Boost_LIBRARIES})






