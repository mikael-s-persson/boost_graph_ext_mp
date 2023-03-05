
include_directories("${CMAKE_CURRENT_SOURCE_DIR}")

set(EDGE_STORAGE_TYPES 
  vecBC poolBC listBC 
  setBC unordered_setBC multisetBC unordered_multisetBC)

foreach(EDGE_STORAGE_TYPE ${EDGE_STORAGE_TYPES})
  set(TEST_TREES_TNAME "test_trees_t${EDGE_STORAGE_TYPE}")
  configure_file("test_trees_t.cpp.in" "${TEST_TREES_TNAME}.cpp")
  add_executable(${TEST_TREES_TNAME} "${TEST_TREES_TNAME}.cpp")
  add_test(NAME ${TEST_TREES_TNAME} COMMAND "./${TEST_TREES_TNAME}")
  target_link_libraries(${TEST_TREES_TNAME} GTest::gtest_main)
  
  set(TEST_TREES_ANAME "test_trees_a${EDGE_STORAGE_TYPE}")
  configure_file("test_trees_a.cpp.in" "${TEST_TREES_ANAME}.cpp")
  add_executable(${TEST_TREES_ANAME} "${TEST_TREES_ANAME}.cpp")
  add_test(NAME ${TEST_TREES_ANAME} COMMAND "./${TEST_TREES_ANAME}")
  target_link_libraries(${TEST_TREES_ANAME} GTest::gtest_main)
  
  set(TEST_GRAPHS_ANAME "test_graphs_a${EDGE_STORAGE_TYPE}")
  configure_file("test_graphs_a.cpp.in" "${TEST_GRAPHS_ANAME}.cpp")
  add_executable(${TEST_GRAPHS_ANAME} "${TEST_GRAPHS_ANAME}.cpp")
  add_test(NAME ${TEST_GRAPHS_ANAME} COMMAND "./${TEST_GRAPHS_ANAME}")
  target_link_libraries(${TEST_GRAPHS_ANAME} GTest::gtest_main)
  
  set(TEST_PROPMAPS_ANAME "test_prop_maps_a${EDGE_STORAGE_TYPE}")
  configure_file("test_prop_maps_a.cpp.in" "${TEST_PROPMAPS_ANAME}.cpp")
  add_executable(${TEST_PROPMAPS_ANAME} "${TEST_PROPMAPS_ANAME}.cpp")
  add_test(NAME ${TEST_PROPMAPS_ANAME} COMMAND "./${TEST_PROPMAPS_ANAME}")
  target_link_libraries(${TEST_PROPMAPS_ANAME} GTest::gtest_main)
  
  set(TEST_BUNDLES_ANAME "test_bundles_a${EDGE_STORAGE_TYPE}")
  configure_file("test_bundles_a.cpp.in" "${TEST_BUNDLES_ANAME}.cpp")
  add_executable(${TEST_BUNDLES_ANAME} "${TEST_BUNDLES_ANAME}.cpp")
  add_test(NAME ${TEST_BUNDLES_ANAME} COMMAND "./${TEST_BUNDLES_ANAME}")
  target_link_libraries(${TEST_BUNDLES_ANAME} GTest::gtest_main)
  
  set(TEST_PROPMAPS_TNAME "test_prop_maps_t${EDGE_STORAGE_TYPE}")
  configure_file("test_prop_maps_t.cpp.in" "${TEST_PROPMAPS_TNAME}.cpp")
  add_executable(${TEST_PROPMAPS_TNAME} "${TEST_PROPMAPS_TNAME}.cpp")
  add_test(NAME ${TEST_PROPMAPS_TNAME} COMMAND "./${TEST_PROPMAPS_TNAME}")
  target_link_libraries(${TEST_PROPMAPS_TNAME} GTest::gtest_main)
  
  set(TEST_BUNDLES_TNAME "test_bundles_t${EDGE_STORAGE_TYPE}")
  configure_file("test_bundles_t.cpp.in" "${TEST_BUNDLES_TNAME}.cpp")
  add_executable(${TEST_BUNDLES_TNAME} "${TEST_BUNDLES_TNAME}.cpp")
  add_test(NAME ${TEST_BUNDLES_TNAME} COMMAND "./${TEST_BUNDLES_TNAME}")
  target_link_libraries(${TEST_BUNDLES_TNAME} GTest::gtest_main)
  
endforeach(EDGE_STORAGE_TYPE)

add_executable(test_trees_bfl "test_trees_bfl.cpp")
add_test(NAME test_trees_bfl COMMAND "./test_trees_bfl")
target_link_libraries(test_trees_bfl GTest::gtest_main)

add_executable(test_prop_maps_bfl "test_prop_maps_bfl.cpp")
add_test(NAME test_prop_maps_bfl COMMAND "./test_prop_maps_bfl")
target_link_libraries(test_prop_maps_bfl GTest::gtest_main)

add_executable(test_bundles_bfl "test_bundles_bfl.cpp")
add_test(NAME test_bundles_bfl COMMAND "./test_bundles_bfl")
target_link_libraries(test_bundles_bfl GTest::gtest_main)

add_executable(test_trees_vebl "test_trees_vebl.cpp")
add_test(NAME test_trees_vebl COMMAND "./test_trees_vebl")
target_link_libraries(test_trees_vebl GTest::gtest_main)

add_executable(test_prop_maps_vebl "test_prop_maps_vebl.cpp")
add_test(NAME test_prop_maps_vebl COMMAND "./test_prop_maps_vebl")
target_link_libraries(test_prop_maps_vebl GTest::gtest_main)

add_executable(test_bundles_vebl "test_bundles_vebl.cpp")
add_test(NAME test_bundles_vebl COMMAND "./test_bundles_vebl")
target_link_libraries(test_bundles_vebl GTest::gtest_main)

add_executable(basic_vebl_tests "basic_vebl_tests.cpp")
add_test(NAME basic_vebl_tests COMMAND "./basic_vebl_tests")

add_executable(dummy_tests "dummy_tests.cpp")
add_test(NAME dummy_tests COMMAND "./dummy_tests")






