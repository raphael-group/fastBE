enable_testing()

add_executable(piecewiselinear_test piecewiselinearf_tests.cpp)

target_link_libraries(piecewiselinear_test gtest gtest_main)

# add_test(NAME piecewiselinear COMMAND piecewiselinear_test)
include(GoogleTest)
gtest_discover_tests(piecewiselinear_test)
