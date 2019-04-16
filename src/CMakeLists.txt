# CMakeList.txt : CMake project for UKF2, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

# Add source to this project's executable.
add_executable (UKF2 "UKF2.cpp" "UKF2.h" "ukf.cpp" "ukf.h" "tools.cpp" "tools.h" "ground_truth_package.h" "measurement_package.h")

# TODO: Add tests and install targets if needed.
target_link_libraries(UKF2)