cmake_minimum_required(VERSION 3.14)
project(optionpricer CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Default to an optimised build (Monte Carlo will be slow in Debug).
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# The library: all pricing code lives here.
add_library(optionpricer src/black_scholes.cpp)
target_include_directories(optionpricer PUBLIC include)
target_compile_options(optionpricer PRIVATE -Wall -Wextra)

# Small demo program for step 1.
add_executable(bs_example apps/bs_example.cpp)
target_link_libraries(bs_example PRIVATE optionpricer)
