cmake_minimum_required(VERSION 3.20)
project(grover_search)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(SDL2 REQUIRED)

add_executable(grover_search main.cpp)

target_include_directories(grover_search PRIVATE ${SDL2_INCLUDE_DIRS})
target_link_libraries(grover_search PRIVATE ${SDL2_LIBRARIES})
