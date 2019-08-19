cmake_minimum_required(VERSION 2.6)
set(PROJ_NAME CPM)
project (${PROJ_NAME})

find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})

#add_definitions(-DUSE_DAISY)
add_definitions(-DWITH_SSE)

include_directories(./include)

add_executable(${PROJ_NAME} main.cpp CPM.cpp)
target_link_libraries(${PROJ_NAME} ${OpenCV_LIBS})
