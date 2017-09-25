cmake_minimum_required(VERSION 2.8)

####
# project settings
####

project(Correlation)

set(sample_VERSION_MAJOR 1)
set(sample_VERSION_MINOR 0)


####
# headers
####

set(HEADERS
)


####
# sources
####

set(SOURCES
	correlation_func.cpp
)


####
# dependencies
####

find_package(GLAMER NO_MODULE REQUIRED)

include_directories(${GLAMER_INCLUDE_DIRS})


####
# build
####

add_executable(Correlation ${HEADERS} ${SOURCES})

target_link_libraries(Correlation ${GLAMER_LIBRARIES})


