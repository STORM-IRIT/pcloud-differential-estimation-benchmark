cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
Project(PointCloudDiffBenchmark 
        VERSION 1.0
        DESCRIPTION "Benchmarking point cloud differential quantities estimators."
        LANGUAGES CXX)


set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/recipes)

include(config)

add_subdirectory(src/generators)
add_subdirectory(src/estimators)
add_subdirectory(src/tools)
