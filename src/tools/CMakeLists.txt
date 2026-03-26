include(recipes/polyscope)
include(recipes/libigl)

# Change the location of the compiled programs
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tools)

# Set enabled the Polyscope headless mode -DPOLYSCOPE_BACKEND_OPENGL3_EGL=ON
set(ENV{POLYSCOPE_BACKEND_OPENGL3_EGL} "1")

find_package(PNG)

add_executable(comparePTS comparePTS.cpp)
target_include_directories(comparePTS  PRIVATE "${CMAKE_SOURCE_DIR}/deps/")
target_include_directories(comparePTS  PRIVATE "${libigl_SOURCE_DIR}/include")
target_link_libraries(comparePTS spdlog::spdlog Eigen3::Eigen)

add_executable(displayPTS2 displayPTS2.cpp)
target_include_directories(displayPTS2  PRIVATE "${CMAKE_SOURCE_DIR}/deps/")
target_include_directories(displayPTS2  PRIVATE "${DGTAL_INCLUDE_DIRS}")
target_include_directories(displayPTS2  PRIVATE "${libigl_SOURCE_DIR}/include")
target_link_libraries(displayPTS2 ${DGTAL_LIBRARIES} Eigen3::Eigen polyscope spdlog::spdlog nanoflann::nanoflann)
if(PNG_FOUND)
    target_link_libraries(displayPTS2 PNG::PNG)
endif()

add_executable(displayPTS_mesh displayPTS_mesh.cpp)
target_include_directories(displayPTS_mesh  PRIVATE "${CMAKE_SOURCE_DIR}/deps/")
target_include_directories(displayPTS_mesh  PRIVATE "${DGTAL_INCLUDE_DIRS}")
target_include_directories(displayPTS_mesh  PRIVATE "${libigl_SOURCE_DIR}/include")
target_link_libraries(displayPTS_mesh ${DGTAL_LIBRARIES} Eigen3::Eigen polyscope spdlog::spdlog nanoflann::nanoflann)


include(recipes/nanoflann)
include(recipes/ponca)
add_executable(rescale_pointCloud rescale_pointCloud.cpp)
target_include_directories(rescale_pointCloud  PRIVATE "${CMAKE_SOURCE_DIR}/deps/")
target_include_directories(rescale_pointCloud  PRIVATE "${DGTAL_INCLUDE_DIRS}")
target_include_directories(rescale_pointCloud  PRIVATE "${libigl_SOURCE_DIR}/include")
target_link_libraries(rescale_pointCloud ${DGTAL_LIBRARIES} Ponca::Fitting Eigen3::Eigen polyscope spdlog::spdlog nanoflann::nanoflann)


add_executable(CADprocessing CADprocessing.cpp)
target_include_directories(CADprocessing  PRIVATE "${CMAKE_SOURCE_DIR}/deps/")
target_include_directories(CADprocessing  PRIVATE "${DGTAL_INCLUDE_DIRS}")
target_include_directories(CADprocessing  PRIVATE "${libigl_SOURCE_DIR}/include")
target_link_libraries(CADprocessing ${DGTAL_LIBRARIES} Ponca::Fitting Eigen3::Eigen spdlog::spdlog nanoflann::nanoflann)

add_executable(addNormalsFromFile addNormalsFromFile.cpp)
target_include_directories(addNormalsFromFile  PRIVATE "${CMAKE_SOURCE_DIR}/deps/")
target_include_directories(addNormalsFromFile  PRIVATE "${libigl_SOURCE_DIR}/include")
target_link_libraries(addNormalsFromFile ${DGTAL_LIBRARIES} Eigen3::Eigen spdlog::spdlog )

add_executable(measure_normal_noise measure_normal_noise.cpp)
target_include_directories(measure_normal_noise  PRIVATE "${CMAKE_SOURCE_DIR}/deps/")
target_include_directories(measure_normal_noise  PRIVATE "${libigl_SOURCE_DIR}/include")
target_link_libraries(measure_normal_noise ${DGTAL_LIBRARIES} Eigen3::Eigen spdlog::spdlog )

add_executable(meshCleaner meshCleaner.cpp)
target_include_directories(meshCleaner  PRIVATE "${CMAKE_SOURCE_DIR}/deps/")
target_include_directories(meshCleaner  PRIVATE "${libigl_SOURCE_DIR}/include")
target_link_libraries(meshCleaner ${DGTAL_LIBRARIES} Ponca::Fitting polyscope Eigen3::Eigen spdlog::spdlog )
