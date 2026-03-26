# Define version of C++
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(recipes/nanoflann)
include(recipes/cgal)
include(recipes/ponca)
include(recipes/libigl)

# set the libigl include directory
set(LIBIGL_INCLUDE_DIR "${FETCHCONTENT_BASE_DIR}/libigl/include/")

# Change the location of the compiled programs
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/estimators)

SET(ESTIM_SRCs
  compute_jetfitting
  compute_ASO
  compute_ponca_estimators
  compute_ponca_neighbors
  compute_PCPNet_error
  )

find_package(OpenMP)

foreach(FILE ${ESTIM_SRCs})
add_executable(${FILE} ${FILE}.cpp)
target_include_directories(${FILE}  PRIVATE "${CMAKE_SOURCE_DIR}/deps/")
target_include_directories(${FILE}  PRIVATE "${DGTAL_INCLUDE_DIRS}")
target_include_directories(${FILE}  PRIVATE "${libigl_SOURCE_DIR}/include")
if(OpenMP_CXX_FOUND)
  target_link_libraries(${FILE} OpenMP::OpenMP_CXX)
endif()
endforeach()
  
  # Specific target deps
  target_link_libraries(compute_jetfitting ${DGTAL_LIBRARIES} CGAL::CGAL_Core Eigen3::Eigen polyscope spdlog::spdlog nanoflann::nanoflann)
  target_link_libraries(compute_ASO ${DGTAL_LIBRARIES} Ponca::Fitting Eigen3::Eigen polyscope spdlog::spdlog nanoflann::nanoflann)
  target_link_libraries(compute_ponca_estimators ${DGTAL_LIBRARIES} Ponca::Fitting Eigen3::Eigen polyscope spdlog::spdlog nanoflann::nanoflann)
  target_link_libraries(compute_ponca_neighbors ${DGTAL_LIBRARIES} Ponca::Fitting Eigen3::Eigen polyscope spdlog::spdlog nanoflann::nanoflann)
  target_link_libraries(compute_PCPNet_error ${DGTAL_LIBRARIES} CGAL::CGAL_Core Eigen3::Eigen polyscope spdlog::spdlog nanoflann::nanoflann)
  
  # GMP
  set(prefix ${FETCHCONTENT_BASE_DIR}/gmp)
  set(gmp_INSTALL ${prefix}/install)
  set(gmp_LIB_DIR ${gmp_INSTALL}/lib)
  set(gmp_LIBRARY 
  ${gmp_LIB_DIR}/${CMAKE_STATIC_LIBRARY_PREFIX}gmp${CMAKE_STATIC_LIBRARY_SUFFIX}
  ${gmp_LIB_DIR}/${CMAKE_STATIC_LIBRARY_PREFIX}gmpxx${CMAKE_STATIC_LIBRARY_SUFFIX}
  )
  set(gmp_LIBRARIES ${gmp_LIBRARY})
  set(gmp_INCLUDE_DIR ${gmp_INSTALL}/include)

  # compute_VCM
  add_executable(compute_VCM compute_VCM.cpp)
  target_include_directories(compute_VCM  PRIVATE "${CMAKE_SOURCE_DIR}/deps/")
  target_include_directories(compute_VCM  PRIVATE "${DGTAL_INCLUDE_DIRS}")
  target_include_directories(compute_VCM  PRIVATE "${gmp_INCLUDE_DIR}")
  target_include_directories(compute_VCM  PRIVATE "${libigl_SOURCE_DIR}/include")
  target_link_libraries(compute_VCM ${DGTAL_LIBRARIES} CGAL::CGAL_Core Eigen3::Eigen polyscope spdlog::spdlog nanoflann::nanoflann ${gmp_LIBRARIES})