add_library(geogram STATIC ${CMAKE_SOURCE_DIR}/deps/Delaunay_psm.cpp)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/generators)


SET(GENERATORS_SRCs
        generateFromImplicit
        generate_adaptator_pcpnet
)

foreach(FILE ${GENERATORS_SRCs})
  add_executable(${FILE} ${FILE}.cpp)

  set_target_properties(${FILE} PROPERTIES
          CXX_STANDARD 17
          CXX_STANDARD_REQUIRED ON
  )

  target_link_libraries(${FILE} ${DGTAL_LIBRARIES} spdlog::spdlog geogram)
  target_include_directories(${FILE} PRIVATE "${CMAKE_SOURCE_DIR}/deps/")
  target_include_directories(${FILE} PRIVATE "${DGTAL_INCLUDE_DIRS}")
endforeach()
