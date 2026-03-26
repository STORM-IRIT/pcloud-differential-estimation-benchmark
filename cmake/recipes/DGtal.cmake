if(TARGET DGtal)
    return()
endif()
include(CPM)

CPMAddPackage(
  NAME DGtal
  GITHUB_REPOSITORY DGtal-team/DGtal
  GIT_TAG 07a845b27c91148cc1459d95e1e4bffc0b512e9f
  OPTIONS
    "BUILD_EXAMPLES OFF"
)
include("${DGtal_BINARY_DIR}/DGtalConfig.cmake")

