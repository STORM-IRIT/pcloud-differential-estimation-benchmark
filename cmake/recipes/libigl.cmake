# -----------------------------------------------------------------------------
# ponca
# -----------------------------------------------------------------------------
# include(CPM)
# CPMAddPackage(
#   NAME ponca
#   GITHUB_REPOSITORY poncateam/ponca
#   GIT_TAG master
#   OPTIONS "PONCA_CONFIGURE_EXAMPLES OFF" "PONCA_CONFIGURE_DOC OFF" "PONCA_CONFIGURE_TESTS OFF"
# )

include(CPM)
CPMAddPackage(
  NAME libigl
  GITHUB_REPOSITORY libigl/libigl
  GIT_TAG main
  DOWNLOAD_ONLY ON
)

# Make the include directory available
# add_library(libigl INTERFACE)
# target_include_directories(libigl INTERFACE ${libigl_SOURCE_DIR}/include)

message (STATUS "=========================+> libigl_SOURCE_DIR: ${libigl_SOURCE_DIR}")