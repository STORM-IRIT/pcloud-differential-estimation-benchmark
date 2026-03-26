# -----------------------------------------------------------------------------
# ponca
# -----------------------------------------------------------------------------

include(CPM)
CPMAddPackage(
        NAME ponca
        GITHUB_REPOSITORY poncateam/ponca
#        GIT_TAG master
        GIT_TAG 19ef58f # commit from ~20/03/2026 just before changing addLocalNeighbor return type to void
        OPTIONS "PONCA_CONFIGURE_EXAMPLES OFF" "PONCA_CONFIGURE_DOC OFF" "PONCA_CONFIGURE_TESTS OFF"
)