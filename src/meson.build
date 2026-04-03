# ===========================================================================
# Build KRONROD C library as a static archive
# ===========================================================================

kronrod_sources = files(
    'kronrod.c',
)

kronrod_inc = include_directories('.')

kronrod_lib = static_library('kronrod',
    kronrod_sources,
    include_directories: kronrod_inc,
    dependencies: m_dep,
)
