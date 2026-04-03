# ===========================================================================
# KRONROD test suite
# ===========================================================================

test_kronrod = executable('test_kronrod',
    'test_kronrod.c',
    include_directories: kronrod_inc,
    link_with: kronrod_lib,
    dependencies: m_dep,
)

test('kronrod_tests', test_kronrod)
