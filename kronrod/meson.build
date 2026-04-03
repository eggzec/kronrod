
_kronrod_module = custom_target('_kronrodmodule',
    input: '_kronrod.pyf',
    output: ['_kronrodmodule.c'],
    command: [py, '-m', 'numpy.f2py', '@INPUT@', '--build-dir', '@OUTDIR@'],
)

py.extension_module('_kronrod',
    [_kronrod_module, '_kronrod_helpers.c', fortranobject_c],
    include_directories: kronrod_inc,
    c_args: np_compile_args,
    link_with: kronrod_lib,
    dependencies: [py_dep, m_dep],
    install: true,
    subdir: 'kronrod',
)

py.install_sources(
    '__init__.py',
    subdir: 'kronrod',
)
