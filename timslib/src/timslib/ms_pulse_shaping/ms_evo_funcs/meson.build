# Dependencies for numpy
incdir_numpy = run_command(py3, [
    '-c',
    '''import os, numpy as np
try:
    print(os.path.relpath(np.get_include()))
except:
    print(np.get_include())
    ''',
], check: true).stdout().strip()

inc_np = include_directories(incdir_numpy)
incdir_f2py = incdir_numpy / '..' / '..' / 'f2py' / 'src'
inc_f2py = include_directories(incdir_f2py)
fortranobject_c = incdir_f2py / 'fortranobject.c'

# Dependencies for Fortran
fortranobject_lib = static_library(
    '_fortranobject',
    fortranobject_c,
    dependencies: py3_dep,
    include_directories: [inc_np, inc_f2py],
)
fortranobject_dep = declare_dependency(
    link_with: fortranobject_lib,
    include_directories: [inc_np, inc_f2py],
)


# ppoly_alpha_and_chi module
ppoly_alpha_and_chi_module = custom_target(
    'ppoly_alpha_and_chimodule',
    output : 'ppoly_alpha_and_chimodule.c',
    input : 'ppoly_alpha_and_chi.pyf',
    command: [py3, '-m', 'numpy.f2py', '@INPUT@', '--build-dir', '@OUTDIR@']
)

py3.extension_module(
    'ppoly_alpha_and_chi',
    [ppoly_alpha_and_chi_module, 'ppoly_alpha_and_chi.f95', 'ppoly_alpha_and_chi_mod.f95', 'expr_for_ppoly_alpha_and_chi.f95'],
    dependencies: fortranobject_dep,
    install : true,
    link_language: 'fortran',
    subdir: 'timslib' / 'ms_pulse_shaping' / 'ms_evo_funcs',
)


ppoly_seg_matrices_module = custom_target(
    'ppoly_seg_matricesmodule',
    output : 'ppoly_seg_matricesmodule.c',
    input : 'ppoly_seg_matrices.pyf',
    command: [py3, '-m', 'numpy.f2py', '@INPUT@', '--build-dir', '@OUTDIR@']
)
py3.extension_module(
    'ppoly_seg_matrices',
    [ppoly_seg_matrices_module, 'ppoly_seg_matrices.f95', 'expr_for_ppoly_alpha_and_chi.f95'],
    dependencies: fortranobject_dep,
    install : true,
    link_language: 'fortran',
    subdir: 'timslib' / 'ms_pulse_shaping' / 'ms_evo_funcs',
)


num_ppoly_alpha_and_chi_module = custom_target(
    'num_ppoly_alpha_and_chimodule',
    output : 'num_ppoly_alpha_and_chimodule.c',
    input : 'num_ppoly_alpha_and_chi.pyf',
    command: [py3, '-m', 'numpy.f2py', '@INPUT@', '--build-dir', '@OUTDIR@']
)

py3.extension_module(
    'num_ppoly_alpha_and_chi',
    [num_ppoly_alpha_and_chi_module, 'num_ppoly_alpha_and_chi.f95', 'ppoly_alpha_and_chi_mod.f95', 'expr_for_ppoly_alpha_and_chi.f95'],
    dependencies: fortranobject_dep,
    install : true,
    link_language: 'fortran',
    subdir: 'timslib' / 'ms_pulse_shaping' / 'ms_evo_funcs',
)


ipowexp_module = custom_target(
    'ipowexpmodule',
    output : 'ipowexpmodule.c',
    input : 'ipowexp.pyf',
    command: [py3, '-m', 'numpy.f2py', '@INPUT@', '--build-dir', '@OUTDIR@']
)


py3.extension_module(
    'ipowexp',
    [ipowexp_module, 'ipowexp.f95', 'expr_for_ppoly_alpha_and_chi.f95'],
    dependencies: fortranobject_dep,
    install : true,
    link_language: 'fortran',
    subdir: 'timslib' / 'ms_pulse_shaping' / 'ms_evo_funcs',
)


# Python files
py3.install_sources([
    '__init__.py',
], subdir: 'timslib' / 'ms_pulse_shaping' / 'ms_evo_funcs')

