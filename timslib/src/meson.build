project(
    'timslib',
    'c',
    version: '1.0.0',
    meson_version: '>= 0.64.0',
)

# https://mesonbuild.com/Python-module.html
py_mod = import('python')
py3 = py_mod.find_installation(pure: false)
py3_dep = py3.dependency()

# Compilers
add_languages('fortran', native: false)

add_project_arguments(
  '-ffree-line-length-none',
  language : 'fortran'
  )

subdir('timslib')

