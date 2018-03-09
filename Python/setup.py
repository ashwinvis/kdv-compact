from numpy.distutils.core import Extension, setup

ext = Extension(name='finite_diff', sources=['finite_diff.f90'])

setup(
    name="kdv",
    description="Python version of the KdV solver",
    install_requires=['scipy', 'matplotlib'],
    ext_modules=[ext],
    script_name='setup.py',
    script_args=['build_ext', '--inplace']
)
