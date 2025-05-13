'''
'''

from setuptools import setup, Extension

# ext = Extension(
#     "testmod",
#     sources=["src/test_pymodule.cpp"],
#     language="c++",
#     extra_compile_args=["-std=c++23"],
# )

# setup(
#     name="testmod",
#     ext_modules=[ext]
# )

setup(
    name="custom",
    version="0.1.0",
    ext_modules=[Extension("custom", ["src/custom_pymodule.cpp"])]
)