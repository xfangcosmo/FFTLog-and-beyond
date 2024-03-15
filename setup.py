import os
from setuptools import setup, find_packages
from setuptools import Distribution as _distribution
from setuptools.command.build_py import build_py as _build
from setuptools.command.develop import develop as _develop
from distutils.command.clean import clean as _clean
from subprocess import check_call, CalledProcessError


def make_fftlogx():
    try:
        check_call(['make'])
    except CalledProcessError as e:
        raise RuntimeError(f"Could not build cfftlog, return code: {e.returncode}")
    check_call(['cp', 'build/libfftlogx.so', 'fftlogx/'])


class Distribution(_distribution):
    global_options = _distribution.global_options

    global_options += [
        ("debug", None, "build on debug mode"),
    ]

    def __init__(self, attr=None):
        self.debug = False
        super().__init__(attr)


class Build(_build):
    def run(self):
        make_fftlogx()
        _build.run(self)


class Develop(_develop):
    def run(self):
        make_fftlogx()
        _develop.run(self)


class Clean(_clean):
    """Remove the copied libfftlogx.so"""

    def run(self):
        if os.path.isfile("fftlogx/libfftlogx.so"):
            os.remove("fftlogx/libfftlogx.so")
            print("Removed fftlogx/libfftlogx.so")
        _clean.run(self)
        check_call(['make', 'clean'])


setup(
    name='fftlogx',
    version='0.1.0',
    description='Python wrapper for the generalized FFTLog algorithm',
    url='https://github.com/xfangcosmo/FFTLog-and-beyond',
    packages=find_packages(exclude=['test']),
    install_requires=['numpy'],
    package_data={
        'fftlogx': ['libfftlogx.so']
    },
    distclass=Distribution,
    cmdclass={'build_py': Build, 'develop': Develop, 'clean': Clean},
)
