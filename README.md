# FFTLog-and-beyond

An extended FFTLog code for efficiently computing integrals containing:

* one Bessel function (i.e. Hankel transform); or

* one spherical Bessel function; or

* one 1st or 2nd-derivative of spherical Bessel function.

v2.0: fftlogx

The code is written in C ([./src/](src)) and provides a python wrapper ([./fftlogx/](fftlogx)). To use it, first run
```shell
make
```
to build a shared library `libfftlogx.so`, then follow the test python scripts provided in [/fftlogx/](fftlogx) to import and use it.

-----

The older version (v1.0):

The code is *independently* written and tested in python ([./python/fftlog.py](python/fftlog.py)) and C ([./cfftlog/](cfftlog)).

See more details in [Notes.pdf](Notes.pdf)

Please cite [Fang et al (2019); arXiv:1911.11947](https://arxiv.org/abs/1911.11947), if you find the algorithm or the code useful to your research.

Please feel free to use and adapt the code for your own purpose, and let me know if you are confused or find a bug (just open an [issue](https://github.com/xfangcosmo/FFTLog-and-beyond/issues)). FFTLog-and-beyond is open source and distributed with the
[MIT license](https://opensource.org/licenses/mit).
