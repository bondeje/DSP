# DSP
Digital signal processing library

# pre-requisites
-C99 compatible compiler. Tested on `gcc` from MSYS2/MinGW64 in Windows and Linux/Ubuntu

-Make. Tested with GNU/make

# build dll
```
~>git clone https://github.com/bondeje/DSP.git
```
or however you are able to copy source from repo.
```
~>cd DSP\src
~\DSP\src>make -f make_dll.mak
```

# python

## Load module
```
~\DSP\src>cd ..\python
~\DSP\python>python
>>>import dsp
```

## Run python example band-pass filter
```
~\DSP\python>python dsp.py
```