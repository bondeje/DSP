#gcc -pedantic Wall -O3 -g -DDEBUG_MALLOC RTFilter.c polynomials.c Lpolys.c legendre.c chebyshev.c ../mallocs/debug_malloc.c test_rtfilters.c -o f.exe
UNAME := $(shell uname)
CC = gcc

MAIN_NAME = dsp

BUILD_DIR = ../bin/
WORK_DIR = 
PYLIB_DIR = 
LIB_DIR = ../lib/

EXT = 
LIB_NAME = 

LFLAGS = #-L./ 
CFLAGS = -std=c99 -O3 -Wall -pedantic -fPIC -shared -DDLL_EXPORT
# really cool, -g creates symbols so that valgrind will actually show you the lines of errors
#CFLAGS += -g # for debugging only
IFLAGS = -I../include
CCFLAGS = #-D DEVELOPMENT # -D DEVELOPMENT turns on internal print statements for diagnostics

ifeq ($(OS),Windows_NT)
	# might have to encapsulate with a check for MINGW. Need this because Windows f-s up printf with size_t and MINGW only handles it with their own implementation of stdio
	CFLAGS += -D__USE_MINGW_ANSI_STDIO
    CCFLAGS += -D WIN32
	EXT = .dll
	LIB_NAME = $(MAIN_NAME)
    COPY_CMD = copy /y
    COPY_SRC = ..\bin\\
    WORK_DIR = ..\workbench\\
    PYLIB_DIR = ..\python\Lib\\
    #ifeq ($(PROCESSOR_ARCHITECTURE),AMD64)
    #    CCFLAGS += -D AMD64
    #else
    #    ifeq ($(PROCESSOR_ARCHITECTURE),AMD64)
    #        CCFLAGS += -D AMD64
    #    endif
    #    ifeq ($(PROCESSOR_ARCHITECTURE),x86)
    #        CCFLAGS += -D IA32
    #    endif
    #endif
else
    UNAME_S := $(shell uname -s)
	# for dynamic memory allocation extensions in posix, e.g. getline()
	#CFLAGS += -D__STDC_WANT_LIB_EXT2__=1
	EXT = .so
	LIB_NAME = lib$(MAIN_NAME)
    ifeq ($(UNAME_S),Linux)
		# needed because linux must link to the separate math library
		LFLAGS += -lm
        CCFLAGS += -D LINUX
    endif
    WORK_DIR = ../workbench/
    PYLIB_DIR = ../python/Lib/
    COPY_CMD = cp
    COPY_SRC = $(BUILD_DIR)
    #ifeq ($(UNAME_S),Darwin)
    #    CCFLAGS += -D OSX
    #endif
    #UNAME_P := $(shell uname -p)
    #ifeq ($(UNAME_P),x86_64)
    #    CCFLAGS += -D AMD64
    #endif
    #ifneq ($(filter %86,$(UNAME_P)),)
    #    CCFLAGS += -D IA32
    #endif
    #ifneq ($(filter arm%,$(UNAME_P)),)
    #    CCFLAGS += -D ARM
    #endif
endif

LFLAGS += -Wl,--out-implib,$(LIB_DIR)$(LIB_NAME).lib

CFLAGS += -o $(BUILD_DIR)$(LIB_NAME)$(EXT)

all: build install

build:
	$(CC) $(CFLAGS) $(IFLAGS) $(CCFLAGS) RTFilter.c polynomials.c Lpolys.c legendre.c chebyshev.c hermite.c laguerre.c filterutils.c $(LFLAGS)

install: install_pylib install_workbench

install_pylib:
	$(COPY_CMD) $(COPY_SRC)$(LIB_NAME)$(EXT) $(PYLIB_DIR)

install_workbench:
	$(COPY_CMD) $(COPY_SRC)$(LIB_NAME)$(EXT) $(WORK_DIR)
