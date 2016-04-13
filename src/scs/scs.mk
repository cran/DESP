
#CC = $(BINPREF)gcc $(M_ARCH)


CC = `${R_HOME}/bin/R CMD config CC`
LDFLAGS = `${R_HOME}/bin/R CMD config LDFLAGS`
CFLAGS = `${R_HOME}/bin/R CMD config CFLAGS`
SHLIB_EXT = $(shell echo `${R_HOME}/bin/R CMD config SHLIB_EXT`)

SHARED = $(SHLIB_EXT)


ifeq ($(OS),Windows_NT)
	SYS = $(shell $(CC) -dumpmachine)
	ifneq (, $(findstring cygwin, $(SYS)))
		ISWINDOWS := 1
	else 
		ifneq (, $(findstring mingw, $(SYS)))
			ISWINDOWS := 1
		else 
			ifneq (, $(findstring msys, $(SYS)))
				ISWINDOWS := 1
			else
				ISWINDOWS := 0
			endif
		endif
	endif
	UNAME = $(OS)
else
	ISWINDOWS := 0
	UNAME = $(shell uname -s)
endif

ifeq ($(ISWINDOWS), 1)
	# we're on windows (cygwin or msys)
	LDFLAGS += -lm
else 
	ifeq ($(UNAME), Darwin)
		# we're on apple
		LDFLAGS += -lm
	else
		# we're on a linux (or solaris) system, use accurate timer provided by clock_gettime()
		LDFLAGS += -lm -lrt
	endif
endif

# shared library installation path (for OS X)
IN_PATH =

ifeq ($(UNAME), Darwin)
	IN_PATH = -install_name @rpath/libscsdir$(SHARED)
endif

# gcc and clang shared libraries compiler option
SHLIB_COMPILER_OPTION = -shared

# Add on default CFLAGS
ADDITIONAL_CFLAGS = -g -DCTRLC=1 $(CC_WARN) -pedantic -O3 -funroll-loops -Wstrict-prototypes -I. -Iinclude
ifneq ($(ISWINDOWS), 1)
	ADDITIONAL_CFLAGS += -fPIC
endif

# Add on default CFLAGS for Oracle Solaris Studio and modify shared libraries compiler option
ifeq ($(UNAME), SunOS)
	ifeq (, $(findstring gcc, $(CC))) 
		ifeq (, $(findstring clang, $(CC))) 
			# we're on solaris (with Oracle Solaris Studio)
			ADDITIONAL_CFLAGS = -g -DCTRLC=1 -I. -Iinclude
			SHLIB_COMPILER_OPTION = -G
		endif
	endif
endif

CFLAGS += $(ADDITIONAL_CFLAGS)


LINSYS = linsys
DIRSRC = $(LINSYS)/direct
DIRSRCEXT = $(DIRSRC)/external
INDIRSRC = $(LINSYS)/indirect

OUT = out
LIB = lib
AR = ar
ARFLAGS = rv
ARCHIVE = $(AR) $(ARFLAGS)
RANLIB = ranlib

########### OPTIONAL FLAGS ##########
# CFLAGS += -DDLONG # use longs rather than ints
# CFLAGS += -DFLOAT # use floats rather than doubles
# CFLAGS += -DNOVALIDATE # remove data validation step
# CFLAGS += -DBLASSUFFIX="" # blas suffix no underscore
# CFLAGS += -DEXTRAVERBOSE # extra verbosity level
# CFLAGS += -DNOTIMER # no timing, times reported as nan
CFLAGS += -DCOPYAMATRIX # if normalize, copy A

############ OPENMP: ############
# set USE_OPENMP = 1 to allow openmp (multi-threaded matrix multiplies):
# set the number of threads to, for example, 4 by entering the command:
# export OMP_NUM_THREADS=4

USE_OPENMP = 0

ifneq ($(USE_OPENMP), 0)
  CFLAGS += -fopenmp -DOPENMP
  LDFLAGS += -lgomp
endif

############ SDPS: BLAS + LAPACK ############
# set USE_LAPACK = 1 below to enable solving SDPs
# NB: point the libraries to the locations where
# you have blas and lapack installed

USE_LAPACK = 0

ifneq ($(USE_LAPACK), 0)
  # edit these for your setup:
  LDFLAGS += -lblas -llapack #-lgfortran
  CFLAGS += -DLAPACK_LIB_FOUND
  # CFLAGS += -DBLAS64 # if blas/lapack lib uses 64 bit ints
endif

# debug to see var values, e.g. 'make print-OBJECTS' shows OBJECTS value
print-%: ; @echo $*=$($*)
