
###########################################################################
# Note: This Makefile is intended for development.  For production
# builds, use build_dynamic/Makefile.
###########################################################################

# Does not currently support -DCPU_CHECK_...

# BASEFLAGS=-flto -ffp-contract=off -DZSTD_MULTITHREAD -DIGNORE_BUNDLED_ZSTD
BASEFLAGS=-g -ffp-contract=off -DZSTD_MULTITHREAD -DIGNORE_BUNDLED_ZSTD
# BASEFLAGS=-flto -mavx2 -mbmi -mbmi2 -mfma -mlzcnt -DZSTD_MULTITHREAD
# BASEFLAGS=-flto -msse4.2 -DZSTD_MULTITHREAD

include Makefile.src

# Respect the environment:
# Use defaults below only if not set in env or make arguments
CC		?= gcc
CXX		?= g++
CFLAGS		?= -O3
CXXFLAGS	?= -O2
ZLIB		?= -lz
# ZLIB		?= -L. ../zlib-${ZLIB_VER}/libz.a

# dynamic zstd: remove -DIGNORE_BUNDLED_ZSTD above, change this to -lzstd,
# replace OBJ with OBJ_NO_ZSTD below
ZSTD		?=

BLASFLAGS64	?= -llapack -lf77blas -latlas
ARCH32		?=

# Mandatory flags added to defaults or env settings
CFLAGS		+= -std=gnu99 $(BASEFLAGS) $(CWARN) $(CINCLUDE)
CXXFLAGS	+= -std=c++17 $(BASEFLAGS) $(CXXWARN) $(CXXINCLUDE)
LDFLAGS		+= -lm -lpthread -flto $(ZLIB) $(ZSTD)

# Installation defaults
MKDIR		?= mkdir
INSTALL		?= install
STRIP_CMD	?= strip
PREFIX		?= /usr/local
DESTDIR		?= .

UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
  # argh
  CFLAGS += -I/opt/homebrew/include -I/usr/local/include
  CXXFLAGS += -I/opt/homebrew/include -I/usr/local/include -DACCELERATE_NEW_LAPACK
  BLASFLAGS=-framework Accelerate
  BLASFLAGS64=-framework Accelerate
  LDFLAGS=-ldl -lpthread -lz -lzstd -flto -L/opt/homebrew/lib -L/usr/local/lib
endif

%.o: %.c
	$(CC) -c $(CFLAGS) $(ARCH32) -o $@ $<

%.o: %.cc
	$(CXX) -c $(CXXFLAGS) $(ARCH32) -o $@ $<

all: adjutant

# for clean build, "make clean" first
adjutant: $(OBJ)
	$(MKDIR) -p bin
	$(CXX) $(ARCH32) $(OBJ) -o bin/adjutant $(BLASFLAGS64) $(LDFLAGS)

pgs_eval_test: adjutant_pgs_eval.o adjutant_pgs_eval_test.o $(OBJ)
	mkdir -p bin
	$(CXX) $(ARCH32) adjutant_pgs_eval_test.o $(filter-out adjutant.o,$(OBJ)) -o bin/pgs_eval_test $(BLASFLAGS64) $(LDFLAGS)

# Optimised build.  Two things to preserve here:
#  * -ffp-contract=off is NOT optional: it stops the compiler fusing
#    multiply-add in the scoring dot product, which keeps results bit-identical
#    to stock plink2 across optimisation levels.
#  * -O3 goes in BASEFLAGS, not CXXFLAGS.  Setting CXXFLAGS on the command line
#    overrides the whole variable, discarding the platform block below that adds
#    -DACCELERATE_NEW_LAPACK and the homebrew include paths.
.PHONY: release
release:
	$(MAKE) clean
	$(MAKE) -j8 \
	  BASEFLAGS="-O3 -DNDEBUG -ffp-contract=off -DZSTD_MULTITHREAD" \
	  LDFLAGS="-ldl -lpthread -lz -flto" adjutant

.PHONY: install-strip install clean

install-strip: install
	$(STRIP_CMD) $(DESTDIR)$(PREFIX)/bin/*

install: all
	$(MKDIR) -p $(DESTDIR)$(PREFIX)/bin
	$(INSTALL) -c bin/* $(DESTDIR)$(PREFIX)/bin

clean:
	rm -f $(CLEAN)
