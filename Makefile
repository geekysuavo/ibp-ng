
# HAVE_PTHREAD: whether to enable *any* multi-threading, cpu or gpu.
# HAVE_CUDA: whether to enable gpu code. requires HAVE_PTHREAD=y.
IBP_PTHREAD=y
IBP_CUDA=n

# CC: compiler binary filename.
CC=gcc
NVCC=nvcc
LEX=flex
YACC=bison

# LD: linkage binary filename.
ifeq ($(IBP_CUDA),y)
LD=nvcc
else
LD=gcc
endif

# CFLAGS, LFLAGS, YFLAGS, LIBS: compilation flags and library linkage flags.
CFLAGS=-ggdb -O3 -std=c99 -D_GNU_SOURCE
CFLAGS+= -Wall -Wformat -Wextra -Wno-unused-parameter
LFLAGS=
YFLAGS=-d
LIBS=-lm

# CFLAGS, LIBS: pthread/cuda-only compilation flags and libraries.
ifeq ($(IBP_PTHREAD),y)
CFLAGS+= -pthread -D__IBP_HAVE_PTHREAD=y
LIBS+= -lpthread
endif
ifeq ($(IBP_CUDA),y)
CFLAGS+= -D__IBP_HAVE_CUDA=y
LIBS+= -lcuda -lcudart
endif

# installation configuration variables.
INSTALL=install
PREFIX=/usr/local
BINDIR=$(PREFIX)/bin

# BIN: binary output filenames(s).
BIN=bin/ibp-ng

# SRC_C: basenames of gcc source files.
SRC_C=str value vector intervals trace opts reorder graph assign
SRC_C+= topol-alloc topol-auto topol-add topol
SRC_C+= param-alloc param-add param-get param
SRC_C+= peptide-alloc peptide-residues peptide-atoms peptide-bonds
SRC_C+= peptide-angles peptide-torsions peptide-impropers
SRC_C+= peptide-graph peptide-field
SRC_C+= enum enum-thread enum-reduce enum-write enum-prune
SRC_C+= enum-prune-ddf enum-prune-taf enum-prune-path
SRC_C+= enum-prune-future enum-prune-energy
SRC_C+= dmdgp dmdgp-hash psf

# SRC_N: basenames of nvcc source files.
SRC_N=enum-gpu

# SRC_L: basenames of flex source files.
SRC_L=  assign-scan topol-scan param-scan reorder-scan
SRC_L+= fasta-scan pdb-scan psf-scan cns-scan
SRC_L_C=$(addsuffix .c,$(addprefix src/,$(SRC_L)))

# SRC_Y: basenames of bison source files.
SRC_Y=  assign-parse topol-parse param-parse reorder-parse
SRC_Y+= fasta-parse pdb-parse psf-parse cns-parse
SRC_Y_C=$(addsuffix .c,$(addprefix src/,$(SRC_Y)))
SRC_Y_H=$(addsuffix .h,$(addprefix src/,$(SRC_Y)))

# OBJ: filenames of all compiled object files.
OBJ=  $(addsuffix .o,$(addprefix src/,$(SRC_C)))
OBJ+= $(SRC_L_C:.c=.o) $(SRC_Y_C:.c=.o)

# OBJ: cuda-only compiled object files.
ifeq ($(IBP_CUDA),y)
OBJ+= $(addsuffix .o,$(addprefix src/,$(SRC_N)))
endif

# TESTS_C: basenames of gcc test-case source files.
TESTS_C=base

# TBIN: filenames of all linked test-case binary executables.
TBIN=intervals-alloc intervals-union intervals-intersect intervals-grid
TBIN+= solve-linear solve-spheres solve-omegak
TESTS_O=$(addsuffix .o,$(addprefix tests/,$(TBIN)))
TESTS_X=$(addsuffix .x,$(addprefix tests/,$(TBIN)))

# TOBJ: filenames of all compiled test-case object files.
TOBJ=  $(addsuffix .o,$(addprefix tests/,$(TESTS_C)))

# DATE: date string for making tarballs.
DATE=$(shell date +%Y%m%d)

# SUFFIXES: registered filename suffixes for implicit make rules.
.SUFFIXES: .c .cu .l .o .x .y

# PRECIOUS: files that make is forced to keep around.
.PRECIOUS: %.c %.cu %.h %.o

# all: global, default make target.
all: $(SRC_Y_C) $(SRC_L_C) $(OBJ) $(BIN)

# BIN: binary linkage make target.
$(BIN): $(OBJ) src/ibp-ng.o
	@echo " LD   $@"
	@$(LD) $^ -o $@ $(LIBS)

# .c => .o: gcc source compilation make target.
.c.o:
	@echo " CC   $^"
	@$(CC) $(CFLAGS) -c $^ -o $@

# .cu => .o: nvcc source compilation make target.
.cu.o:
	@echo " NVCC $^"
	@$(NVCC) -c $^ -o $@

# .l => .c: flex compilation make target.
.l.c:
	@echo " LEX  $^"
	@$(LEX) $(LFLAGS) -P$(@:src/%-scan.c=%_io_) -o$@ $^

# .y => .c: bison compilation make target.
.y.c:
	@echo " YACC $^"
	@$(YACC) $(YFLAGS) -p $(@:src/%-parse.c=%_io_) -o $@ $^

# .o => .x: gcc test-case binary linkage make target.
.o.x:
	@echo " LD   $@"
	@$(LD) $(OBJ) $(TOBJ) $^ -o $@ $(LIBS)

# test: target to execute all generated test programs.
test: $(OBJ) $(TOBJ) $(TESTS_X)
	@for tx in $(TESTS_X); do echo " TEST $$tx"; ./$$tx; done

# install: target to install all generated output files.
install: install-bin

# install-bin: target to install all binary files.
install-bin: $(BIN)
	@echo " INSTALL $(BIN)"
	@$(INSTALL) -d $(BINDIR)
	@$(INSTALL) $(BIN) $(BINDIR)

# clean: target to remove all generated intermediate and output files.
clean:
	@echo " CLEAN"
	@rm -f $(SRC_L_C) $(SRC_Y_C) $(SRC_Y_H) $(OBJ) $(TOBJ)
	@rm -f $(TESTS_X) $(TESTS_O) $(BIN)

# again: target to fully recompile all sources and rebuild all binaries.
again: clean all

# target to count lines of all source files.
lines: clean
	@echo " WC"
	@wc -l src/*.[chly] tests/*.[ch]

# target to search all source files for 'fixme' statements.
fixme:
	@echo " FIXME"
	@grep \
	   --recursive --with-filename \
	   --line-number --ignore-case --color \
	   fixme src/*.[chly] tests/*.[ch] || \
	 echo " No statements found"

# target to wrap all source files up into a dated tarball.
dist: clean
	@echo " DIST $(DATE)"
	@rm -rf ibp-ng-$(DATE)
	@rm -f ../ibp-ng-$(DATE).tgz
	@$(INSTALL) -d ibp-ng-$(DATE)
	@cp -ar bin lib src tests Makefile ibp-ng-$(DATE)/
	@tar czf ../ibp-ng-$(DATE).tgz ibp-ng-$(DATE)/
	@rm -rf ibp-ng-$(DATE)

