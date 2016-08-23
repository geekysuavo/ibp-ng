
# ibp-ng

An implementation (and experimentation sandbox) of the interval
Branch-and-Prune (iBP) structure determination algorithm.

## Introduction

The iBP algorithm enumerates the set of molecular structures that conform
to a user-specified set of distance, angle and dihedral constraints. For
a proof of principle in biomolecules, the reader is referred to:

> Cassioli, et al., _An algorithm to enumerate all possible protein
> conformations verifying a set of distance constraints_,
> BMC Bioinformatics, 2015, 16: 23.

The mathematics used to embed each atom into three-dimensional space are
discussed in:

> Lavor et al., _Clifford algebra and the discretizable molecular distance
> geometry problem_, Advances in Applied Clifford Algebras,
> 2015, 25: 925.

A recursive description of the iBP algorithm is published in:

> Lavor et al., _The interval Branch-and-Prune algorithm for the
> discretizable molecular distance geometry problem with inexact
> distances_, Journal of Global Optimization, 2013, 56: 855.

The current implementation uses both breadth-first and depth-first iterative
tree searches based on a multidimensional index data structure, but this may
change in future versions.

### Compilation

The iBP source code should be compilable on any decently modern
GNU/Linux distribution. The development configuration is:

 * Linux 2.6.18
 * Binutils 2.17.50
 * GCC 4.1.2
 * Flex 2.5.4
 * Bison 2.3

There are two compilation options that must be set in the
[Makefile] (Makefile):

 * **IBP_PTHREAD**: enable support for multiple parallel (CPU) threads.
 * **IBP_CUDA**: enable (experimental, incomplete) support for GPU threads.

Once these options are set, **ibp-ng** may be compiled like so:

```bash
make
```

A simple pull-and-build for the lazy:

```bash
git clone git://github.com/geekysuavo/ibp-ng.git
cd ibp-ng
make
```

### Basic Usage

Successful compilation will produce a file in [bin] (bin) called **ibp-ng**.
For basic usage information, just run:

```bash
bin/ibp-ng -h
```

A toy example is provided in [data/test] (data/test). The script that invokes
**ibp-ng** is named [data/test/run] (data/test/run). It will use four threads
to generate 100 structures using a very fine granularity and only direct
distance feasibility pruning. To execute the script, run:

```bash
cd data/test
./run
```

More examples will be placed in the [data] (data) directory as the source
code progresses.

## Licensing

This project is released under the
[MIT license] (https://opensource.org/licenses/MIT). See the
[LICENSE.md] (LICENSE.md) file for the complete license terms.

And as always, enjoy!
*~ Brad.*

