PhiID -- Integrated Information Decomposition
=============================================

This library contains Octave (v5.0.0) and Matlab (>R2016a) functions to compute
Integrated Information Decomposition (PhiID) in continuous and discrete data.
All relevant details can be found in the original paper:

* Mediano*, Rosas*, Luppi, Carhart-Harris, Bor, Seth and Barrett (2019).
  Towards an extended taxonomy of information dynamics via Integrated
  Information Decomposition. arxiv:2109:13186

If you haven't read this paper, or the Williams & Beer paper below, you are
strongly encouraged to do so before using this software.


Usage
-----

This library implements PhiID in bipartite systems. To compute PhiID on the
time-delayed mutual information of a time series with 2 dimensions and 1000
samples, run:

```octave
X = randn([2, 1000]);
atoms = PhiIDFull(X);
```

The resulting `struct` contains all the PhiID atoms for the input data (see
`help PhiIDFull.m` for details on atom names). In our view, single atoms may
not always be particularly interpretable, but quantities of interest (such as
transfer entropy or integrated information) can be built from them. In general,
some atoms can be negative, but certain combinations of them are guaranteed to
be non-negative, e.g. transfer entropy:

```octave
TE = atoms.xtr + atoms.str + atoms.xty + atoms.sty
```

At the moment, only a few PID functions are implemented: Minimum Mutual
Information (MMI) by Barrett (2015), and Common Change in Surprisal (CCS) by
Ince (2017).

Most of the functions in this library depend on the Java Information Dynamics
Toolbox (JIDT) by Lizier, Ozdemir & Mediano:

https://www.github.com/jlizier/jidt

For convenience, all JIDT-related information and an up-to-date jar are kept in
the `private` folder of this repo.

Feature requests and bug reports are warmly welcome. Email Pedro Mediano (see
email in the paper above) for any questions or comments.

This code is distributed under a BSD-3 licence (i.e. you can use it and modify
it, as long as you acknowledge the original authors, but you cannot distribute
a modified version under the same name).


Further reading
---------------

* Williams PL, Beer RD (2010). Nonnegative decomposition of multivariate
  information. arXiv:1004.2515.

* Barrett AB (2015). Exploration of synergistic and redundant information
  sharing in static and dynamical Gaussian systems. Physical Review E.

* Ince RA (2017). Measuring multivariate redundant information with pointwise
  common change in surprisal. Entropy.


(C) Pedro Mediano, Fernando Rosas and Andrea Luppi, 2019-21

