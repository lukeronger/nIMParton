nIMParton v1.0

1. Description
This package gives nuclear parton distribution functions (nPDFs) of various
nuclei starting from low Q^2 ~ 0.07 GeV^2, which are based on analysis to
deep inelastic scattering (DIS) data applying DGLAP equations with nonlinear
corrections. The parton recombination effect for nuclear shoadowing, nucleon
swelling for the EMC effect, and the Fermi motion and the off-shell effect
are taken to model the complicated x-dependence of nuclear modification.
Refs. Xurong Chen et al, Int. J. Mod. Phys. E 23 (2014) 1450058
[arXiv:1306.1874v2]; website, http://www.escience.cn/people/ronger/index.html

2. Method
We provide two data sets of nPDFs obtained from the global analyses to
experimental data. Set A is from the global analysis to the data of only
isospin-scalar nuclei, and set B is from analysis to all the measured nuclear
DIS data so far. This package contains a data folder named grid_data which
stores a lot of table grids for various nuclei. The table grids are generated
in the kinematic range of 10^-6 < x < 1 and 0.125 < Q^2 < 2.68435e+08 (GeV^2).
The number of the grid points is 32 for the variable Q^2 and 61 for the
variable x. For nPDF values in small x region (x<10^-4), we use function form
A*x^B to do interpolation. In large x region (x>0.5), we use function form
A*(1-x)^B to do the interpolation. In middle x region (10^-4<x<0.5), we use
quadratic interpolation. The nPDF values outside of the grid range are given
using extrapolation method. The sophisticated extrapolation method is expected
to be effective.

3. Usage
The library consists of a C++ class named nIMParton (look ./nIMParton.h
and ./nIMParton.cpp for details). The construction function
nIMParton::nIMParton(Z, A) is used to choose a nuclear target.
(eg. nIMParton Calcium(20, 40);) nIMParton has two important methods
nIMParton::getRToN(Iparton, X, Q2) and  nIMParton::getRToD(Iparton, X, Q2),
which are suggested to be called in users' programs. RToN is the distibution
ratio to free nucleon, and RToD is the distribution ratio to deuteron. Iparton
set as -3, -2, -1, 0, 1, 2, 3 correspond to getting nuclear modification
factors R of sbar, dbar, ubar, gluon, u, d, s quark/gluon distributions
respectively. The modification factor of charm distribution is not provided,
but we suggest using R_c = R_g if it is needed. Because charm quark
distribution is mainly from the gluon splitting. Another important method of
nIMParton is setDataSet(int). setDataSet(1) corresponds to use the data set A,
and setDataSet(2) corresponds to use the data set B. We also provide methods
getF2RToN(X, Q2) and getF2RToD(X, Q2) to get structure function F2 ratios to
free nucleon and to deuteron respectively.

./test.cpp gives an example to get nuclear modifications of different types of
PDFs of lead target at high Q^2. ./test.cpp can be modified as users' wants.
To run the example,
>tar -zvxf nIMParton-v1.0.tar.gz
>cd nIMParton-v1.0
>make
>./test

4. Questions
If you have detailed questions concerning these nuclear corrections of PDFs,
or if you find problems/bugs using this package, direct inquires to
wzhu@phy.ecnu.edu.cn or rwang@impcas.ac.cn (rwangcn8@gmail.com).

