nIMParton v1.1
-- nIMParton = nuclear (version of) I'mParton

1. Description
This package gives nuclear parton distribution functions (nPDFs) of various
nuclei starting from low Q^2 ~ 0.07 GeV^2, which are based on the analysis
to deep inelastic scattering data applying DGLAP equations with nonlinear
corrections. The extended dynamical parton model is used, which has zero
gluon and zero sea quark distributions at the input scale Q_0^2 = 0.0671 GeV^2.
The parton recombination effect for nuclear shadowing, nucleon swelling
for the EMC effect, and the Fermi motion and the off-shell effect
are taken to model the complicated x-dependence of nuclear modification.
Refs. Rong Wang, Xurong Chen, and Qiang Fu, Nucl. Phys. B 920 (2017) 1 [arXiv:1611.03670v4];
Xurong Chen et al, Int. J. Mod. Phys. E 23 (2014) 1450058 [arXiv:1306.1874v2];

2. Usage
The library consists of a C++ class named nIMParton (read ./nIMParton.h
and ./nIMParton.cpp for details). The constructor function
nIMParton::nIMParton(Z, A) is used to choose a nuclear target.
( eg. nIMParton Calcium(20, 40); ) nIMParton has four important methods:
nMParton::getRToN_p(int Iparton, double x, double Q2),
nMParton::getRToD_p(int Iparton, double x, double Q2),
nMParton::getRToN_n(int Iparton, double x, double Q2),
and nMParton::getRToD_n(int Iparton, double x, double Q2),
which are suggested to be called in users' programs.
RToN_p is parton distribution ratio of the bound proton in nucleus to
the free proton; RToD_p is parton distribution ratio of the bound proton
in nucleus to the bound proton in Deuterium; RToN_n is parton distribution
ratio of the bound neutron in nucleus to the free neutron; RToD_n is parton
distribution ratio of the bound neutron in nucleus to the bound neutron
in Deuterium. Iparton set as -3, -2, -1, 0, 1, 2, 3 correspond to getting
nuclear modification factors R of sbar, dbar, ubar, gluon, u, d, s
distributions respectively. The modification factor R of charm distribution
is not provided, but we suggest using R_c = R_g if it is needed. Because
charm quark distribution is mainly from the gluon splitting. Another important
method of nIMParton is setDataSet(int). setDataSet(1) corresponds to use the
data set A, and setDataSet(2) corresponds to use the data set B. There is a
little difference between set A and set B. (For isospin-scalar nuclei, set A
is suggested, and for isospin non-scalar nuclei, set B is suggested.)


To access the nuclear parton distribution (per nucleon), the formula should be,
f_i^A = ( Z*R_i^{p in A}*f_i^p + N*R_i^{n in A}*f_i^n ) / A,
where i is the index for flavor, f_i^p and f_i^n are the free proton
PDF and the free neutron PDF respectively. Z, N and A are the proton number,
neutron number and the atomic number respectively.

./test.cpp gives an example to get nuclear modifications of different types of
PDFs of ^3Helium target at high Q^2. ./test.cpp can be modified as users' wants.
To run the example,
>tar -zvxf nIMParton.tar.gz
>cd nIMParton
>make
>./test

3. Method
We provide two data sets of nPDFs obtained from the global analyses to
experimental data. Set A is from the global analysis to the data of only
isospin-scalar nuclei, and set B is from analysis to all the measured nuclear
DIS data so far. This package contains a data folder named grid_data which
stores a lot of table grids for various nuclei. The table grids are generated
in the kinematic range of 10^-6 < x < 1 and 0.125 < Q^2 < 2.68435e+08 (GeV^2).
The number of the grid points is 32 for the variable Q^2 and 61 for the
variable x (in version1.1, we add more data points in the large x region).
For nPDF values in small x region (x<10^-4), we use function form A*x^B
to do interpolation. In large x region (x>0.5), we use function form A*(1-x)^B
to do the interpolation. In middle x region (10^-4<x<0.5), we use quadratic
interpolation. The nPDF values outside of the grid range are given using
extrapolation method. The sophisticated extrapolation method is expected
to be effective.

4. Questions
If you have detailed questions concerning these nuclear corrections of PDFs,
or if you find problems/bugs using this package, direct inquires to
wangrong11@mails.ucas.ac.cn (rwangcn8@gmail.com), xchen@impcas.ac.cn,
or wzhu@phy.ecnu.edu.cn.
