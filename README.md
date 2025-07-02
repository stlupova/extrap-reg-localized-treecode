# extrap-reg-localized-treecode
Code for https://arxiv.org/abs/2507.00156

For questions/comments contact Svetlana Tlupova (svetlana.tlupova@farmingdale.edu) 

Reference:
J. Siebor and S. Tlupova. Localized evaluation and fast summation in the extrapolated regularization method for integrals in Stokes flow. arXiv; Cornell University Library, 2025. https://arxiv.org/abs/2507.00156

This code was developed as part of work supported by the National Science Foundation grant DMS-2012371.

## Testing instructions

### Test: Stokeslet integral past a translating spheroid:

1.	Code subdirectory: `extrap_tree_stokeslet`
2.	In utilities.h, set the following parameters:
    * grid size (h=1/N_gridlines): N_gridlines = 32, 64, 128, 256, etc.
    * type of regularization del_h:
        * del_h = 1:  delta=c*h
        * del_h = 2 (anything not equal 1):  delta=c*h^(⅘)
3.  In KITC.h, set the following treecode parameters:
    * degree of interpolating polynominal P (recommended: P=6 for h=1/64, P=8 for h=1/128, P=10 for h=1/256, etc.)
    * leaf size N0 (recommended: N0=2000 for h=1/64, N0=4000 for h=1/128, N0=8000 for h=1/256, etc.)
    * MAC parameter theta in sq_theta=theta^2 (recommended: theta=0.6)
5.  Compile using `make`. The makefile is set up for a Mac OS.
6.  Run using `./main.out`
   
