SLRA: package for structured low-rank approximation 
==============================================================================
Overview
--------

This Matlab/Octave/R package is designed for solving
structured low-rank approximation problems 
\f[
\text{minimize over} \quad \widehat p \quad \|p - \widehat p\| \quad
\text{subject to} \quad \text{rank}({\mathscr{S}}(\widehat p)) \le r,
\f]
where \f${\mathscr{S}}: \mathbb{R}^{n_p} \to \mathbb{R}^{m \times n}\f$
is an *affine matrix structure*,  \f$p \in \mathbb{R}^{n_p}\f$
is a given parameter vector, and *r* is a specified bound on the rank.

### Considered structures
 * A *general affine structure* 
\f$\mathscr{S}(p) = S_0 + \sum\limits_{k=1}^{n_p} p_k S_k \f$
where \f$S_k\f$ are matrices of zeros and ones.
 * A special case of *mosaic-Hankel-like structure*
\f$\mathscr{S}(p) = \Phi \mathscr{H}_{\bf m, \bf n}(p) \f$,
where \f$\Phi\f$ is full row rank and 
\f$\mathscr{H}_{\bf m, \bfn}\f$ is a mosaic Hankel structure.

### Approximation criteria
 * General weighted semi-norm 
\f$\|p\|^2_W := p^{\top} W p\f$
defined by a positive semidefinite matrix \fW \in \f
 * Elementwise extended weighted 2-norm
\f$\|p\|^2_w := \sum\limits_{k=1}^{n_p} w_k p_k^2\f$,
defined by a weight vector \f$w \in [0,\infty]^{n_p}\f$, where
    * \f$w_k = \infty\f$ is equivalent to the constraint
      \f$\widehat{p}_k = p_k\f$ (fixed values)
    * \f$w_k = 0\f$ implies that 
      \f$\widehat{p}_k, p_k\f$ are not used (missing values)

### Constraint on the kernel of approximating matrix
The left kernel of the approximating matrix is determined by the matrix
\f$R \in \mathbb{R}^{(m-r)\times m}\f$ such that \f$R \mathscr{S}(\widehat{p}) = 0 \f$.
Additional constraints can be imposed on the matrix $R$:
  * *General linear constraint* \f$R = \text{vec}_d^{-1} \theta^{\top} \Psi\f$,
	where \f$\Psi \in \mathbb{R}^{n_{\theta}\times md}\f$, and \f$\theta \in n_{\theta}\f$.
  * *Matrix-product linear constraint*  \f$R = \Theta \Psi\f$,  where 
    \f$\Theta in \mathbb{R}^{d \times m''}\f$ and 
    \f$\Psi \in \bbR^{m'' \times m}\f$ is a full row rank matrix.
	The matrix-product linear constraint is a special case of the general linear constraint
    since \f$ \text{vec}^{\top} (\Theta \Psi) = \text{vec}^{\top}(\Theta) (\Psi \otimes I)\f$

The problem formulation and definitions of the objects can be found
in \cite slra-software. (Note that in \cite slra-software, the general affine structure 
and general weight matrix are not considered. The definitions for th general
affine structure and weign matrix can be found in \cite slra-efficient or \cite slra-ext.)
     
Algorithms
----------
This package contains implementation of the following methods for 
 1. fast C++ implementation of the variable projection (VARPRO)
    method for mosaic Hankel matrices \cite slra-efficient
 2. an implementation of the VARPRO method for SLRA with 
    missing data \cite slra-ext.
    This method is also called ''experimental Matlab solver'' in \cite slra-ext.
 3. an implementation of the factorization approach to SLRA
    based on a penalty method \cite rslra

The following table contains a summary of the features 

  Feature                                   | 1.   | 2.   | 3.
--------------------------------------------|------|------|------ 
  General affine structure                  | -    | +    | +  
  mosaic-Hankel-like structure              | ++   | +    | +
  weight matrix *W*                         | -    | +    | + 
  elementwise weights \f$w_k =(0,\infty]\f$ | +    | +    | +      
  missing data  \f$w_k = 0\f$               | -/+  | +    | +
  General linear constraint on the kernel   | -    | +    | + 
  matrix-product constraint on the kernel   | +    | +    | +
  MATLAB implementation/interface           | +    | +    | +
  R      implementation/interface           | +    | -    | -
 
Note: in 1., missing data can be approximated by small weights

Using the package
-----------------

The package consists of a single \ref slra_doc function.
The standard help for the function is also available by typing:
* `help slra` in MATLAB/Octave
* `?slra` in R
In the HTML version of the documentation, the help is available only for MATLAB.



Directories `test_m` and `test_r` contain demo files for MATLAB/Octave and R.


### Citing the package
If you use the package in your research, please cite \cite slra-software :

    @Article{slra-software,
        author = {I. Markovsky and K. Usevich},
    	title = {Software for weighted structured low-rank approximation},
    	journal = {J. Comput. Appl. Math.},
        volume = {256},
    	pages = {278--292},
    	year = {2014},
    }
	
	

