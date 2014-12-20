%% SLRA_MEX_OBJ - an efficient solver for mosaic-Hankel-like structure
%                 
%  The solver is based on C++ object-oriented implementation.
% 
%  obj = SLRA_MEX_OBJ('new', p, s, r) - creates an SLRA object, based on the 
%  parameters p, s, r described in the documentation of the slra function. 
%  Only mosaic-Hankel-like structure Phi * H and non-zero weights are allowed.
%
%  The created object allows evaluation of the VARPRO cost function f(R),  
%   
%  Internally, the object contains an instance of the a C++ object SLRAObject, 
%  and is created by calling the constructor SLRAObject::SLRAObject.
%
%  SLRA_MEX_OBJ('delete', obj) - deletes the SLRA object obj.
%  Internally, the destructor SLRAObject::~SLRAObject() is called
%
%% Fast cost function and derivatives evaluation
%
%  f = SLRA_MEX_OBJ('func', obj, R) - for a given SLRA object obj,
%  evaluates the VARPRO cost  function at a given (m-r) x m argument R.
%
%  g = SLRA_MEX_OBJ('grad', obj, R) - for a given SLRA object obj,
%  computes the VARPRO cost function matrix gradient at a given  R. 
%  Returns an (m-r) x m matrix g.
%
%% Built-in optimization:
%  [ph, info] = SLRA_MEX_OBJ('optimize', obj, opt) - runs optimization.
%   
%  Input: 
%      opt - optimization options, which include
%        * The options described in the slra.m documentation (Rini, disp)
%        * opt.method - optimization method, for example,
%              'l' - GSL implementation of Levenberg-Marquardt
%              'q' - GSL quasi-Newton methods
%              'n' - GSL Nelder-Mead derivative-free optimization method
%              'p' - own implementation of Levenberg-Marquardt based 
%                    on computing pseudoinverse
%              a complete description of opt.method possible values is
%              contained in 
% 
%        * other optimization options:
%          - advanced options
%              opt.avoid_xi,  opt.ls_correction, opt.reggamma
%          - stopping criteria 
%              opt.epsabs, opt.epsrel, opt.epsgrad, opt.epsx, opt.maxx
%          - method-specific minor parameters
%              opt.step, opt.tol, opt.epscov
%          the complete description and default values are decribed in 
%          the documentation of the OptimizationOptions class.
%
%   Output:
%      ph - the approximating parameter vector (see the documentation in slra.m) 
%      info - additional information on optimization process, which include
%        * the options described in the documentation in slra.m
%        * info.iterinfo - a structure which contains information on each iteration
%  
%% See also
%   slra