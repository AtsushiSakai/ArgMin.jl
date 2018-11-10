# argmin.jl
[![Build Status](https://travis-ci.org/AtsushiSakai/argmin.jl.svg?branch=master)](https://travis-ci.org/AtsushiSakai/argmin.jl)

A numerical solver by pure Julia

Mathematical equations are rendered by [Online LaTeX Equation Editor](https://www.codecogs.com/latex/eqneditor.php)

# API

## solve_least_square

<img src="https://latex.codecogs.com/gif.latex?\hat{x}=argmin(|Ax-b|^2)" />

or

<img src="https://latex.codecogs.com/gif.latex?\hat{x}=argmin(|Ax&space;=&space;b|^2)\&space;s.t.\&space;Cx&space;=&space;d" title="\hat{x}=argmin(|Ax = b|^2)\ s.t.\ Cx = d" />

## solve_multi_objective_least_square

<img src="https://latex.codecogs.com/gif.latex?\hat{x}=argmin(&space;\lambda&space;_1|Ax&space;=&space;b|^2&plus;\lambda_2|Ax&space;=&space;b|^2&plus;...)" title="\hat{x}=argmin( \lambda _1|Ax = b|^2+\lambda_2|Ax = b|^2+...)" />

## solve_nonlinear_least_square_with_newton_raphson

<img src="https://latex.codecogs.com/gif.latex?\hat{x}&space;=&space;argmin(|f(x)|^2)" title="\hat{x} = argmin(||f(x)||^2)" />

## solve_nonlinear_least_square_with_gauss_newton

<img src="https://latex.codecogs.com/gif.latex?\hat{x}&space;=&space;argmin(|f(x)|^2)" title="\hat{x} = argmin(||f(x)||^2)" />

## solve_nonlinear_least_square_with_levenberg_marquardt

<img src="https://latex.codecogs.com/gif.latex?\hat{x}&space;=&space;argmin(|f(x)|^2)" title="\hat{x} = argmin(||f(x)||^2)" />

## solve_constrained_nonlinear_least_square_with_augmented_lagragian

<img src="https://latex.codecogs.com/gif.latex?\hat{x}&space;=&space;argmin(||f(x)||^2)\&space;s.t.\&space;g(x)=0" title="\hat{x} = argmin(||f(x)||^2)\ s.t.\ g(x)=0" />

## solve_quadratic_programing

<a href="https://www.codecogs.com/eqnedit.php?latex=\hat{x}&space;=&space;argmin(\frac{1}{2}xPx&plus;q'x)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\hat{x}&space;=&space;argmin(\frac{1}{2}xPx&plus;q'x)" title="\hat{x} = argmin(\frac{1}{2}xPx+q'x)" /></a>

or

<a href="https://www.codecogs.com/eqnedit.php?latex=\hat{x}&space;=&space;argmin(\frac{1}{2}xPx&plus;q'x)\&space;s.t.\&space;Ax=b" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\hat{x}&space;=&space;argmin(\frac{1}{2}xPx&plus;q'x)\&space;s.t.\&space;Ax=b" title="\hat{x} = argmin(\frac{1}{2}xPx+q'x)\ s.t.\ Ax=b" /></a>


# License

MIT

# Author

- Atsushi Sakai
