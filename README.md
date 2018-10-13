# argmin.jl
[![Build Status](https://travis-ci.org/AtsushiSakai/argmin.jl.svg?branch=master)](https://travis-ci.org/AtsushiSakai/argmin.jl)

A numerical solver by pure Julia

Mathematical equations are rendered by [Online LaTeX Equation Editor](https://www.codecogs.com/latex/eqneditor.php)

# API

## solve least square

<img src="https://latex.codecogs.com/gif.latex?\hat{x}=argmin(|Ax-b|^2)" />

## solve_multi_objective_least_square

<img src="https://latex.codecogs.com/gif.latex?\hat{x}=argmin(&space;\lambda&space;_1|Ax&space;=&space;b|^2&plus;\lambda_2|Ax&space;=&space;b|^2&plus;...)" title="\hat{x}=argmin( \lambda _1|Ax = b|^2+\lambda_2|Ax = b|^2+...)" />

## solve_constrained_least_square

<img src="https://latex.codecogs.com/gif.latex?\hat{x}=argmin(|Ax&space;=&space;b|^2)\&space;s.t.\&space;Cx&space;=&space;d" title="\hat{x}=argmin(|Ax = b|^2)\ s.t.\ Cx = d" />

# License

MIT

# Author

- Atsushi Sakai
