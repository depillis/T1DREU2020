clc
% Initial guess
a = 2; b = 3;
k=[a b];
[k lest_squares]=fminsearch('simple_least_squares',k)