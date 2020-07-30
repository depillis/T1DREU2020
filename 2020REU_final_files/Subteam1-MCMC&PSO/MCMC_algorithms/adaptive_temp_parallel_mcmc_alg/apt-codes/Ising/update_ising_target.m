function p = update_ising_target(x,k,l)
%update log density after  swap the pixel with index (k,l) 

% update_ising_target.m
% 
% Copyright (c) Blazej Miasojedow, Eric Moulines and Matti Vihola 2012
% 
% This file is part of the implementation of the Adaptive Parallel 
% tempering algorithm (hereafter "APT"); see http://arxiv.org/abs/1205.1076.
% 
% APT is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% APT is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with APT.  If not, see <http://www.gnu.org/licenses/>.

global IceFloe;
global alpha_Ising;
global beta_Ising;
p_=0;
p1=0;
p2=0;
if x(k,l) == IceFloe(k,l);
    p_=p_-alpha_Ising;
else
    p_= p_+alpha_Ising;
end
ls=max(1,k-1);
us=min(40,k+1);
lt=max(1,l-1);
ut=min(40,l+1);
for s = ls:us
          
    for t= lt:ut
            if s~=k|t~=l
            if x(k,l) == x(s,t);
                p1=p1+beta_Ising;
            else
                p2=p2+beta_Ising;
            end
            end
    end
end

p1=p1;
p_=p_-2*p1+2*p2;
p=p_;

