function p=log_ising_target(x)

% log_ising_target.m
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

% Evaluate the log-density 
global IceFloe;
global alpha_Ising;
global beta_Ising;
p_=0;
for k= 1:40;
    for l= 1:40;
        if x(k,l) == IceFloe(k,l);
            p_=p_+alpha_Ising;
        end
        for s = max(1,k-1):min(40,k+1)
            for t= max(1,l-1):min(40,l+1)
            if x(k,l) == x(s,t);
                p_=p_+beta_Ising;
            end
            end
        end
       end
end
p=p_-beta_Ising;


