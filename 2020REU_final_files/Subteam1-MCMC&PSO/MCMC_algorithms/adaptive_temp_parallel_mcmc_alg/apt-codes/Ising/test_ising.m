%Spatial imaging simualtion in the same setting as in he paper
%load data

% test_ising.m
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

load('icefloe.mat');

%parameters of model
global alpha_Ising;
alpha_Ising=1;
global beta_Ising;
beta_Ising=.7;
global IceFloe;

%initial paramters of alghorithm
 param=pt_default_param_ising(10);

%initialise outputs from 3 levels of temperatures
image_out=zeros(40,40,100);
mean_posterior=zeros(40,40);

image_out_mid=zeros(40,40,100);
mean_posterior_mid=zeros(40,40);

image_out_high=zeros(40,40,100);
mean_posterior_high=zeros(40,40);

%main loop
for j=1:100
    [X, beta, stats] = adaptive_pt_ising(param,10^5);
for k=1:40
    for l=1:40
        image_out(k,l,j)=X(k,l,1);
    
        image_out_mid(k,l,j)=X(k,l,4);
     
        image_out_high(k,l,j)=X(k,l,8);
      
    end
end
end

%computing means over replications
for k=1:40
    for l=1:40
        mean_posterior(k,l)=mean(image_out(k,l,:));
        mean_posterior_mid(k,l)=mean(image_out_mid(k,l,:));
        mean_posterior_high(k,l)=mean(image_out_high(k,l,:));
       end
end


%plots
cmap=contrast(mean_posterior,128);
figure('Position',[1 1 400 400])
subplot(2,2,1)
image(IceFloe','CDataMapping','scaled')
axis xy
axis image
colormap(cmap)
subplot(2,2,2)

image(mean_posterior','CDataMapping','scaled')
axis xy
axis image
colormap(cmap)
subplot(2,2,3)
cmap=contrast(mean_posterior_mid);
image(mean_posterior_mid','CDataMapping','scaled')
axis xy
axis image
colormap(cmap)
subplot(2,2,4)
cmap=contrast(mean_posterior_high);
image(mean_posterior_high','CDataMapping','scaled')
axis xy
axis image
colormap(cmap)
