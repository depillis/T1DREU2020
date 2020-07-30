%   FACTOR SLICE SAMPLING: Generate random samples from a target 
%   distribution with a global covariance matrix V. Function syntax is:
%
%   G = FSS(X,OPT)
%  
%   in which:
%
%   X:           The initial state vector
%
%   OPT:         Configurations of the MCMC sampler
%   OPT.Mmax:    The number of samples that will be generated
%   OPT.logpdf:  Handle to the log-density function. 
%   OPT.V:       The initial covariance matrix to approximate the topology
%                of the target distribution
%   OPT.Dims:    The number of dimensions which is equal to the number of
%                elements in X
%              
%   Written by Khoa T. Tran, School of EE & CS
%                            University of Newcastle
%        		             Australia.
%
%   Copyright (C) Khoa T. Tran
%
%   Definitions of more terms in the codes
%   ---------------------------------
%   burnin    % Burn-in time
%   best_X    % Highest joint density state
%   best_P    % Highest joint density
%   X_bar     % Sample mean
%   V         % Sample Covariance Matrix
%   delta     % History of adaptation steps
%   neval     % Density evaluation counter in Slice sampling
%   width     % Initial width in Slice sampling
function G = FSS(X,OPT)
%% Initialization
Dims = OPT.Dims;                % Number of dimensions
logpdf = OPT.logpdf;            % Handle to the log-density function
V = OPT.V;                      % Sample covariance matrix
[V2,D] = eig(V);d=diag(D);      % V2 is the eigen-direction basis
width = 0.1*sqrt(abs(d));       % Initial slice width
eps = -1;                       % Could be -0.5 to -1
opt_eval = 5;                   % Can't be smaller than 5
maxeval = 100*opt_eval;         % Maximum density evaluations per update
chk_stop=1000;                  % Periodically check to stop burn-in
best_X=X;                       % Record highest density state
best_P = logpdf(X);             % Best seen density
P0 = logpdf(X);                 % Last density
prerun = 100;                   % Generate 'prerun' samples before burnin
pcom = 0;                       % Percentage complete count initialised to zero;
burnin = 1e5;                   % Burn-in time
neval = zeros(Dims,burnin);     % Density evaluation counter
delta = zeros(Dims+1,burnin);   % History of adaptation steps   
rnd_vert = log(rand(1,burnin)); % Random vertical position of the slice
rnd_horz = rand(Dims,burnin);   % Random width position
rnd_draw = rand(Dims,burnin);   % Random draw within the slice
X_2 = zeros(Dims,prerun);       % Pre-burnin samples
%% Pre-burnin run - generates 'prerun' samples with the initial covariance matrix
mark = cputime;mark_0 = cputime;% Used to keep track of elapsed time
for i= 1:prerun 
    if OPT.dsp          % Display completion status
        if mod(i,floor(burnin/20))==0
            pcom=pcom+5;
            fprintf('Complete = %d%% burn-in, Time since last update = %f s\n',pcom,cputime-mark); 
            if i>2
                remaining = (20 - i/floor(burnin/20))*(cputime-mark);
                hrs = floor(remaining/3600); remaining = rem(remaining,3600);
                mins = floor(remaining/60); 
                secs = floor(rem(remaining,60));
                fprintf('Predicted completion in %d hrs:%d mins:%d secs\n',hrs,mins,secs)
            end
            mark=cputime;
        end
    end 
    max_try = 100;      % Maximun trials to define the slice
    e=0;                % Counting failures in defining the slice
    z = P0 + rnd_vert(i);      % z = log(pdf(X)*U) is a random 'vertical' level from (0,pdf(X))
    while e<=max_try
        for k=1:Dims    % Perform slice sampling in k directions
            r = width(k)*rnd_horz(k,i); 
            X_lw = X - r*V2(:,k); 
            X_up = X_lw + width(k)*V2(:,k); 
            evals = 1;  % Counting the density evaluations
            % Expansion procedure
            while (logpdf(X_lw)>z) && evals<maxeval
                X_lw = X_lw - (width(k)/10)*V2(:,k);
                evals = evals +1;
            end       
            if evals>=maxeval
                width(k) = width(k)*1.1;
                e=e+1;  % Increase failure counter
                break
            end
            evals = evals +1;
            % Expand in the opposite direction
            while (logpdf(X_up)>z) && evals<maxeval
                X_up = X_up + (width(k)/10)*V2(:,k);
                evals = evals+1;        
            end
            if evals>=maxeval
                width(k) = width(k)*1.1;
                e=e+1;  % Increase failure counter
                break
            end
            %Update the initial width by Robbins-Monro approximation
            delta(k+1,i) = (1/i)*(1-opt_eval/evals);
            width(k) = exp(log(width(k))+delta(k+1,i));
            % Uniform sampling in the defined slice
            X_pk = rnd_draw(k,i)*(X_up-X_lw) + X_lw;
            P0 = logpdf(X_pk);
            evals = evals+1;        
            % shrinking the slice if sample falls outside the true slice
            while(P0<z) && evals<maxeval 
                if sign(X_up - X) == sign(X_pk - X)
                    X_up = X_pk;
                else
                    X_lw = X_pk;
                end
                X_pk = rand*(X_up-X_lw) + X_lw; % draw again
                P0 = logpdf(X_pk);
                evals = evals+1;
            end
            neval(k,i) = evals;     
            if evals>=maxeval
                 width(k) = width(k)*0.9;
                 e=e+1;
                 break
            end
            X = X_pk;       % New sample found
            if P0>best_P    % Compare with the highest density and update if new sample is better
                best_P= P0;
                best_X = X; 
            end
        end
        if k == Dims        % Finish this iteration
            max_try=e-1; 
        else                % Repeat this iteration
            rnd_draw(:,i) = rand(Dims,1);
            rnd_horz(:,i) = rand(Dims,1);
        end
    end
    X_2(:,i) = X;           % Store prerun samples
end
X_bar=mean(X_2,2);          % Sample mean
V=cov(X_2');                % Sample covariance matrix
L = tril(V,-1);V=diag(diag(V))+L+L'; % Check for symmetry
[V2,~] = eig(V);            % Eigendecomposition
chk_flg = logical([0 ones(1,Dims)]); % Adaptation flags for Robbins-Monro approximation
i=prerun+1;                 % Keep the iteration counter for next while loop
%% Burnin run
while i <= burnin 
    if OPT.dsp          % Display completion status
        if mod(i,floor(burnin/20))==0
            pcom=pcom+5;
            fprintf('Complete = %d%% burn-in, Time since last update = %f s\n',pcom,cputime-mark); 
            if i>2
                remaining = (20 - i/floor(burnin/20))*(cputime-mark);
                hrs = floor(remaining/3600); remaining = rem(remaining,3600);
                mins = floor(remaining/60); 
                secs = floor(rem(remaining,60));
                fprintf('Predicted completion in %d hrs:%d mins:%d secs\n',hrs,mins,secs)
            end
            mark=cputime;
        end
    end 
    max_try = 100;      % Maximun trials to define the slice
    e=0;                % Counting failures in defining the slice
    z = P0 + rnd_vert(i);      % z = log(pdf(X)*U) is a random 'vertical' level from (0,pdf(X))
    while e<=max_try
        for k=1:Dims    % Perform slice sampling in k directions
            r = width(k)*rnd_horz(k,i); 
            X_lw = X - r*V2(:,k); 
            X_up = X_lw + width(k)*V2(:,k); 
            evals = 1;  % Counting the density evaluations
            % Expansion procedure
            while (logpdf(X_lw)>z) && evals<maxeval
                X_lw = X_lw - (width(k)/10)*V2(:,k);
                evals = evals +1;
            end       
            if evals>=maxeval
                width(k) = width(k)*1.1;
                e=e+1;  % Increase failure counter
                break
            end
            evals = evals +1;
            % Expand in the opposite direction
            while (logpdf(X_up)>z) && evals<maxeval
                X_up = X_up + (width(k)/10)*V2(:,k);
                evals = evals+1;        
            end
            if evals>=maxeval 
                width(k) = width(k)*1.1;
                e=e+1;  % Increase failure counter
                break
            end
            %Update the initial width by Robbins-Monro approximation
            delta(k+1,i) = (1/i)*(1-opt_eval/evals);
            width(k) = exp(log(width(k))+delta(k+1,i));
            % Uniform sampling in the defined slice
            X_pk = rnd_draw(k,i)*(X_up-X_lw) + X_lw;
            P0 = logpdf(X_pk);
            evals = evals+1;        
            % shrinking the slice if sample falls outside the true slice
            while(P0<z) && evals<maxeval 
                if sign(X_up - X) == sign(X_pk - X)
                    X_up = X_pk;
                else
                    X_lw = X_pk;
                end
                X_pk = rand*(X_up-X_lw) + X_lw; % draw again
                P0 = logpdf(X_pk);
                evals = evals+1;
            end
            neval(k,i) = evals;     
            if evals>=maxeval
                 width(k) = width(k)*0.9;
                 e=e+1;
                 break
            end
            X = X_pk;       % New sample found
            if P0>best_P    % Compare with the highest density and update if new sample is better
                best_P= P0;
                best_X = X; 
            end
        end
        if k == Dims        % Finish this iteration
            max_try=e-1; 
        else                % Repeat this iteration
            rnd_draw(:,i) = rand(Dims,1);
            rnd_horz(:,i) = rand(Dims,1);
        end
    end
    %Update Sample Covariance Matrix
    gamma = i^eps;      % Adaptation step size
    delta_X=(X-X_bar); X_bar = X_bar+gamma*delta_X;
    Sdelta=gamma*(delta_X*delta_X')-((i-1)^eps)*V;
    V = V+Sdelta;L = tril(V,-1);V=diag(diag(V))+L+L';
    delta(1,i) = norm(Sdelta,1);
    [V2,~] = eig(V);
    % Check stopping rule at every chk_stop counts
    if ~mod(i-prerun,chk_stop)
        Flags=chk_flg;
        chk_flg(Flags) = ~(log(mean(abs(delta(Flags,i-chk_stop+1:i-chk_stop/2)),2))<log(mean(abs(delta(Flags,i-chk_stop/2:i-1)),2)));
        if all(~chk_flg) && (log(mean(abs(delta(1,i-chk_stop+1:i-chk_stop/2)),2))<log(mean(abs(delta(1,i-chk_stop/2:i-1)),2)))
            burnin = i;
            delta = delta(:,1:burnin);   % Adaptation delta
            neval = neval(:,1:burnin);
        end
    end
    i = i+1;
end
G.cputime.burnin = cputime-mark_0; % Record burnin time
if OPT.dsp % Display burnin results
    figure;plot( 1:burnin,delta(1,1:burnin),'.k');
    title('Adaptation process of the covariance matrix $$\mathbf{\Sigma}$$');set(gca,'YScale','log');
    xlabel('Number of iterations');ylabel('$$\Delta_m = \Vert\Sigma_m - \Sigma_{m-1}\Vert_1$$');
    axis tight;drawnow
    figure;plot( 1:burnin,abs(delta(2,1:burnin)),'.k');
    title('Adaptation process of the initial width $$\mathbf{w}$$');set(gca,'YScale','log');
    xlabel('Number of iterations');ylabel('$$\Delta_m = \Vert w_m - w_{m-1}\Vert$$');
    axis tight;drawnow
    fprintf('Burn-in completed after %d samples\n',burnin);
    fprintf('Each sample costs %d density evaluations on average\n and %d evaluations maximum\n'...
        ,(mean(neval(length(neval(:))-chk_stop*Dims:end)))*Dims,(max(neval(:)))*Dims);
    fprintf('Optimum width component#%d is %d\n',[1:Dims; width']);
end
%% Final run
div = OPT.div;                  % Generate 'div' consecutive set of samples
mrun = ceil(OPT.Mmax/div);      % with 'mrun' samples in each set
his = zeros(Dims,1.2*OPT.hgrd+1);       % Histogram grid points
step_width = bsxfun(@times,V2,width');  % Step width in eigendirections
X_bar = zeros(size(X));         % Sample mean vector
pcom = 0;                       % Completion percentage 
mark = cputime;mark_0 = mark;   % Used to keep track of elapsed time
for m=1:div                     % Devide run into 'div' number of pieces to save RAM
    TH = zeros(Dims,mrun);      % Store all MCMC sample in TH
    rnd_vert = log(rand(1,mrun));     % Random vertical position of the slice
    rnd_horz = rand(Dims,mrun);       % Random width position
    rnd_draw = rand(Dims,mrun);       % Random draw within the slice
    neval = zeros(Dims,1);            % Density evaluation counter
    for i = 1:mrun
        z = P0 + rnd_vert(i);         % z = log(pdf(X)*U) is a random 'vertical' level from (0,pdf(X))
        for k=1:Dims    % Perform slice sampling in k directions
            r = width(k)*rnd_horz(k,i);
            X_lw = X - r*V2(:,k); 
            X_up = X_lw + step_width(:,k); 
            evals = 1;  % Counting the density evaluations
            % Expansion procedure
            while (logpdf(X_lw)>z) && evals<maxeval
                X_lw = X_lw - width(k)*V2(:,k);
                evals = evals +1;
            end       
            if evals>=maxeval 
                error('It takes too many iterations to step out, expanded width at iter#%d, dimension#%d',i,k)
            end
            % Expand in the opposite direction
            evals = evals + 1;
            while (logpdf(X_up)>z) && evals<maxeval
                X_up = X_up + width(k)*V2(:,k);
                evals = evals+1;        
            end
            if evals>=maxeval
                error('It takes too many iterations to step out, expanded width at iter#%d, dimension#%d',i,k) 
            end
            % Uniform sampling in the defined slice
            X_pk = rnd_draw(k,i)*(X_up-X_lw) + X_lw;
            P0 = logpdf(X_pk);
            evals = evals + 1;
            % shrinking the slice if sample falls outside the true slice
            while(P0<z) && evals<maxeval 
                if sign(X_up - X) == sign(X_pk - X)
                    X_up = X_pk;
                else
                    X_lw = X_pk;
                end
                X_pk = rand*(X_up-X_lw) + X_lw; % draw again
                P0 = logpdf(X_pk);
                evals = evals+1;
            end
            neval(k) = neval(k)*(i-1)/i + evals/i; 
            if evals>=maxeval 
                 error('It takes too many iterations to shrink in, stopped at iter#%d, dimension#%d',i,k) 
            end
            X = X_pk;       % New sample found
            if P0>best_P    % Compare with the highest density and update if new sample is better
                best_P= P0;
                best_X = X; 
            end
        end
        TH(:,i) = X;
    end
    %% Produce sample histograms
    % Define grid points
    h_max = max(TH,[],2);    h_min = min(TH,[],2);    
    if m==1 % First grid
        range = 1.2;
        h_range = h_max-h_min;
        hist_grid = my_linspace(h_min-0.1*h_range,h_max+0.1*h_range,range*OPT.hgrd+1);
        mark_low = min(hist_grid,[],3); mark_high = max(hist_grid,[],3); 
    else    % Subsequence grid
        while any(h_max>mark_high) || any(h_min<mark_low)
            % Expand grid points
            range = range + 0.2;
            hist_grid = my_linspace(mark_low-0.1*h_range,mark_high+0.1*h_range,range*OPT.hgrd+1);
            mark_low = min(hist_grid,[],3); mark_high = max(hist_grid,[],3); 
            % Added zeros into old histogram
            his = [zeros(Dims,0.1*OPT.hgrd) his zeros(Dims,0.1*OPT.hgrd)];
        end
    end
    %% Compute raw sample histogram
    for k=1:Dims
        [memo,~] = hist(TH(k,:),hist_grid(k,1,:)); 
        his(k,:)= his(k,:)*(m-1)/m+memo/m;
    end
    %% Compute sample mean
    X_bar = X_bar*(m-1)/m+mean(TH,2)/m;
    %% Save MCMC samples to hard disk and free RAM
    save([OPT.filename '_' num2str(m)],'TH','-v7.3')
    if OPT.dsp  % Display completion status
        pcom=pcom+100/div;
        fprintf('Complete = %d%%, Time since last update = %f s\n',pcom,cputime-mark); 
        remaining = (div - m)*(cputime-mark);
        hrs = floor(remaining/3600); remaining = rem(remaining,3600);
        mins = floor(remaining/60); 
        secs = floor(rem(remaining,60));
        fprintf('Predicted completion in %d hrs:%d mins:%d secs\n',hrs,mins,secs)
        fprintf('Each sample costs %d density evaluations on average\n'...
        ,(mean(neval(:)))*Dims);
        mark=cputime;
    end
end
G.cputime.run = cputime-mark_0;         % Record time taken to run MCMC
G.best_TH = best_X;                     % Record highest density state
G.mean = X_bar;                         % Sample mean state
G.filename = OPT.filename;              % Passing back the filename
% Normalise sample histograms in each dimension
for k=1:Dims
    h.p = his(k,:);h.x = hist_grid(k,1,:);h.x = h.x(:)';
    h.p=h.p/(sum(h.p)*(h.x(2)-h.x(1)));
    G.hist.p(k,:) = h.p; G.hist.x(k,:) = h.x;
end
if OPT.dsp  % Display final results
    fprintf('Successfully generated %d samples in %d seconds\n',div*mrun, G.cputime.run);
end