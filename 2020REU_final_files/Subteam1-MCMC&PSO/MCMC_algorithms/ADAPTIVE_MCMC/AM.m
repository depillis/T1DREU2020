%   ADAPTIVE METROPOLIS: Generate random samples from a target 
%   distribution with a global covariance matrix V. Function syntax is:
%
%   G = AM(X,OPT)
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
%   s         % Metropolis scaling factor
%   acc_rate  % Metropolis acceptance rate
function G = AM(X,OPT)
%% Initiation
Dims = OPT.Dims;                % Number of dimensions
logpdf = OPT.logpdf;            % Handle to the log-density function
V = OPT.V;                      % Sample covariance matrix
best_X=X;                       % Record highest density state
best_P = logpdf(X);             % Best seen density
P0 = logpdf(X);                 % Last density
eps = -1;                       % Could be -0.5 to -1
pcom = 0;                       % Percentage complete count initialised to zero;
chk_stop=2000;                  % Periodically check to stop burn-in
burnin = 1e6;                   % Initial burnin set to 100% sample
delta = zeros(2,burnin);        % Step change in adaptation
% Take the Cholesky factor of V
[V2,errV]=chol(V,'lower');
if errV
    warning('Bad covariance matrix at start, perform reconditioning..')
    [V3,D] = eig(V);d=diag(D);d(d<=0)=10*realmin;
    V= V3*diag(d)*V3';
    L = tril(V,-1);V=diag(diag(V))+L+L';
    V2=chol(V,'lower');
end
% Define optimal acceptance rate
switch Dims 
    case 1
        optimal=0.441;
    case 2
        optimal=0.352;
    case 3
        optimal=0.316;
    case 4
        optimal=0.285;
    case 5
        optimal=0.275;
    case 6
        optimal=0.273;
    case 7
        optimal=0.270;
    case 8
        optimal=0.268;
    case 9
        optimal=0.267;
    case 10
        optimal=0.266;
    case 11
        optimal=0.265;
    case 12
        optimal=0.264;
    otherwise
        optimal=0.255;
end
% Set initial scaling factor
s = 2.38^2/Dims;
%% Adaptive Metropolis Algorithm.
U = log(rand(1,burnin));        % Random acceptance numbers 
RD = randn(Dims,burnin);        % Random distubance proposal
prerun = 1000;                  % Generate 'prerun' samples before burnin
X_2 = zeros(Dims,prerun);       % Pre-burnin samples
avg_acc = 0;                    % Average acceptance rate
%% Pre-burnin run - generates 'prerun' samples with the initial covariance matrix
mark = cputime;mark_0 = cputime;% Used to keep track of elapsed time
for i = 1:prerun
    if OPT.dsp % Display completion status
        if mod(i,floor(burnin/20))==0
            pcom=pcom+5;
            fprintf('Complete = %d%% burn-in, Time since last update = %f s\n',pcom,cputime-mark); 
            if i>2
            remaining = (20 - i/floor(burnin/20))*(cputime-mark);
            hrs = floor(remaining/3600); remaining = rem(remaining,3600);
            mins = floor(remaining/60); 
            secs = floor(rem(remaining,60));
            fprintf('Predicted burn-in completion in %d hrs:%d mins:%d secs\n',hrs,mins,secs)
            fprintf('Acceptance rate %f %%\n',100*avg_acc);
            end;
        mark=cputime;
        end;
    end; 
    X_new = X + sqrt(s)*V2*RD(:,i); % propose a new sample
    P_new = logpdf(X_new);          % new density
    if P_new>best_P                 % new density is highest 
        best_X=X_new;               % record best seen parameters
        best_P = P_new;
        X = X_new;                  % accept X_new 
        P0 = P_new;
        acc_rate=1;                 % record acceptance rate
    else 
        rho = P_new - P0;
        acc_rate=exp(min(rho,0));   % record acceptance rate
        if U(i)<= rho
        X = X_new;                  % accept X_new 
        P0 = P_new;
        end
    end
    X_2(:,i) = X;                   % Store prerun samples
    %Update scale factor by Robbins-Monro approximation and record everage acceptance rate
    delta(1,i) = (i^eps)*(acc_rate-optimal);
    s = exp(log(s)+delta(1,i));
    avg_acc = avg_acc*(i-1)/i + acc_rate/i;
end
X_bar=mean(X_2,2);          % Sample mean
V=cov(X_2');                % Sample covariance matrix
L = tril(V,-1);V=diag(diag(V))+L+L'; % Check for symmetry
[V2,errV]=chol(V,'lower');
if errV
    warning('Bad covariance matrix at prerun, perform reconditioning..')
    [V3,D] = eig(V);d=diag(D);d(d<=0)=10*realmin;
    V= V3*diag(d)*V3';
    L = tril(V,-1);V=diag(diag(V))+L+L'; % Check for symmetry again
    V2=chol(V,'lower');
end
i=prerun+1;
%% Burnin run
while i <= burnin
    if OPT.dsp % Display completion status
        if mod(i,floor(burnin/20))==0
            pcom=pcom+5;
            fprintf('Complete = %d%% burn-in, Time since last update = %f s\n',pcom,cputime-mark); 
            if i>2
            remaining = (20 - i/floor(burnin/20))*(cputime-mark);
            hrs = floor(remaining/3600); remaining = rem(remaining,3600);
            mins = floor(remaining/60); 
            secs = floor(rem(remaining,60));
            fprintf('Predicted burn-in completion in %d hrs:%d mins:%d secs\n',hrs,mins,secs)
            fprintf('Acceptance rate %f %%\n',100*avg_acc);
            end;
        mark=cputime;
        end;
    end; 
    X_new = X + sqrt(s)*V2*RD(:,i);     % propose a new sample
    P_new = logpdf(X_new);              % new density
    if P_new>best_P                     % new density is highest 
        best_X=X_new;                   % record best seen parameters
        best_P = P_new;
        X = X_new;                      % accept X_new
        P0 = P_new;
        acc_rate=1;                     % record acceptance rate
    else 
        rho = P_new - P0;
        acc_rate=exp(min(rho,0));       % record acceptance rate
        if U(i)<= rho
        X = X_new;                      % accept X_new 
        P0 = P_new;
        end
    end
    % Check stopping rule at every chk_stop counts
    if ~mod(i-prerun,chk_stop)
        if (log(mean(abs(delta(1,i-chk_stop+1:i-chk_stop/2)),2))<log(mean(abs(delta(1,i-chk_stop/2:i-1)),2))) &&...
           (log(mean(abs(delta(2,i-chk_stop+1:i-chk_stop/2)),2))<log(mean(abs(delta(2,i-chk_stop/2:i-1)),2)))
            burnin = i;
            delta = delta(:,1:burnin);   % Adaptation cost
        end
    end
    %Update scale
    gamma = i^eps;
    delta(1,i) = gamma*(acc_rate-optimal);
    s = exp(log(s)+delta(1,i));
    avg_acc = avg_acc*(i-1)/i + acc_rate/i;
    %Update Sigma
    delta_X=(X-X_bar); X_bar = X_bar+gamma*delta_X;
    Sdelta=gamma*(delta_X*delta_X')-((i-1)^eps)*V;
    delta(2,i) = norm(Sdelta,1);
    V = V+Sdelta; L = tril(V,-1);V=diag(diag(V))+L+L'; % This help keeping V symmetrical
    [V2,errV]=chol(V,'lower');
    if errV
        warning('Bad covariance matrix at burnin, perform reconditioning..')
        [V3,D] = eig(V);d=diag(D);d(d<=0)=10*realmin;
        V= V3*diag(d)*V3';
        L = tril(V,-1);V=diag(diag(V))+L+L';
        V2=chol(V,'lower');
    else
        i=i+1;
    end
end
G.cputime.burnin = cputime-mark_0; % Record time taken to adapt
if OPT.dsp % If requested, give feedback on completion status
    figure;plot(1:burnin,delta(2,:));
    title('Adaptation process of the covariance matrix $$\mathbf{\Sigma}$$');set(gca,'YScale','log');
    xlabel('Number of iterations');ylabel('$$\Delta_m = \Vert\Sigma_m - \Sigma_{m-1}\Vert_1$$');
    axis tight;drawnow
    figure;plot(1:burnin,abs(delta(1,:)));
    title('Adaptation process of the scaling factor $$s$$');set(gca,'YScale','log');
    xlabel('Number of iterations');ylabel('$$\Delta_m = \Vert s_m - s_{m-1}\Vert$$');
    axis tight;drawnow
    fprintf('Burn-in completed after %d samples\n',burnin);
    fprintf('Acceptance rate %f %%\n',100*avg_acc);
end
%% Final run
div = OPT.div;                  % Generate 'div' consecutive set of samples
mrun = ceil(OPT.Mmax/div);      % with 'mrun' samples in each set
his = zeros(Dims,1.2*OPT.hgrd+1);   % Histogram grid points
X_bar = zeros(size(X));         % Sample mean vector
pcom = 0;                       % Completion percentage 
mark = cputime;mark_0 = mark;   % Used to keep track of elapsed time
for m=1:div                     % Devide run into 'div' number of pieces to save RAM
    U = log(rand(1,mrun));      % Random acceptance numbers 
    RD = V2*randn(Dims,mrun);   % Random distubance proposal
    TH = zeros(Dims,mrun);      % Final MCMC output
    avg_acc = 0;                % Average acceptance rate
    for i=1:mrun
        X_new = X + sqrt(s)*RD(:,i);    % propose a new sample
        P_new = logpdf(X_new);          % new density
        if P_new>best_P                 % new density is highest 
            best_X=X_new;               % record best seen parameters
            best_P = P_new;
            X = X_new;                  % accept X_new 
            P0 = P_new;
            acc_rate=1;                 % record acceptance rate
        else 
            rho = P_new - P0;
            acc_rate=exp(min(rho,0));   % record acceptance rate
            if U(i)<= rho
            X = X_new;                  % accept X_new 
            P0 = P_new;
            end
        end
        TH(:,i) = X;
        %Update scale
        idx = (i+(m-1)*mrun+burnin);    % Global index
        gamma=idx^eps;
        s =exp(log(s)+gamma*(acc_rate-optimal));
        avg_acc = avg_acc*(i-1)/i + acc_rate/i;
        %Update Sigma
        delta_X=X-X_bar;    X_bar = X_bar+gamma*delta_X;
        Sdelta=gamma*(delta_X*delta_X')-((idx-1)^eps)*V;
        V = V+Sdelta;   L = tril(V,-1);V=diag(diag(V))+L+L';
        [V2,errV]=chol(V,'lower');
        if errV
            warning('Bad covariance matrix at generation state, perform reconditioning..')
            [V3,D] = eig(V);d=diag(D);d(d<=0)=10*realmin;
            V= V3*diag(d)*V3';
            L = tril(V,-1);V=diag(diag(V))+L+L';
            V2=chol(V,'lower');
        end
    end;
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
        fprintf('Acceptance rate %f %%\n',100*avg_acc);
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