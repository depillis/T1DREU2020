%This function is needed for parameter estimation

function negloglik= Lor_loglik(initPar,y,x0,Px,InfDS,pNoise,oNoise)
   % initPar %Printing happens here for the parameter estimates
  initPar
  %oNoise.cov = diag([initPar]);
  oNoise.cov = diag([initPar(4:6)]);
  InfDS.par = initPar;
    % Full Unscented Kalman Filter step
  % =================================
  [xhat, Px, xh_, Px_, yh_, inov, Py, KG, negloglik,Pxall, Pyall, Px_all]= ukfaddmultiP(x0,Px,pNoise, oNoise, y,[],[],InfDS);
  %xhat(:,1)
  
