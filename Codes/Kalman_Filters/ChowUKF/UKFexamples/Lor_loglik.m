function negloglik= Lor_loglik(initPar,y,x0,Px,InfDS,pNoise,oNoise)
    initPar
  oNoise.cov = diag([initPar]);
  InfDS.par = initPar;
    % Full Unscented Kalman Filter step
  % =================================
  [xhat, Px, xh_, Px_, yh_, inov, Py, KG, negloglik,Pxall, Pyall, Px_all]= ukfaddmultiP(x0,Px,pNoise, oNoise, y,[],[],InfDS); 
  
