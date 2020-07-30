function negloglik= Lor_loglikMultiP2(initPar,y,x0,Px,InfDS,pNoise,oNoise)
 
  initPar
  oNoise.cov = diag([initPar]);
  InfDS.par = initPar;
 %if (initPar(3) >60 ) error('It is probably time to stop.'); end

  
    % Full Unscented Kalman Filter step
  % =================================
[xhat, Px, xh_, Px_, yh_, inov, Py, KG, negloglik,Pxall]= ukfaddmultiP(x0,Px,pNoise, oNoise, y,[],[],InfDS); 


