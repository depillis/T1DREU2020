%--------------------------------------------------------------------------
% N. Tania (May 30, 2018)
%
% ODE Solver

function ydot = modelBasic(t,y,pars)

y1  = y(1);     % M
y2  = y(2);     % TC
y3  = y(3);     % N 
y4  = y(4);     % Tr

% load parameters
sM = 0.001;	rM = pars(1);	KM = pars(2);	deltaM = pars(3);
aNM = pars(4);	bNM = pars(5);	aCM = pars(6);	bCM = pars(7);
aCNM = pars(8);	aMM = pars(9);	bMM = pars(10);	aRM = pars(11);
bRM = pars(12);	rC = pars(13);	KC = pars(14);	deltaC = pars(15);
aMC = pars(16);	bMC = pars(17);	aNC = pars(18);	bNC = pars(19);
sN = pars(20);	rN = pars(21);	KN = pars(22);	deltaN = pars(23);
aCN = pars(24);	bCN = pars(25);	rR = pars(26);	KR = pars(27);
deltaR = pars(28);	aMR = pars(29);	bMR = pars(30);	


% differential equations
dy1 = sM + rM * (1 - y1/KM) * y1 - deltaM * (1 + (aNM * y3/(bNM + y3) + aCM * y2/(bCM + y2) + ...
    (aCNM * y3/(bNM + y3)) * (y2/(bCM + y2))) * ...
    (1 - aMM * y1/(bMM + y1) - aRM * y4/(bRM + y4))) * y1;

dy2 = rC * (1 - y2/KC) * (1 + aMC * y1/(bMC + y1) + aNC * y3/(bNC + y3)) * y2 - deltaC * y2;

dy3 = sN + rN * (1 - y3/KN) * (1 + aCN *y2/(bCN +y2)) * y3 - deltaN * y3;
dy4 = rR * (1 - y4/KR) * (1 + aMR * y1/(bMR + y1)) * y4 -deltaR * y4;

ydot = [dy1; dy2; dy3; dy4];

%%%%%%%%%%%%%%%
% Below is where we check for stopping conditions
MINSTEP = 1e-15;  %Minimum step
persistent tprev elapsedtime
if isempty(tprev)
    tprev = -inf;
end

timestep = t - tprev;
tprev = t;
if (t > 0.01) && (timestep > 0) && (timestep < MINSTEP)    
    error(['Stopped. Time step is too small: ' num2str(timestep)])
end