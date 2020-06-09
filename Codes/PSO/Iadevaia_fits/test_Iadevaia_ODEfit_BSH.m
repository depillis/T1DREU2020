
% %Author the ODE solver
%function test_Iadevaia_ODEfit
clear all;
close all;

N=65; %Number of reactions
t0= 0;
scale_factor=1;
tfinal = scale_factor*60;
global stimulated;

t_interpolate= linspace(0,tfinal,15);


%Solve the ODE, with non-stimulated conditions. 


    stimulated=0;
    
    reduced_model_parameters
   
   [t1,y1] = ode15s(@IGFRODE,[0 tfinal], init_vec);
    
   [t1_perturb,y1_perturb] = ode15s(@IGFRODE,t_interpolate, init_vec);
   
   y1_perturb = y1_perturb + 0.005*randn(size(y1_perturb));


%Solve the ODE, with stimulated conditions of IGF-1

stimulated=1;

reduced_model_parameters;

[t2,y2] = ode15s(@IGFRODE,[0 tfinal],init_vec);

[t2_perturb,y2_perturb] = ode15s(@IGFRODE,t_interpolate, init_vec);

y2_perturb = y2_perturb + 0.005*randn(size(y2_perturb));


% t_interpolate=t_interpolate./scale_factor;
% t=t./scale_factor;


 %We need to plot:
        %p-AKT-x(22)
        %p-GSK3-x(28)
        %p-MAPK-x(19)
        %p-mTOR-x(35)
        %p-P70S6K-x(37)
        %p-TSC2-x(31)
 ipAKT=22;
 ipGSK3 =28;
 ipMAPK=19;
 ipmTOR=35;
 ipP70S6K=37;
 ipTSC2=31;
 
 %Plot the starved results
 y=y1;
 t=t1;
 yp=y1_perturb;
 tp=t1_perturb;
 
 figure(1)
subplot(2,3,1)
plot(t,y(:,ipAKT),'-r','Linewidth',3)
hold on;
plot(tp,yp(:,ipAKT),'*b','Linewidth',3)
hold off;
title('p-AKT:starved');
xlabel('Time t');
ylabel('Solution y');
legend('p-AKT')

subplot(2,3,2)
plot(t,y(:,ipGSK3),'-r','Linewidth',3)
hold on;
plot(tp,yp(:,ipGSK3),'*b','Linewidth',3)
hold off;
title('p-GSK3: starved');
xlabel('Time t');
ylabel('Solution y');
legend('p-GSK3')

subplot(2,3,3)
plot(t,y(:,ipMAPK),'-r','Linewidth',3) 
hold on
plot(tp,yp(:,ipMAPK),'*b','Linewidth',3) 
hold off;
title('p-MAPK:starved');
xlabel('Time t');
ylabel('Solution y');
legend('p-MAPK')

subplot(2,3,4)
plot(t,y(:,ipmTOR),'-r','Linewidth',3) 
hold on;
plot(tp,yp(:,ipmTOR),'*b','Linewidth',3)
hold off;
hold off;
title('p-mTOR:starved');
xlabel('Time t');
ylabel('Solution y');
legend('p-mTOR')

subplot(2,3,5)
plot(t,y(:,ipP70S6K),'-r','Linewidth',3)
hold on;
plot(tp,yp(:,ipP70S6K),'*b','Linewidth',3)
hold off;
title('p-p7SK06:starved');
xlabel('Time t');
ylabel('Solution y');
legend('p-p7SK06')

subplot(2,3,6)
plot(t,y(:,ipTSC2),'-r','Linewidth',3)
hold on;
plot(tp,yp(:,ipTSC2),'*b','Linewidth',3)
hold off;
title('p-TSC2:starved');
xlabel('Time t');
ylabel('Solution y');
legend('p-TSC2')

%Plot the stimulated results
 y=y2;
 t=t2;
 yp=y2_perturb;
 tp=t2_perturb;
 
figure(2)
subplot(2,3,1)
plot(t,y(:,ipAKT),'-r','Linewidth',3)
hold on;
plot(tp,yp(:,ipAKT),'*b','Linewidth',3)
hold off;
title('p-AKT:starved');
xlabel('Time t');
ylabel('Solution y');
legend('p-AKT','in silico p-AKT')

subplot(2,3,2)
plot(t,y(:,ipGSK3),'-r','Linewidth',3)
hold on;
plot(tp,yp(:,ipGSK3),'*b','Linewidth',3)
hold off;
title('p-GSK3: starved');
xlabel('Time t');
ylabel('Solution y');
legend('p-GSK3','in silico p-GSK3')

subplot(2,3,3)
plot(t,y(:,ipMAPK),'-r','Linewidth',3) 
hold on
plot(tp,yp(:,ipMAPK),'*b','Linewidth',3) 
hold off;
title('p-MAPK:starved');
xlabel('Time t');
ylabel('Solution y');
legend('p-MAPK','in silico p-MAPK')

subplot(2,3,4)
plot(t,y(:,ipmTOR),'-r','Linewidth',3) 
hold on;
plot(tp,yp(:,ipmTOR),'*b','Linewidth',3)
hold off;
hold off;
title('p-mTOR:starved');
xlabel('Time t');
ylabel('Solution y');
legend('p-mTOR','in silico p-mTOR')

subplot(2,3,5)
plot(t,y(:,ipP70S6K),'-r','Linewidth',3)
hold on;
plot(tp,yp(:,ipP70S6K),'*b','Linewidth',3)
hold off;
title('p-p7SK06:starved');
xlabel('Time t');
ylabel('Solution y');
legend('p-p7SK06','in silico p-p7SK06')

subplot(2,3,6)
plot(t,y(:,ipTSC2),'-r','Linewidth',3)
hold on;
plot(tp,yp(:,ipTSC2),'*b','Linewidth',3)
hold off;
title('p-TSC2:starved');
xlabel('Time t');
ylabel('Solution y');
legend('p-TSC2','in silico p-TSC2')



    function dxdt = IGFRODE(t,x)
        
        dxdt = zeros (65,1);
       
       reduced_model_parameters;

        %Units for the ODE system go as follow:
        %Molecule weight: moles scaled by -9, micromoles
        %Volume: litre scaled by -6, microlitre
        %Time: 60 second intervals
        %Model reactions:
        
       
        
        %IGF-1 reactions
        
        dxdt(1) = -k1r1*x(1)*x(2) + k2r1*x(3); % IGF 
        
        dxdt(2) = -k1r1*x(1)*x(2) + k2r1*x(3);  %IGFR 
        
        dxdt(3) = k1r1*x(1)*x(2) - k2r1*x(3)  - k1r2*x(3)  + k2r2*x(4)...
                 -k1r3*x(5)*x(3) + k2r3*x(40) + k1r4*x(40) - k1r11*x(8)*x(3)...
                 +k2r11*x(43)    + k1r12*x(43); %IGFR*
        
        dxdt(4) = k1r2*x(3) - k2r2*x(4);%IGFRi
        
        dxdt(5) = -k1r3*x(5)*x(3) + k2r3*x(40) + k1r5*x(6) - k1r6*x(5)*x(19)...
                  +k2r6*x(41)- k1r8*x(5)*x(37)+k2r8*x(42)+k1r10*x(7);%IRS1
              
        dxdt(6) = k1r4*x(40) - k1r5*x(6) - k1r13*x(8)*x(6) + k2r13*x(44)...
                  +k1r14*x(44) - k1r32*x(21)*x(6) + k2r32*x(52)+ k1r33*x(52);%IRS1*
              
        dxdt(7) = k1r7*x(41)+k1r9*x(42)-k1r10*x(7);%IRS1(u)
        
        dxdt(8) = -k1r11*x(8)*x(3) + k2r11*x(43) - k1r13*x(8)*x(6)...
                  + k2r13*x(44) + k1r15*x(9);%RasGDP
        
        dxdt(9) = k1r12*x(43) + k1r14*x(44) - k1r15*x(9)...
                  - k1r23*x(9)*x(11) + k2r23*x(15);%RasGTP
        
        dxdt(10) = -k1r16*x(10)*x(13) + k2r16*x(45) + k1r19*x(46)...
                   -k1r20*x(10)*x(22) + k2r20*x(47) + k1r22*x(12);%Raf
        
        dxdt(11) = k1r17*x(45) - k1r18*x(11)*x(14) + k2r18*x(46) ...
                  -k1r23*x(9)*x(11) + k2r23*x(15);%Raf*
        
        dxdt(12) = k1r21*x(47) - k1r22*x(12);%Raf(u)
        
        dxdt(13) = -k1r16*x(10)*x(13) + k2r16*x(45) + k1r17*x(45);%PKC
        
        dxdt(14) = -k1r18*x(11)*x(14)+k2r18*x(46)+k1r19*x(46)...
                   -k1r26*x(17)*x(14)+k2r26*x(49)+k1r27*x(49);%PP2A
        
        dxdt(15) = k1r23*x(9)*x(11)-k2r23*x(15)-k1r24*x(16)*x(15)...
                   +k2r24*x(48) +k1r25*x(48);%RasRaf
               
        dxdt(16) = -k1r24*x(16)*x(15) -k2r24*x(48) +k1r27*x(49);%MEK
        
        dxdt(17) = k1r25*x(48) -k1r26*x(17)*x(14) +k2r26*x(49)...
                  -k1r28*x(18)*x(17) + k2r28*x(50) +k1r29*x(50);%MEK*
        
        dxdt(18) = -k1r28*x(18)*x(17) +k2r28*x(50) +k1r31*x(51);%MAPK
        
        dxdt(19) = -k1r6*x(5)*x(19) +k2r26*x(41) +k1r7*x(41) +k1r29*x(50)...
                   -k1r30*x(19)*x(20) +k2r30*x(51) -k1r48*x(30)*x(19)...
                   +k2r48*x(58) +k1r49*x(58) -k1r57*x(36)*x(19)...
                   +k2r57*x(61) +k1r58*x(61);%MAPK*
        
        dxdt(20) = -k1r30*x(19)*x(20) +k2r30*x(51) +k1r31*x(51);%MPK1
        
        dxdt(21) = -k1r32*x(21)*x(6) +k2r32*x(52) +k1r34*x(22)...
                   -k1r35*x(21)*x(24) +k2r35*x(53) +k1r37*x(23);%AKT
        
%         dxdt(22) = -k1r20*x(10)*x(22) +k2r20*x(47) +k1r21*x(47)...
%                    +k1r33*x(52) -k1r34*x(22) -k1r38*x(25)*x(22)...
%                    +k2r38*x(54) +k1r39*x(54) -k1r41*x(27)*x(22)...
%                    +k2r41*x(55) +k1r42*x(55) +k1r50*x(30)*x(22)...
%                    +k2r50*x(59) +k1r51*x(59) -k1r65*x(38)*x(22)...
%                    +k2r65*x(65) +k1r66*x(65);%AKT*
               
      dxdt(22) = -k1r20*x(10)*x(22) +k2r20*x(47) +k1r21*x(47)...
                   +k1r33*x(52) -k1r34*x(22) -k1r38*x(25)*x(22)...
                   +k2r38*x(54) +k1r39*x(54) -k1r41*x(27)*x(22)...
                   +k2r41*x(55) +k1r42*x(55) -k1r50*x(30)*x(22)...
                   +k2r50*x(59) +k1r51*x(59) -k1r65*x(38)*x(22)...
                   +k2r65*x(65) +k1r66*x(65);%AKT*         
        
        dxdt(23) = k1r36*x(53) -k1r37*x(23);%AKT(u)
        
        dxdt(24) = -k1r35*x(21)*x(24) +k2r35*x(53) +k1r36*x(53);%PTEN
        
        dxdt(25) = -k1r38*x(25)*x(22) +k2r38*x(54) +k1r40*x(26)...
                   -k1r44*x(29)*x(25) +k2r44*x(56) +k1r45*x(56);%AMPK
        
        dxdt(26) = k1r39*x(54)-k1r40*x(26);%AMPK(p)
        
        dxdt(27) = -k1r41*x(27)*x(22) +k2r41*x(55) +k1r43*x(28)...
                   -k1r46*x(29)*x(27) +k2r46*x(57) +k1r47*x(57);%GSK3
        
        dxdt(28) = k1r42*x(55)-k1r43*x(28);%GSK3(p)
        
        dxdt(29) = -k1r44*x(29)*x(25) +k2r44*x(56) -k1r46*x(29)*x(27)...
                   +k2r46*x(57) +k1r49*x(58);%TSC2
        
        dxdt(30) = k1r45*x(56) +k1r47*x(57) -k1r48*x(30)*x(19)+k2r48*x(58)...
                  -k1r50*x(30)*x(22) +k2r50*x(59) +k1r52*x(31)...
                  -k1r53*x(32)*x(30) +k2r53*x(60) +k1r54*x(60);%TSC2*
        
        dxdt(31) = k1r51*x(59) -k1r52*x(31);%TSC2(p)
        
        dxdt(32) = -k1r53*x(32)*x(30) +k2r53*x(60) +k1r55*x(33)...
                   -k1r56*x(34)*x(32)+k2r56*x(35);%Rheb
        
        dxdt(33) = k1r54*x(60)-k1r55*x(33);%Rheb(u)
        
        dxdt(34) = -k1r56*x(34)*x(32) +k2r56*x(35);%mTOR
        
        dxdt(35) = k1r56*x(34)*x(32) -k2r56*x(35) -k1r59*x(36)*x(35)...
                  +k2r59*x(62) +k1r60*x(62) -k1r63*x(38)*x(35)...
                  +k2r63*x(64)+k1r64*x(64);%mTOR*
        
        dxdt(36) = -k1r57*x(36)*x(19) +k2r57*x(61) -k1r59*x(36)*x(35)...
                   +k2r59*x(62) +k1r62*x(63);%S6K
        
        dxdt(37) = -k1r8*x(5)*x(37) +k2r8*x(42) +k1r9*x(42) +k1r58*x(61)...
                   +k1r60*x(62)-k1r61*x(37)*x(38)+k2r61*x(63);%S6K*
        
        dxdt(38) = -k1r61*x(37)*x(38) +k2r61*x(63) +k1r62*x(63)...
                   -k1r63*x(38)*x(35) +k2r63*x(64) -k1r65*x(38)*x(22)...
                   +k2r65*x(65) +k1r67*x(39);%S6P
        
        dxdt(39) = k1r64*x(64) +k1r66*x(65) -k1r67*x(39);%S6P(u)
        
        dxdt(40) = k1r3*x(5)*x(3) -k2r3*x(40) -k1r4*x(40);%CIRS
        
        dxdt(41) = k1r6*x(5)*x(19) -k2r6*x(41) -k1r7*x(41);%CIRSiM
        
        dxdt(42) = k1r8*x(5)*x(37) -k2r8*x(42) -k1r9*x(42);%CIRSiS
        
        dxdt(43) = k1r11*x(8)*x(3)-k2r11*x(43)-k1r12*x(43);%CRasIG
        
        dxdt(44) = k1r13*x(8)*x(6) -k2r13*x(44)-k1r14*x(44);%CRasIR
        
        dxdt(45) = k1r16*x(10)*x(13)-k2r16*x(45)-k1r17*x(45);%CRaf
       
        dxdt(46) = k1r18*x(11)*x(14)-k2r18*x(46)-k1r19*x(46);%CRafd
        
        dxdt(47) = k1r20*x(10)*x(22)-k2r20*x(47)-k1r21*x(47);%CRafiA
        
        dxdt(48) = k1r24*x(16)*x(15)-k2r24*x(48)-k1r25*x(48);%CMEK
        
        dxdt(49) = k1r26*x(17)*x(14)-k2r26*x(49)-k1r27*x(49);%CMEKd
        
        dxdt(50) = k1r28*x(18)*x(17)-k2r28*x(50)-k1r29*x(50);%CMAPK
        
        dxdt(51) = k1r30*x(19)*x(20)-k2r30*x(51)-k1r31*x(51);%CMAPKd
        
        dxdt(52) = k1r32*x(21)*x(6)-k2r32*x(52)-k1r33*x(52);%CAKT
        
        dxdt(53) = k1r35*x(21)*x(24)-k2r35*x(53)-k1r36*x(53);%CAKTiP
        
        dxdt(54) = k1r38*x(25)*x(22)-k2r38*x(54)-k1r39*x(54);%CAMPKiA
        
        dxdt(55) = k1r41*x(27)*x(22)-k2r41*x(55)-k1r42*x(55);%CGSK3iA
        
        dxdt(56) = k1r44*x(29)*x(25)-k2r44*x(56)-k1r45*x(56);%CTSC2aA
        
        dxdt(57) = k1r46*x(29)*x(27)-k2r46*x(57)-k1r47*x(57);%CTSC2aG
        
        dxdt(58) = k1r48*x(30)*x(19)-k2r48*x(58)-k1r49*x(58);%CTSC2iM
        
        dxdt(59) = k1r50*x(30)*x(22)-k2r50*x(59)-k1r51*x(59);%CTSC2iA
        
        dxdt(60) = k1r53*x(32)*x(30)-k2r53*x(60)-k1r54*x(60);%CRhebiT
        
        dxdt(61) = k1r57*x(36)*x(19)-k2r57*x(61)-k1r58*x(61);%CS6KM
        
        dxdt(62) = k1r59*x(36)*x(35)-k2r59*x(62)-k1r60*x(62);%CS6KT
        
        dxdt(63) = k1r61*x(37)*x(38)-k2r61*x(63)-k1r62*x(63);%CS6Kd
        
        dxdt(64) = k1r63*x(38)*x(35)-k2r63*x(64)-k1r64*x(64);%Cs6PiT
        
        dxdt(65) = k1r65*x(38)*x(22)-k2r65*x(65)-k1r66*x(65);%CS6PiA
        
 
        
    end
