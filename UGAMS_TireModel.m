function [FY,Muy,dfz,Fz0_,dpi,By,Cy,Ey,Dy,Kya_,SVy,SHy] = UGAMS_TireModel(SA,FZ,P,IA)
% UGA Motorsports Tire Model Magic Formula 6.1.2
% Coefficients imported from OptimumTire for Hoosier R20 43075 16 x 7.5
% Equations from Tire and Vehicle Dynamics Third Edition, Hans Pacejka

% Lateral Coefficients
    PCY1 = 1.47938;
    PDY1 = 2.30264;
    PDY2 = -0.181517;
    PDY3 = 12.2537;
    PEY1 = 0.258174;
    PEY2 = 0.134876;
    PEY3 = 0.964615;
    PEY4 = 0.0467532;
    PEY5 = -156.583;
    PKY1 = 33.3683;
    PKY2 = 3.92018;
    PKY3 = 0.963332;
    PKY4 = 5.07208;
    PKY5 = 49.6379;
    PKY6 = -3.67812;
    PKY7 = -1.32934;
    PHY1 = -0.000546807;
    PHY2 = 0.00152636;
    PVY1 = -0.00246833;
    PVY2 = -0.0198409;
    PVY3 = -0.0718706;
    PVY4 = 1.50852;
    PPY1 = 0.290986;
    PPY2 = 0.971999;
    PPY3 = -0.185848;
    PPY4 = 0.19433;
    PPY5 = -0.351119;
% Nominal Values
    FNOMIN = -1110;
    NOMPRES = 14 * 6894.76;

% Scaling Coefficients
    LFZO = 1;
    LCY = 1;
    LMUY = .66;
    LKY = 1;
    LKYC = 1;
    LVY = 1;
    LHY = 1;
    LEY = 1;

% Spin (set = 1 if turn slip neglected/ camber remains small)
    ZETA0 = 1;
    ZETA2 = 1;
    ZETA3 = 1;
    ZETA4 = 1;
   
% (4.E4)
    gammaAST = sin(IA);

% (4.E1)
    Fz0_ = LFZO.*FNOMIN;
   
% (4.E2a)
    dfz = (FZ-Fz0_)./Fz0_;

% (4.E2b)
    pi0 = NOMPRES;
    dpi = (P*6894.76-pi0)./pi0;

% (4.E21)
    Cy = PCY1.*LCY;

% (4.E23)
    Muy = (PDY1+PDY2.*dfz).*(1+PPY3.*dpi+PPY4.*dpi.^2).*(1-PDY3.*gammaAST.^2).*LMUY;

% (4.E22)
    Dy = Muy.*FZ.*ZETA2;

% (4.E25)
    Kya = PKY1.*Fz0_.*(1+PPY1.*dpi).*(1-PKY3.*abs(gammaAST)).*sin(PKY4.*atan((FZ./Fz0_)./((PKY2+PKY5.*gammaAST.^2).*(1+PPY2.*dpi)))).*ZETA3.*LKY;

% (4.E24)
    Ey = (PEY1+PEY2.*dfz).*(1+PEY5.*gammaAST.^2-(PEY3+PEY4.*gammaAST).*sign(SA)).*LEY;

% (4.E39)
    signKya = sign(Kya);
    signKya(~signKya) = 1;
    Kya_ = Kya + eps(Kya).*signKya;

% (4.E28)
    SVyg = FZ.*(PVY3+PVY4.*dfz).*gammaAST*ZETA2.*LKYC.*LMUY;

% (4.E30)
    Kyy0 = FZ.*(PKY6+PKY7.*dfz).*(1+PPY5.*dpi).*LKYC;
    
% (4.E29)
    SVy = FZ.*(PVY1+PVY2.*dfz).*LVY.*LMUY.*ZETA2+SVyg;

% (4.E27)
    SHy = (PHY1+PHY2.*dfz).*LHY+(Kyy0.*gammaAST-SVyg)./Kya_.*ZETA0+ZETA4-1;
% (4.E20) 
    alphay = SA+SHy;

% (4.E26)
    signCy = sign(Cy);
    signCy(~signCy) = 1;
    By = Kya./(Cy.*Dy+eps(Cy).*signCy);

% (4.E19)
    FY  = Dy.*sin(Cy.*atan(By.*alphay-Ey.*(By.*alphay-atan(By.*alphay))))+SVy;




