% SLAC-BIR-4 magnitude and phase
% Edward Green 2019

% parameters
noPoints    = 8001; % choose number of points


% simulation parameters
w1max = 200*2*pi; % rad/s
xiRef = w1max*0.95;

DOMax   = 5000; % Hz

zeta = 20;
kappa = atan(20);

Tp = 0.01; % pulse length (s)
moduleLength = Tp/4; % s

fAng = 45; % degrees
phi1 = pi + fAng/2*pi/180;
phi2 = -pi + fAng/2*pi/180;

tp1 = moduleLength;
tp2 = 2*moduleLength;
tp3 = 3*moduleLength;
tp4 = 4*moduleLength;

tEnv = linspace(0,Tp,noPoints);

%% Module 1

% calculate amplitude and frequency offsets
w1Amplitude1               = @(t,xi) xi*w1max*tanh(zeta*(tp1 - t)./moduleLength);
w1Frequency1             = @(t) 2*pi*DOMax*tan(kappa*(t/moduleLength))./tan(kappa);
w1Phase1             = @(t) -DOMax*2*pi*moduleLength/kappa/tan(kappa)*(log(cos(kappa*(t/moduleLength))./cos(kappa)) - ...
    log(cos(kappa)));

alphaMod1 = @(t,xi) atan(w1Amplitude1(t,xi)./w1Frequency1(t));
epsilon = 1e-12;
syms t xi;
dalphadtMod1 = matlabFunction(diff(sym(alphaMod1),t),'Vars',[t xi]);
dalphadtMod1 = @(t,xi) max(dalphadtMod1(t,xi),dalphadtMod1(t+epsilon,xi));

xAdiabaticPulseEnv1  = @(t,w1) 0.0 + w1/w1max.*w1Amplitude1(t,1).*cos(w1Phase1(t)).*((t>=0) & (t<tp1));
yAdiabaticPulseEnv1  = @(t,w1) 0.0 + w1/w1max.*w1Amplitude1(t,1).*sin(w1Phase1(t)).*((t>=0) & (t<tp1));

w1SLXEnv1     = @(t,w1,slac) 0.0 + slac*w1/xiRef.*dalphadtMod1(t,w1max\xiRef).*cos(w1Phase1(t)-pi/2).*((t>=0) & (t<tp1));
w1SLYEnv1     = @(t,w1,slac) 0.0 + slac*w1/xiRef.*dalphadtMod1(t,w1max\xiRef).*sin(w1Phase1(t)-pi/2).*((t>=0) & (t<tp1));

%% Module 2
% calculate amplitude and frequency offsets
w1Amplitude2               = @(t,xi) xi*w1max*tanh(zeta*(t-tp1)./moduleLength);
w1Frequency2             = @(t) -2*pi*DOMax*tan(kappa*((tp2-t)/moduleLength))./tan(kappa);
w1Phase2             = @(t) DOMax*2*pi*moduleLength/kappa/tan(kappa)*(log(cos(kappa)) - ...
       log(cos(kappa*((tp2 - t)/moduleLength))./cos(kappa)))+phi1;

alphaMod2 = @(t,xi) atan(w1Amplitude2(t,xi)./w1Frequency2(t));
epsilon = 1e-12;
syms t xi;
dalphadtMod2 = matlabFunction(diff(sym(alphaMod2),t),'Vars',[t xi]);
dalphadtMod2 = @(t,xi) max(dalphadtMod2(t,xi),dalphadtMod2(t+epsilon,xi));

xAdiabaticPulseEnv2  = @(t,w1) 0.0 + w1/w1max.*w1Amplitude2(t,1).*cos(w1Phase2(t)).*((t>=tp1) & (t<tp2));
yAdiabaticPulseEnv2  = @(t,w1) 0.0 + w1/w1max.*w1Amplitude2(t,1).*sin(w1Phase2(t)).*((t>=tp1) & (t<tp2));

w1SLXEnv2     = @(t,w1,slac) 0.0 + slac*w1/xiRef.*dalphadtMod2(t,w1max\xiRef).*cos(w1Phase2(t)-pi/2).*((t>=tp1) & (t<tp2));
w1SLYEnv2     = @(t,w1,slac) 0.0 + slac*w1/xiRef.*dalphadtMod2(t,w1max\xiRef).*sin(w1Phase2(t)-pi/2).*((t>=tp1) & (t<tp2));

%% Module 3
% calculate amplitude and frequency offsets
w1Amplitude3               = @(t,xi) xi*w1max*tanh(zeta*(tp3 - t)./moduleLength);
w1Frequency3             = @(t) 2*pi*DOMax*tan(kappa*((t-tp2)/moduleLength))./tan(kappa);
w1Phase3             = @(t) -DOMax*2*pi*moduleLength/kappa/tan(kappa)*(log(cos(kappa*((t-tp2)/moduleLength))./cos(kappa)) - ...
       log(cos(kappa)))+phi1;

alphaMod3 = @(t,xi) atan(w1Amplitude3(t,xi)./w1Frequency3(t));
epsilon = 1e-12;
syms t xi;
dalphadtMod3 = matlabFunction(diff(sym(alphaMod3),t),'Vars',[t xi]);
dalphadtMod3 = @(t,xi) max(dalphadtMod3(t,xi),dalphadtMod3(t+epsilon,xi));

xAdiabaticPulseEnv3  = @(t,w1) 0.0 + w1/w1max.*w1Amplitude3(t,1).*cos(w1Phase3(t)).*((t>=tp2) & (t<tp3));
yAdiabaticPulseEnv3  = @(t,w1) 0.0 + w1/w1max.*w1Amplitude3(t,1).*sin(w1Phase3(t)).*((t>=tp2) & (t<tp3));

w1SLXEnv3     = @(t,w1,slac) 0.0 + slac*w1/xiRef.*dalphadtMod3(t,w1max\xiRef).*cos(w1Phase3(t)-pi/2).*((t>=tp2) & (t<tp3));
w1SLYEnv3     = @(t,w1,slac) 0.0 + slac*w1/xiRef.*dalphadtMod3(t,w1max\xiRef).*sin(w1Phase3(t)-pi/2).*((t>=tp2) & (t<tp3));
%% Module 4

% calculate amplitude and frequency offsets
w1Amplitude4               = @(t,xi) xi*w1max*tanh(zeta*(t-tp3)./moduleLength);
w1Frequency4             = @(t) -2*pi*DOMax*tan(kappa*((tp4-t)/moduleLength))./tan(kappa);
w1Phase4             = @(t) DOMax*2*pi*moduleLength/kappa/tan(kappa)*(log(cos(kappa)) - ...
     log(cos(kappa*((tp4-t)/moduleLength))./cos(kappa)))+phi1-phi2;

alphaMod4 = @(t,xi) atan(w1Amplitude4(t,xi)./w1Frequency4(t));
epsilon = 1e-12;
syms t xi;
dalphadtMod4 = matlabFunction(diff(sym(alphaMod4),t),'Vars',[t xi]);
dalphadtMod4 = @(t,xi) max(dalphadtMod4(t,xi),dalphadtMod4(t+epsilon,xi));

xAdiabaticPulseEnv4  = @(t,w1) 0.0 + w1/w1max.*w1Amplitude4(t,1).*cos(w1Phase4(t)).*((t>=tp3) & (t<=tp4));
yAdiabaticPulseEnv4  = @(t,w1) 0.0 + w1/w1max.*w1Amplitude4(t,1).*sin(w1Phase4(t)).*((t>=tp3) & (t<=tp4));

w1SLXEnv4     = @(t,w1,slac) 0.0 + slac*w1/xiRef.*dalphadtMod4(t,w1max\xiRef).*cos(w1Phase4(t)-pi/2).*((t>=tp3) & (t<=tp4));
w1SLYEnv4     = @(t,w1,slac) 0.0 + slac*w1/xiRef.*dalphadtMod4(t,w1max\xiRef).*sin(w1Phase4(t)-pi/2).*((t>=tp3) & (t<=tp4));

xEnv = @(t,w1,slac) xAdiabaticPulseEnv1(t,w1) + xAdiabaticPulseEnv2(t,w1) + xAdiabaticPulseEnv3(t,w1) + xAdiabaticPulseEnv4(t,w1)...
                               + w1SLXEnv1(t,w1,slac) + w1SLXEnv2(t,w1,slac) + w1SLXEnv3(t,w1,slac) + w1SLXEnv4(t,w1,slac);
yEnv = @(t,w1,slac) yAdiabaticPulseEnv1(t,w1) + yAdiabaticPulseEnv2(t,w1) + yAdiabaticPulseEnv3(t,w1) + yAdiabaticPulseEnv4(t,w1)...
                              + w1SLYEnv1(t,w1,slac) + w1SLYEnv2(t,w1,slac) + w1SLYEnv3(t,w1,slac) + w1SLYEnv4(t,w1,slac);

%% 
% generate magnitude and phase points for BIR4, unscaled slac and sar-matched 
% SLAC-BIR-4
                          

% BIR-4 magnitude and phase
BIR4Mag             = sqrt(xEnv(tEnv,w1max,0).^2 + yEnv(tEnv,w1max,0).^2);
BIR4Phase           = atan2(yEnv(tEnv,w1max,0),xEnv(tEnv,w1max,0));

% SLAC-BIR-4  magnitude and phase (unscaled)
UnscaledSLACBIR4Mag = sqrt(xEnv(tEnv,w1max,1).^2 + yEnv(tEnv,w1max,1).^2);

% Calculate scaling factor for SAR-matching
SARrescalefactor = sqrt(sum(BIR4Mag.^2)/sum(UnscaledSLACBIR4Mag.^2));

% SLAC magnitude and phase (SAR-matched)
SARmatchedSLACBIR4Mag = sqrt(xEnv(tEnv,w1max*SARrescalefactor,1).^2 +...
                              yEnv(tEnv,w1max*SARrescalefactor,1).^2);
SLACBIR4Phase         = atan2(yEnv(tEnv,w1max,1),xEnv(tEnv,w1max,1));


% Normalise amplitude

normalisedBIR4Mag     = BIR4Mag/max(BIR4Mag(:));
normalisedSLACBIR4Mag = SARmatchedSLACBIR4Mag/max(SARmatchedSLACBIR4Mag(:));

% Scaling factor for max B1 amplitude (required for voltage setting on
% scanner)
SLACVoltageScale           = max(UnscaledSLACBIR4Mag(:))/max(BIR4Mag(:));
SARmatchedSLACVoltageScale = max(SLACBIR4MagEnv(:))/max(BIR4Mag(:));


