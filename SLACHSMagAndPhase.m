% SLAC-HS magnitude and phase
% Edward Green 2019

% parameters
noPoints = 8001; % choose number of points

beta = 3.4516;
mu   = 975.5;
Tp = 0.01024; % pulse length (s)
hsParams = [beta mu];

tEnv = linspace(0,Tp,noPoints);

% RF excitation
w1max    = 500*2*pi; % rad/s
xiRef    = w1max*0.7;

%% Create anonymous functions

% generating functions
w1Amplitude     = @(t,xi) xi.*w1max.*sech(beta*(2*t./Tp-1)); % rad/s
w1Frequency     = @(t) mu*beta*tanh(beta*(2*t./Tp-1));       % rad/s
w1Phase         = @(t) -Tp/2*mu*log(sech(beta*(2*t./Tp-1))); % rad/s

alpha = @(t,xi) atan(w1Amplitude(t,xi)./w1Frequency(t));
syms t xi;
dalphadt = matlabFunction(diff(sym(alpha),t),'Vars',[t xi]);

xAdiabaticPulseEnv  = @(t,w1) 0.0 + w1/w1max.*w1Amplitude(t,1).*cos(w1Phase(t));
yAdiabaticPulseEnv  = @(t,w1) 0.0 + w1/w1max.*w1Amplitude(t,1).*sin(w1Phase(t));

w1SLXEnv            = @(t,w1,slac) 0.0 + slac*w1/xiRef.*dalphadt(t,w1max\xiRef).*cos(w1Phase(t)-pi/2);
w1SLYEnv            = @(t,w1,slac) 0.0 + slac*w1/xiRef.*dalphadt(t,w1max\xiRef).*sin(w1Phase(t)-pi/2);

xEnv                = @(t,w1,slac) xAdiabaticPulseEnv(t,w1) + w1SLXEnv(t,w1,slac);
yEnv                = @(t,w1,slac) yAdiabaticPulseEnv(t,w1) + w1SLYEnv(t,w1,slac);


%% 
% generate magnitude and phase points for hs, unscaled slac and sar-matched 
% SLAC HS

% HS magnitude and phase
HSMag             = sqrt(xEnv(tEnv,w1max,0).^2 + yEnv(tEnv,w1max,0).^2);
HSPhase           = rem(w1Phase(tEnv)+pi,2*pi)-pi; % move phase into ±pi range

% SLAC HS  magnitude and phase (unscaled)
UnscaledSLACHSMag = sqrt(xEnv(tEnv,w1max,1).^2 + yEnv(tEnv,w1max,1).^2);

% Calculate scaling factor for SAR-matching
SARrescalefactor = sqrt(sum(HSMag.^2)/sum(UnscaledSLACHSMag.^2));

% SLAC magnitude and phase (SAR-matched)
SARmatchedSLACHSMag = sqrt(xEnv(tEnv,w1max*SARrescalefactor,1).^2 +...
                              yEnv(tEnv,w1max*SARrescalefactor,1).^2);
SLACHSPhase         = atan2(yEnv(tEnv,w1max,1),xEnv(tEnv,w1max,1));

% Normalise amplitude

normalisedHSMag     = HSMag/max(HSMag(:));
normalisedSLACHSMag = SARmatchedSLACHSMag/max(SARmatchedSLACHSMag(:));

% Scaling factor for max B1 amplitude (required for voltage setting on
% scanner)
SLACVoltageScale           = max(UnscaledSLACHSMag(:))/max(HSMag(:));
SARmatchedSLACVoltageScale = max(SLACHSMagEnv(:))/max(HSMag(:));



