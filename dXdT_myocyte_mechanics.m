%%%% Differential equations for myocyte model

function [dYdT, dSL, F_XB, F_passive] = dXdT_myocyte_mechanics(t, y, para, Ca_fraction, flag, Ca_flag)
%% Model parameters
% Metabolite concentrations based on mean sham rat (Lopez et al. 2020)
MgATP = para(1); % mM
MgADP = para(2); % mM
Pi = para(3); % mM

% XB parameters
kstiff1 = para(4); % Stiffness constant due to myosin-actin interaction (kPa/um)
kstiff2 = para(5); % Stiffness constant due to working stroke of XBs (kPa/um)
k_passive = para(6); % Passive stiffness constant (kPa/um)
Kse = para(7); % Series element stiffness (mmHg/um) 
K_coop = para(8); % Strength of thin filament cooperativity
k_on = para(9); % Rate constant of Ca2+ binding to troponin C (uM^-1s^-1)
k_off = para(10); % Rate constant of Ca2+ unbinding from troponin C (s^-1)
kSR = para(11); % Forward rate constant of force dependent super-relaxed state transition (USR to N or P) (s^-1)
krecruit = para(12); % Force dependence of transition to super-relaxed state (N^-1 m^-1)
k_SR = para(13); % Reverse rate constant of force dependent super-relaxed transition (s^-1)
ka      = para(14); % Myosin actin associaiton rate (P to A1) (s^-1)
kd      = para(15); % Myosin actin dissociation rate (A1 to P) (s^-1) 
k1      = para(16); % A1 to A2 transition forward rate constant (s^-1)
k_1     = para(17); % A2 to A1 transition reverse rate constant (s^-1) 
k2      = para(18); % A2 to A3 transition forward rate constant (s^-1)
k_2     = para(19); % A3 to A2 transition reverse rate constant (s^-1) 
k3      = para(20); % A3 to P transition forward rate constant (s^-1)
visc = para(21); % Viscosity (mmHg*s/um)
stim_period = para(22); % Hz

Lsref = 1.8; % Sarcomere length at which passive force is zero (um)
alpha1  = 10; % Stretch sensing parameter for k1 and k_1 (1/um)
alpha2  = 9.1; % Stretch sensing parameter for k2 and k_2
alpha3  = 5.93; % Stretch sensing parameter for k3
s3      = 9.9e-3; % Strain at which k3 is minimum (um)
K_Pi = 4.00; % Pi dissociation constant (mM)
K_T = 0.4897 ; % MgATP dissociation constant (mM)
K_D = 0.194;  % MgADP dissociation constant (mM)
Lthin = 1.200; % Length of thin filament (nm)
Lthick = 1.670; % Length of thick filament (nm)
Lbare = 0.100; % Bare length of thin filament (nm)
dr = 0.01; % Power-stroke size (um)
SLcollagen = 2.25; % Threshold for collagen activation (um)
PConcollagen = 0.01; % Scale factor for passive force contributed by collagen
PExpcollagen = 70; % Expresion of collagen (1/um)

% Defining the metabolite dependent coeficient, rapid equilibrium of the
% cross bridge sub-states
g1 = (MgADP/K_D)/(1 + MgADP/K_D + MgATP/K_T);
g2 = (MgATP/K_T)/(1 + MgATP/K_T + MgADP/K_D);
f1 = (Pi/K_Pi)/(1 + Pi/K_Pi);
f2 = 1/(1 + Pi/K_Pi);

% Adjust for metabolite concentrations
kd = kd * f1; 
k1 = k1 * f2; 
k_2 = k_2 * g1; 
k3 = k3 * g2;

if flag == 0 % force pCa
SLset = 2.25; 
else if flag == 1 % Ktr
    if t <= 4.488 % Fixed length, isometric contraction
    SLset = 2.25;
    end
    if t > 4.488 && t <= 4.5469 % t > 4.488 && t <= 4.494 %length reduced by 54%, allowed to shorten
    SLset = 0.54*2.25; 
    end
    if t > 4.5469 % Return to L0
    SLset = 2.25;
    end
else if flag == 2 % Twitch
SLset = 1.84;
end
end
end

if flag == 0 || flag == 1 % force pCa or Ktr
    Ca_i = Ca_fraction;
       
else if flag == 2 % Twitch
        if Ca_flag == 0 % dATP = 0
            a = 0.106;
            b = 0.5635;
            c = 1.8017;
            Ca0 = 0;
        else
            a = 0.0534;
            b = 0.5484;
            c = 2.4732;
            Ca0 = 0;
        end
        phi = mod(t+0.001,stim_period)/stim_period;
        Ca_i = ((a/phi)*exp(-b*(log(phi)+c)^2) + Ca0);

end
end

%% State variables
% Moments of state variables
p1_0 = y(1);
p1_1 = y(2);
p1_2 = y(3);
p2_0 = y(4);
p2_1 = y(5);
p2_2 = y(6);
p3_0 = y(7);
p3_1 = y(8);
p3_2 = y(9);
N = y(10); % Non-permissible state 
P = 1.0 - N - p1_0 - p2_0 - p3_0; % Permssible state 
SL = y(11); % Sarcomere length (um)
U_NR = y(12); % Non-relaxed state 
U_SR = 1.0 - U_NR; % Super-relaxed state 

%% Stretch-sensitive rates
f_alpha1o = (p1_0 - alpha1*p1_1 + 0.5*(alpha1*alpha1)*p1_2);
f_alpha1i = (p1_1 - alpha1*p1_2);

f_alpha0o = (p2_0 + alpha1*p2_1 + 0.5*alpha1*alpha1*p2_2);
f_alpha0i = (p2_1 + alpha1*p2_2);

f_alpha2o = (p2_0 - alpha2*p2_1 + 0.5*(alpha2*alpha2)*p2_2);
f_alpha2i = (p2_1 - alpha2*p2_2);

f_alpha3o = (p3_0 + alpha3*(s3*s3*p3_0 + 2.0*s3*p3_1 + p3_2)); 
f_alpha3i = (p3_1 + alpha3*(s3*s3*p3_1 + 2.0*s3*p3_2));

%% Compute active and passive force
% Overlap function 
OV_Zaxis = min(Lthick/2,SL/2); % Overlap region closest to Z-axis (nm)
OV_Mline = max(SL/2-(SL-Lthin),Lbare/2); % Overal region closest to M-line (nm)
LOV = OV_Zaxis - OV_Mline; % Length of overlap (nm)
N_overlap_thick = LOV*2/(Lthick - Lbare); % Fraction of thick filament overlap

% Active Force
B_process = kstiff2 * dr * p3_0;   % Force due to XB ratcheting
C_process = kstiff1 * (p2_1 + p3_1 );% Force due to stretching of XBs
F_XB = N_overlap_thick*(B_process + C_process); % Active force

% (Linear) passive force model
% Collagen force
sigma_collagen  = PConcollagen*(exp(PExpcollagen*(SL - SLcollagen)) - 1)*(SL > SLcollagen);

F_passive  = k_passive*(SL/2-Lsref/2) + sigma_collagen; % Passive Force 

Ftotal = F_XB + F_passive; % Total force

%% Simulate various experimental protocols
if flag == 0 %force pCa
%dSL = 0; % Fully isometric
intf = (- Ftotal  + Kse*(SLset - SL));
dSL = intf/visc; % Muscle isometric but not sarcomere isometric

else if flag == 1 % Ktr
if t <= 4.488
    dSL = 0; % Isometric
else
    intf = (- Ftotal  + Kse*(SLset - SL));
    dSL = intf/visc;
end

else if flag == 2 % Twitch
intf = (-Ftotal + Kse*(SLset-SL));
%intf = -Ftotal; % Fully unloaded shortening
dSL = intf/visc;
%dSL = 0; % Fully isometric
end
end
end

dp1_0 = ka*P*N_overlap_thick*U_NR   - kd*p1_0 - k1*f_alpha1o + k_1*f_alpha0o;
dp1_1 = 1*dSL*p1_0 - kd*p1_1 - k1*f_alpha1i + k_1*f_alpha0i;
dp1_2 = 2*dSL*p1_1 - kd*p1_2 - k1*p1_2 + k_1*p2_2;
dp2_0 =         - k_1*f_alpha0o - k2*f_alpha2o + k_2*p3_0 + k1*f_alpha1o;
dp2_1 = 1*dSL*p2_0 - k_1*f_alpha0i - k2*f_alpha2i + k_2*p3_1 + k1*f_alpha1i;
dp2_2 = 2*dSL*p2_1 - k_1*p2_2       - k2*p2_2 + k_2*p3_2 + k1*p1_2;
dp3_0 =         + k2*f_alpha2o - k_2*p3_0 - k3*f_alpha3o;
dp3_1 = 1*dSL*p3_0 + k2*f_alpha2i - k_2*p3_1 - k3*f_alpha3i;
dp3_2 = 2*dSL*p3_1 + k2*p2_2       - k_2*p3_2 - k3*p3_2;


%% Campbell Ca activation model (Campbell et al. 2018)
Jon = k_on*Ca_i*N*(1 + K_coop*(1 - N));
Joff = k_off*P*(1 + K_coop*N);
dNp = - Jon + Joff;   

%% Transitions between super relaxed state and non relaxed state
dU_NR = kSR * ( 1 + krecruit * F_XB ) * U_SR - k_SR * U_NR;
dYdT = [dp1_0; dp1_1; dp1_2; dp2_0; dp2_1; dp2_2; dp3_0; dp3_1; dp3_2; dNp; dSL; dU_NR];

end