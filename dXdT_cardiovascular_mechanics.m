%%%% Differential equations for ventricular model

function dXdT = dXdT_cardiovascular_mechanics(t, x, para)

%% Parameters
% XB parameters
ka = para(1); % Myosin actin associaiton rate (P to A1) (s^-1)
kd = para(2); % Myosin actin dissociation rate (A1 to P) (s^-1) 
k1 = para(3); % A1 to A2 transition forward rate constant (s^-1)
k_1 = para(4); % A2 to A1 transition reverse rate constant (s^-1) 
k2 = para(5); % A2 to A3 transition forward rate constant (s^-1)
k_2 = para(6); % A3 to A2 transition reverse rate constant (s^-1) 
k3 = para(7); % A3 to P transition forward rate constant (s^-1)
krecruit = para(8); % Force dependence of transition to super-relaxed state (N^-1 m^-1)
k_on = para(9); % Rate constant of Ca2+ binding to troponin C (uM^-1s^-1)
k_off = para(10); % Rate constant of Ca2+ unbinding from troponin C (s^-1)
k_coop = para(11); % Strength of thin filament cooperativity
visc = para(12); % Viscosity (mmHg s/um)
kstiff1 = para(13); % Stiffness constant due to myosin-actin interaction (kPa/um)
kstiff2 = para(14); % Stiffness constant due to working stroke of XBs (kPa/um)
kSR = para(15); % Forward rate constant of force dependent super-relaxed state transition (USR to N or P) (s^-1)
k_SR = para(16); % Reverse rate constant of force dependent super-relaxed transition (s^-1)
k_passive = para(17); % Passive stiffness constant (kPa/um)
Lthin = para(18); % Length of thin filament (nm)
Lthick = para(19); % Length of thick filament (nm)
Lbare  = para(20); % Bare length of thin filament (nm)
dr = para(21); % Power-stroke size (um)
SLcollagen = para(22); % Threshold for collagen activation (um)
PConcollagen = para(23); % Scale factor for passive force contributed by collagen
PExpcollagen = para(24); % um^-1
Lsref = para(25); % Resting sarcomere length (um)
L_rest_pas = para(26); % Length at which passive force = 0 (um)
K_T = para(27); % MgATP dissociation constant (mM)
K_D = para(28); % MgADP dissociation constant (mM)
K_Pi = para(29); % Pi dissociation constant (mM)
alpha1 = para(30); % Stretch sensing parameter for k1 and k_1 (1/um)
alpha2 = para(31); % Stretch sensing parameter for k2 and k_2 (1/um)
alpha3 = para(32); % Stretch sensing parameter for k3 (1/um)
s3 = para(33); % Strain at which k3 is minimum (um)
scale = para(34); % Scale factor for scaling myocyte force to ventricular level
kse = para(35); % Series element elastance (mmHg/um) 

% Heart model parameters
Vw_LV = para(36); % LV wall volume (mL)
Vw_SEP = para(37); % Septal wall volume (mL)
Vw_RV = para(38); % RV wall volume (mL)
Amref_LV = para(39); % LV midwall reference surface area (cm^2)
Amref_SEP = para(40); % Septal midwall reference surface area (cm^2)
Amref_RV = para(41); % RV midwall reference surface area (cm^2)

% Lumped circulatory parameters
C_Ao = para(42); % Proximal aortic compliance (mL/mmHg)
C_SA = para(43); % Systemic arterial compliance (mL/mmHg)
C_SV = para(44); % Systemic venous compliance (mL/mmHg) 
C_PV = para(45); % Pulmonary venous compliance (mL/mmHg)
C_PA = para(46); % Pulmonary arterial compliance (mL/mmHg)
R_Ao = para(47); % Resistance of aorta (mmHg*s/mL)
R_SA = para(48); % Systemic vasculature resistance (mmHg*s/mL)
R_PA = para(49); % Pulmonary vasculature resistance, mmHg*sec/mL
R_SV = para(50); % Resistance of systemic veins (mmHg*s/mL)
R_PV = para(51); % Resistance of pulmonary veins (mmHg*s/mL)
R_vlv = para(52); %  Valve resistance (mmHg*s/mL)
R_AV = para(53); % Resistance across aortic valve (mmHg/s*mL)
R_tAo = para(54); % Transmural aortic resistance (mmHg*s/mL)
R_tSA = para(55); % Transmural systemic arterial resistance (mmHg*s/mL)

% Metabolite concentrations
MgATP_LV  = para(56); % LV cytosolic MgATP concentration (mM)
MgADP_LV  = para(57); % LV cytosolic MgADP concentration (mM)
Pi_LV     = para(58); % LV cytosolic Pi concentration (mM)
MgATP_SEP = para(59); % Septal cytosolic MgATP concentration (mM)
MgADP_SEP = para(60); % Septal cytosolic MgADP concentration (mM)
Pi_SEP    = para(61); % Septal cytosolic Pi concentration (mM)
MgATP_RV  = para(62); % RV cytosolic MgATP concentration (mM)
MgADP_RV  = para(63); % RV cytosolic MgADP concentration (mM)
Pi_RV     = para(64); % RV cytosolic Pi concentration (mM)

flag_swap_metabolite = para(65); % 0 = healthy, 1 = HF
Ca_flag = para(66); % 0 = ATP, 1 = dATP
%XB_protocol = para(67); % 0 = Single simulation, 1 = Comparing ATP/dATP just Ca effects/dATP just XB effects/both simultaneously, 2 = Varying dATP percentages
r = para(67);
stim_period = para(68);
dATP = para(69);

% Correcting rate constants for metabolite levels in LV, SEP, and RV
kd_LV  = kd*(Pi_LV/K_Pi)/(1.0 + Pi_LV/K_Pi);
k1_LV  = k1/(1.0 + Pi_LV/K_Pi);
km2_LV = k_2*(MgADP_LV/K_D)/(1.0 + MgADP_LV/K_D + MgATP_LV/K_T);
k3_LV  = k3*(MgATP_LV/K_T)/(1.0 + MgATP_LV/K_T + MgADP_LV/K_D);

kd_SEP  = kd*(Pi_SEP/K_Pi)/(1.0 + Pi_SEP/K_Pi);
k1_SEP  = k1/(1.0 + Pi_SEP/K_Pi);
km2_SEP = k_2*(MgADP_SEP/K_D)/(1.0 + MgADP_SEP/K_D + MgATP_SEP/K_T);
k3_SEP  = k3*(MgATP_SEP/K_T)/(1.0 + MgATP_SEP/K_T + MgADP_SEP/K_D);

kd_RV  = kd*(Pi_RV/K_Pi)/(1.0 + Pi_RV/K_Pi);
k1_RV  = k1/(1.0 + Pi_RV/K_Pi);
km2_RV = k_2*(MgADP_RV/K_D)/(1.0 + MgADP_RV/K_D + MgATP_RV/K_T);
k3_RV  = k3*(MgATP_RV/K_T)/(1.0 + MgATP_RV/K_T + MgADP_RV/K_D);


%% Calcium interpolation function
if Ca_flag == 0 % Always use ATP Ca transient
    a = 0.3694;
    b = 1.6882;
    c = 1.6705;
    Ca0 = 0.1644;

else if Ca_flag == 1 % Use dATP Ca transient if dATP concentration > 0
if dATP == 0 % dATP = 0
    a = 0.3694;
    b = 1.6882;
    c = 1.6705;
    Ca0 = 0.1644;
else % dATP > 0
    a = 0.2026;
    b = 1.9203;
    c = 2.2884;
    Ca0 = 0.1644;
end
end
end

phi = mod(t+0.0001,stim_period)/stim_period;
Ca_i = (a/phi)*exp(-b*(log(phi)+c)^2) + Ca0;

%% State variables
xm_LV  = x(1); %?Maximal axial distance from LV midwall surface to origin (cm)
xm_SEP = x(2); %?Maximal axial distance from septal midwall surface to origin (cm)
xm_RV  = x(3); %?Maximal axial distance from RV midwall surface to origin (cm)
ym     = x(4); % Radius of midwall junction circle (cm)
SL_LV  = x(5); % LV sarcomere length (um)
SL_SEP = x(6); % Septal sarcomere length (um)
SL_RV  = x(7); % RV sarcomere length (um)
V_LV   = x(8); % LV volume (mL)
V_RV   = x(9); % RV volume (mL)

P1_0_LV = x(10); % 0th moment state A1, LV
P1_1_LV = x(11); % 1st moment state A1, LV
P1_2_LV = x(12); % 2nd moment state A1, LV
P2_0_LV = x(13); % 0th moment state A2, LV
P2_1_LV = x(14); % 1st moment state A2, LV
P2_2_LV = x(15); % 2nd moment state A2, LV
P3_0_LV = x(16); % 0th moment state A3, LV
P3_1_LV = x(17); % 1st moment state A3, LV
P3_2_LV = x(18); % 2nd moment state A3, LV
N_LV    = x(19); % Non-permissive state, LV
U_NR_LV = x(20); %  Non-super-relaxed state, LV

P1_0_SEP = x(21); % 0th moment state A1, SEP
P1_1_SEP = x(22); % 1st moment state A1, SEP
P1_2_SEP = x(23); % 2nd moment state A1, SEP
P2_0_SEP = x(24); % 0th moment state A2, SEP
P2_1_SEP = x(25); % 1st moment state A2, SEP
P2_2_SEP = x(26); % 2nd moment state A2, SEP
P3_0_SEP = x(27); % 0th moment state A3, SEP
P3_1_SEP = x(28); % 1st moment state A3, SEP
P3_2_SEP = x(29); % 2nd moment state A3, SEP
N_SEP    = x(30); % Non-permissive state, SEP
U_NR_SEP = x(31); % Non-super-relaxed state, SEP

P1_0_RV = x(32); % 0th moment state A1, RV
P1_1_RV = x(33); % 1st moment state A1, RV
P1_2_RV = x(34); % 2nd moment state A1, RV
P2_0_RV = x(35); % 0th moment state A2, RV
P2_1_RV = x(36); % 1st moment state A2, RV
P2_2_RV = x(37); % 2nd moment state A2, RV
P3_0_RV = x(38); % 0th moment state A3, RV
P3_1_RV = x(39); % 1st moment state A3, RV
P3_2_RV = x(40); % 2nd moment state A3, RV
N_RV    = x(41); % Non-permissive state, RV
U_NR_RV = x(42); % Non-super-relaxed state, RV

V_SV   = x(43); % Volume of systemic veins (mL)
V_PV   = x(44); % Volume of pulmonary veins (mL)
V_SA   = x(45); % Volume of systemic arterys (mL)
V_PA   = x(46); % Volume of pulmonary arterys (mL)
V_Ao   = x(47); % Volume of proximal aorta (mL)

%% Heart and Sarcomere Model
Vm_LV  = (pi/6)*xm_LV*(xm_LV^2 + 3*ym^2);
Vm_SEP = (pi/6)*xm_SEP*(xm_SEP^2 + 3*ym^2);
Vm_RV  = (pi/6)*xm_RV*(xm_RV^2 + 3*ym^2);
Am_LV  = pi*(xm_LV^2 + ym^2); % LV midwall surface area (cm^2)
Am_SEP = pi*(xm_SEP^2 + ym^2); % Septal midwall surface area (cm^2)
Am_RV  = pi*(xm_RV^2 + ym^2); % RV midwall surface area (cm^2)
Cm_LV  = 2*xm_LV/(xm_LV^2 + ym^2); % Curvature of midwall surface, LV
Cm_SEP = 2*xm_SEP/(xm_SEP^2 + ym^2); % Curvature of midwall surface, septum
Cm_RV  = 2*xm_RV/(xm_RV^2 + ym^2); % Curvature of midwall surface, RV
z_LV   = 3*Cm_LV*Vw_LV/(2*Am_LV);
z_SEP  = 3*Cm_SEP*Vw_SEP/(2*Am_SEP);
z_RV   = 3*Cm_RV*Vw_RV/(2*Am_RV);

epsf_LV = 0.5*log(Am_LV/Amref_LV) - 0.083333*z_LV^2 - 0.019*z_LV^4; % LV fiber strain
epsf_SEP = 0.5*log(Am_SEP/Amref_SEP) - 0.083333*z_SEP^2 - 0.019*z_SEP^4; % Septal fiber strain
epsf_RV = 0.5*log(Am_RV/Amref_RV) - 0.083333*z_RV^2 - 0.019*z_RV^4; % RV fiber strain

SLo_LV = Lsref*exp(epsf_LV); % LV sarcomere length (um)
SLo_SEP = Lsref*exp(epsf_SEP); % Sepal sarcomere length (um)
SLo_RV = Lsref*exp(epsf_RV); % RV sarcomere length (um)

% Collagen force
sigma_collagen_LV  = PConcollagen*(exp(PExpcollagen*(SLo_LV - SLcollagen)) - 1).*(SLo_LV > SLcollagen);
sigma_collagen_SEP = PConcollagen*(exp(PExpcollagen*(SLo_SEP - SLcollagen)) - 1).*(SLo_SEP > SLcollagen);
sigma_collagen_RV  = PConcollagen*(exp(PExpcollagen*(SLo_RV - SLcollagen)) - 1).*(SLo_RV > SLcollagen);

sigmapas_LV  = k_passive*(SLo_LV/2-L_rest_pas)  + sigma_collagen_LV ;
sigmapas_SEP = k_passive*(SLo_SEP/2-L_rest_pas) + sigma_collagen_SEP;
sigmapas_RV  = k_passive*(SLo_RV/2-L_rest_pas)   + sigma_collagen_RV;

% Sarcomere geometry
sovr_ze = min(Lthick*0.5, SL_LV*0.5); % Overlap region closest to Z-axis (nm), LV
sovr_cle = max(SL_LV*0.5 - (SL_LV-Lthin),Lbare*0.5); % Overal region closest to M-line (nm), LV
L_sovr = sovr_ze - sovr_cle; % Length of overlap (nm), LV
N_overlap_LV = L_sovr*2/(Lthick - Lbare); % Fraction of thick filament overlap, LV

sovr_ze = min(Lthick*0.5, SL_SEP*0.5); % Overlap region closest to Z-axis (nm), septum
sovr_cle = max(SL_SEP*0.5 - (SL_SEP-Lthin),Lbare*0.5); % Overal region closest to M-line (nm), septum
L_sovr = sovr_ze - sovr_cle; % Length of overlap (nm), septum
N_overlap_SEP = L_sovr*2/(Lthick - Lbare); % Fraction of thick filament overlap, septum

sovr_ze = min(Lthick*0.5, SL_RV*0.5); % Overlap region closest to Z-axis (nm), RV
sovr_cle = max(SL_RV*0.5 - (SL_RV-Lthin),Lbare*0.5); % Overal region closest to M-line (nm), RV
L_sovr = sovr_ze - sovr_cle; %Length of overlap (nm), RV
N_overlap_RV = L_sovr*2/(Lthick - Lbare); % Fraction of thick filament overlap, RV

% Active forces
sigmaact_LV  = N_overlap_LV*(kstiff2*dr*(P3_0_LV) + kstiff1*(P2_1_LV + P3_1_LV)); % mmHg * normalized force
sigmaact_SEP = N_overlap_SEP*(kstiff2*dr*(P3_0_SEP) + kstiff1*(P2_1_SEP+P3_1_SEP)); % mmHg * normalized force
sigmaact_RV  = N_overlap_RV*(kstiff2*dr*(P3_0_RV) + kstiff1*(P2_1_RV+P3_1_RV)); % mmHg * normalized force

% Total forces
sigmaf_LV = -kse*(SL_LV - SLo_LV); % LV force
sigmaf_LV = sigmaf_LV*scale; % LV force, scaled
sigmaf_SEP = -kse*(SL_SEP - SLo_SEP); % Septal foce
sigmaf_SEP = sigmaf_SEP*scale; % Septal force, scaled
sigmaf_RV = -kse*(SL_RV - SLo_RV); % RV force
sigmaf_RV = sigmaf_RV*scale; % RV force, scaled

Tm_LV  = (Vw_LV*sigmaf_LV/(2*Am_LV))*(1 + (z_LV^2)/3 + (z_LV^4)/5); % Tension in LV midwall
Tm_SEP = (Vw_SEP*sigmaf_SEP/(2*Am_SEP))*(1 + (z_SEP^2)/3 + (z_SEP^4)/5); % Tension in septal midwall
Tm_RV  = (Vw_RV*sigmaf_RV/(2*Am_RV))*(1 + (z_RV^2)/3 + (z_RV^4)/5); % Tension in RV midwall

Tx_LV  = Tm_LV*2*xm_LV*ym/(xm_LV^2 + ym^2); % Axial component of LV tension
Tx_SEP = Tm_SEP*2*xm_SEP*ym/(xm_SEP^2 + ym^2); % Axial component of septal tension
Tx_RV  = Tm_RV*2*xm_RV*ym/(xm_RV^2 + ym^2); % Axial component of RV tension

Ty_LV  = Tm_LV*(-xm_LV^2 + ym^2)/(xm_LV^2 + ym^2); % Radial component of LV tension
Ty_SEP = Tm_SEP*(-xm_SEP^2 + ym^2)/(xm_SEP^2 + ym^2); % Radial component of septal tension
Ty_RV  = Tm_RV*(-xm_RV^2 + ym^2)/(xm_RV^2 + ym^2); % Radial component of RV tension

% Ventricular pressure
ptrans_LV = 2*Tx_LV/ym; % LV transmural pressure (mmHg)
ptrans_RV = 2*Tx_RV/ym; % RV transmural pressure (mmHg)
P_LV = -ptrans_LV; % LV pressure (mmHg)
P_RV = ptrans_RV; % RV pressure (mmHg)

%% Lumped circulatory model
P_SV = V_SV/C_SV; % Systemic venous pressure (mmHg)
P_PV = V_PV/C_PV; % Pulmonary venous pressure (mmHg)
P_PA = V_PA/C_PA; % Pulmonary arterial pressure (mmHg)

% Ao valves closed equations
QOUT_LV = 0;
P_Ao = (C_SA*R_Ao*R_SA*V_Ao + C_SA*R_Ao*R_tSA*V_Ao + C_SA*R_SA*R_tSA*V_Ao + C_Ao*R_SA*R_tAo*V_SA + C_Ao*C_SA*P_SV*R_tSA*R_tAo)/(C_Ao*C_SA*(R_Ao*R_SA + R_Ao*R_tSA + R_SA*R_tSA + R_SA*R_tAo + R_tSA*R_tAo));
P_SA = (C_Ao*R_Ao*R_SA*V_SA + C_SA*R_SA*R_tSA*V_Ao + C_Ao*R_SA*R_tAo*V_SA + C_Ao*C_SA*P_SV*R_Ao*R_tSA + C_Ao*C_SA*P_SV*R_tSA*R_tAo)/(C_Ao*C_SA*(R_Ao*R_SA + R_Ao*R_tSA + R_SA*R_tSA + R_SA*R_tAo + R_tSA*R_tAo));
Q_Ao = -(C_Ao*R_SA*V_SA - C_SA*R_SA*V_Ao - C_SA*R_tSA*V_Ao + C_Ao*C_SA*P_SV*R_tSA)/(C_Ao*C_SA*(R_Ao*R_SA + R_Ao*R_tSA + R_SA*R_tSA + R_SA*R_tAo + R_tSA*R_tAo));
if (P_Ao < P_LV)*(V_LV>0) 
% Ao valve open equations 
P_SA    = (C_Ao*R_Ao*R_SA*R_tAo*V_SA + C_Ao*R_Ao*R_SA*R_AV*V_SA + C_SA*R_SA*R_tSA*R_AV*V_Ao + C_Ao*R_SA*R_tAo*R_AV*V_SA + C_Ao*C_SA*P_SV*R_Ao*R_tSA*R_tAo + C_Ao*C_SA*P_SV*R_Ao*R_tSA*R_AV + C_Ao*C_SA*P_LV*R_SA*R_tSA*R_tAo + C_Ao*C_SA*P_SV*R_tSA*R_tAo*R_AV)/(C_Ao*C_SA*(R_Ao*R_SA*R_tAo + R_Ao*R_SA*R_AV + R_Ao*R_tSA*R_tAo + R_Ao*R_tSA*R_AV + R_SA*R_tSA*R_tAo + R_SA*R_tSA*R_AV + R_SA*R_tAo*R_AV + R_tSA*R_tAo*R_AV));
QOUT_LV = -(C_SA*R_Ao*R_SA*V_Ao + C_SA*R_Ao*R_tSA*V_Ao + C_SA*R_SA*R_tSA*V_Ao + C_Ao*R_SA*R_tAo*V_SA - C_Ao*C_SA*P_LV*R_Ao*R_SA - C_Ao*C_SA*P_LV*R_Ao*R_tSA - C_Ao*C_SA*P_LV*R_SA*R_tSA - C_Ao*C_SA*P_LV*R_SA*R_tAo - C_Ao*C_SA*P_LV*R_tSA*R_tAo + C_Ao*C_SA*P_SV*R_tSA*R_tAo)/(C_Ao*C_SA*(R_Ao*R_SA*R_tAo + R_Ao*R_SA*R_AV + R_Ao*R_tSA*R_tAo + R_Ao*R_tSA*R_AV + R_SA*R_tSA*R_tAo + R_SA*R_tSA*R_AV + R_SA*R_tAo*R_AV + R_tSA*R_tAo*R_AV));
Q_Ao    = -(C_Ao*R_SA*R_tAo*V_SA + C_Ao*R_SA*R_AV*V_SA - C_SA*R_SA*R_AV*V_Ao - C_SA*R_tSA*R_AV*V_Ao - C_Ao*C_SA*P_LV*R_SA*R_tAo - C_Ao*C_SA*P_LV*R_tSA*R_tAo + C_Ao*C_SA*P_SV*R_tSA*R_tAo + C_Ao*C_SA*P_SV*R_tSA*R_AV)/(C_Ao*C_SA*(R_Ao*R_SA*R_tAo + R_Ao*R_SA*R_AV + R_Ao*R_tSA*R_tAo + R_Ao*R_tSA*R_AV + R_SA*R_tSA*R_tAo + R_SA*R_tSA*R_AV + R_SA*R_tAo*R_AV + R_tSA*R_tAo*R_AV));
end
QIN_LV  = max((P_PV - P_LV)/R_PV,0); 
QIN_RV  = max((P_SV - P_RV)/R_SV,0); 
QOUT_RV = max((P_RV - P_PA)/R_vlv,0)*(V_RV>0); 

% TriSeg equations
dXdT(1) = (-V_LV - 0.5*Vw_LV - 0.5*Vw_SEP + Vm_SEP - Vm_LV)/V_LV; % Xm_LV, maximal axial distance from LV midwall surface to origin (cm)
dXdT(2) = (Tx_LV + Tx_SEP + Tx_RV);  % Xm_RV, maximal axial distance from septal midwall surface to origin (cm)
dXdT(3) = (V_RV + 0.5*Vw_RV + 0.5*Vw_SEP + Vm_SEP - Vm_RV)/V_LV; % Xm_SEP, maximal axial distance from septal midwall surface to origin (cm)
dXdT(4) = (Ty_LV + Ty_SEP + Ty_RV);  % Ym, radius of midwall junction circle (cm)

%% Myofiber Mechanics: SL_LV
% Calculations for stretch-senstive rates    
f_alpha1o = (P1_0_LV - alpha1*P1_1_LV + 0.5*(alpha1*alpha1)*P1_2_LV);
f_alpha1i = (P1_1_LV - alpha1*P1_2_LV);

f_alpha0o = (P2_0_LV + alpha1*P2_1_LV + 0.5*alpha1*alpha1*P2_2_LV);
f_alpha0i = (P2_1_LV + alpha1*P2_2_LV);

f_alpha2o = (P2_0_LV - alpha2*P2_1_LV + 0.5*(alpha2*alpha2)*P2_2_LV);
f_alpha2i = (P2_1_LV - alpha2*P2_2_LV);

f_alpha3o = (P3_0_LV + alpha3*(s3*s3*P3_0_LV + 2.0*s3*P3_1_LV + P3_2_LV));
f_alpha3i = (P3_1_LV + alpha3*(s3*s3*P3_1_LV + 2.0*s3*P3_2_LV));

dSL_LV = (sigmaf_LV - sigmapas_LV - sigmaact_LV)/visc;

P0_LV = 1.0 - N_LV - P1_0_LV - P2_0_LV - P3_0_LV; 
dXdT(10) = ka*P0_LV*U_NR_LV*N_overlap_LV - kd_LV*P1_0_LV - k1_LV*f_alpha1o + k_1*f_alpha0o; 
dXdT(11) = dSL_LV*P1_0_LV - kd_LV*P1_1_LV - k1_LV*f_alpha1i + k_1*f_alpha0i;
dXdT(12) = 2*dSL_LV*P1_1_LV - kd_LV*P1_2_LV - k1_LV*P1_2_LV + k_1*P2_2_LV;

dXdT(13) = -k_1*f_alpha0o - k2*f_alpha2o + km2_LV*P3_0_LV + k1_LV*f_alpha1o;
dXdT(14) = dSL_LV*P2_0_LV - k_1*f_alpha0i - k2*f_alpha2i + km2_LV*P3_1_LV + k1_LV*f_alpha1i;
dXdT(15) = 2*dSL_LV*P2_1_LV - k_1*P2_2_LV - k2*P2_2_LV + km2_LV*P3_2_LV + k1_LV*P1_2_LV;

dXdT(16) = +k2*f_alpha2o - km2_LV*P3_0_LV - k3_LV*f_alpha3o;
dXdT(17) = dSL_LV*P3_0_LV + k2*f_alpha2i - km2_LV*P3_1_LV - k3_LV*f_alpha3i;
dXdT(18) = 2*dSL_LV*P3_1_LV + k2*P2_2_LV - km2_LV*P3_2_LV - k3_LV*P3_2_LV;


U_SR_LV = 1 - U_NR_LV;
Jon = k_on*Ca_i*N_LV*(1 + k_coop*(1 - N_LV));
Joff = k_off*P0_LV*(1 + k_coop*N_LV);
dXdT(19) = - Jon + Joff; 
dXdT(20) = kSR * (1 + krecruit * sigmaact_LV) * U_SR_LV - k_SR*U_NR_LV;

%% Myofiber Mechanics: SL_SEP
% Calculations for stretch-senstive rates    
f_alpha1o = (P1_0_SEP - alpha1*P1_1_SEP + 0.5*(alpha1*alpha1)*P1_2_SEP);
f_alpha1i = (P1_1_SEP - alpha1*P1_2_SEP);

f_alpha0o = (P2_0_SEP + alpha1*P2_1_SEP + 0.5*alpha1*alpha1*P2_2_SEP);
f_alpha0i = (P2_1_SEP + alpha1*P2_2_SEP);

f_alpha2o = (P2_0_SEP - alpha2*P2_1_SEP + 0.5*(alpha2*alpha2)*P2_2_SEP);
f_alpha2i = (P2_1_SEP - alpha2*P2_2_SEP);

f_alpha3o = (P3_0_SEP + alpha3*(s3*s3*P3_0_SEP + 2.0*s3*P3_1_SEP + P3_2_SEP));
f_alpha3i = (P3_1_SEP + alpha3*(s3*s3*P3_1_SEP + 2.0*s3*P3_2_SEP));

% XB ODEs
dSL_SEP = (sigmaf_SEP - sigmapas_SEP - sigmaact_SEP)/visc;
 
P0_SEP = 1 - N_SEP - P1_0_SEP - P2_0_SEP - P3_0_SEP; 
dXdT(21) = ka*P0_SEP*U_NR_SEP*N_overlap_SEP - kd_SEP*P1_0_SEP - k1_SEP*f_alpha1o + k_1*f_alpha0o; 
dXdT(22) = dSL_SEP*P1_0_SEP - kd_SEP*P1_1_SEP - k1_SEP*f_alpha1i + k_1*f_alpha0i;
dXdT(23) = 2*dSL_SEP*P1_1_SEP - kd_SEP*P1_2_SEP - k1_SEP*P1_2_SEP + k_1*P2_2_SEP;

dXdT(24) = -k_1*f_alpha0o - k2*f_alpha2o + km2_SEP*P3_0_SEP + k1_SEP*f_alpha1o;
dXdT(25) = dSL_SEP*P2_0_SEP - k_1*f_alpha0i - k2*f_alpha2i + km2_SEP*P3_1_SEP + k1_SEP*f_alpha1i;
dXdT(26) = 2*dSL_SEP*P2_1_SEP - k_1*P2_2_SEP - k2*P2_2_SEP + km2_SEP*P3_2_SEP + k1_SEP*P1_2_SEP;

dXdT(27) = +k2*f_alpha2o - km2_SEP*P3_0_SEP - k3_SEP*f_alpha3o;
dXdT(28) = dSL_SEP*P3_0_SEP + k2*f_alpha2i - km2_SEP*P3_1_SEP - k3_SEP*f_alpha3i;
dXdT(29) = 2*dSL_SEP*P3_1_SEP + k2*P2_2_SEP - km2_SEP*P3_2_SEP - k3_SEP*P3_2_SEP;

U_SR_SEP = 1 - U_NR_SEP;
Jon = k_on*Ca_i*N_SEP*(1 + k_coop*(1 - N_SEP));
Joff = k_off*P0_SEP*(1 + k_coop*N_SEP);
dXdT(30) = - Jon + Joff; 
dXdT(31) = kSR*(1 + krecruit * sigmaact_SEP) * U_SR_SEP - k_SR*U_NR_SEP;

%% Myofiber Mechanics: SL_RV
% Calculations for stretch-senstive rates    
f_alpha1o = (P1_0_RV - alpha1*P1_1_RV + 0.5*(alpha1*alpha1)*P1_2_RV);
f_alpha1i = (P1_1_RV - alpha1*P1_2_RV);

f_alpha0o = (P2_0_RV + alpha1*P2_1_RV + 0.5*alpha1*alpha1*P2_2_RV);
f_alpha0i = (P2_1_RV + alpha1*P2_2_RV);

f_alpha2o = (P2_0_RV - alpha2*P2_1_RV + 0.5*(alpha2*alpha2)*P2_2_RV);
f_alpha2i = (P2_1_RV - alpha2*P2_2_RV);

f_alpha3o = (P3_0_RV + alpha3*(s3*s3*P3_0_RV + 2.0*s3*P3_1_RV + P3_2_RV));
f_alpha3i = (P3_1_RV + alpha3*(s3*s3*P3_1_RV + 2.0*s3*P3_2_RV));

% XB ODEs
dSL_RV = (sigmaf_RV - sigmapas_RV - sigmaact_RV)/visc; 

P0_RV = 1.0 - N_RV - P1_0_RV - P2_0_RV - P3_0_RV;  
dXdT(32) = ka*P0_RV*U_NR_RV*N_overlap_RV - kd_RV*P1_0_RV - k1_RV*f_alpha1o + k_1*f_alpha0o;
dXdT(33) = dSL_RV*P1_0_RV - kd_RV*P1_1_RV - k1_RV*f_alpha1i + k_1*f_alpha0i;
dXdT(34) = 2*dSL_RV*P1_1_RV - kd_RV*P1_2_RV - k1_RV*P1_2_RV + k_1*P2_2_RV;

dXdT(35) = -k_1*f_alpha0o - k2*f_alpha2o + km2_RV*P3_0_RV + k1_RV*f_alpha1o;
dXdT(36) = dSL_RV*P2_0_RV - k_1*f_alpha0i - k2*f_alpha2i + km2_RV*P3_1_RV + k1_RV*f_alpha1i;
dXdT(37) = 2*dSL_RV*P2_1_RV - k_1*P2_2_RV - k2*P2_2_RV + km2_RV*P3_2_RV + k1_RV*P1_2_RV;

dXdT(38) = +k2*f_alpha2o - km2_RV*P3_0_RV - k3_RV*f_alpha3o;
dXdT(39) = dSL_RV*P3_0_RV + k2*f_alpha2i - km2_RV*P3_1_RV - k3_RV*f_alpha3i;
dXdT(40) = 2*dSL_RV*P3_1_RV + k2*P2_2_RV - km2_RV*P3_2_RV - k3_RV*P3_2_RV;

U_SR_RV = 1.0 - U_NR_RV;
Jon = k_on*Ca_i*N_RV*(1 + k_coop*(1 - N_RV));
Joff = k_off*P0_RV*(1 + k_coop*N_RV);
dXdT(41) = - Jon + Joff; 
dXdT(42) = kSR*(1 + krecruit * sigmaact_RV) * U_SR_RV - k_SR * U_NR_RV; 

% Myofiber Mechanics
dXdT(5) = dSL_LV;
dXdT(6) = dSL_SEP;
dXdT(7) = dSL_RV;

% Lumped circulation model variables
dXdT(8) = QIN_LV - QOUT_LV; % V_LV
dXdT(9) = QIN_RV - QOUT_RV; % V_RV
dXdT(43) = (P_SA - P_SV)/R_SA - QIN_RV;  % V_SV
dXdT(44) = (P_PA - P_PV)/R_PA - QIN_LV;  % V_PV
dXdT(45) = Q_Ao - (P_SA - P_SV)/R_SA; % V_SA 
dXdT(46) = QOUT_RV - (P_PA - P_PV)/R_PA; % V_PA 
dXdT(47) = QOUT_LV - Q_Ao; % V_Ao

dXdT = dXdT(:);
