%%%% Main code to run biventricular mechanoenergetics and circulatory model

% Adapted by Abigail Teitgen based on code from: 
% Lopez, R. Marzban, B., Gao, X. Lauinger, E., Van den Bergh, F.,
% Whitesall, S. E., Converso-Baran, K., Burant, C. F., Michele, D. E.,
% Beard, D. A. (2020). Impaired myocardial energetics causes mechanical dysfunction
% in decompensated failing hearts. Function, 1(2): zqaa018

% Marzban, B., Lopez, R., Beard, D. A. (2020). Computational Modeling of Coupled 
% Energetics and Mechanics in the Rat Ventricular Myocardium. Physiome. (supplement, has thorough description of model equations and structure)

function [V_LV_store_F, P_LV_store_F, max_force, FS, CO, EF, LVDP, plus_dPdt_F, minus_dPdt_F, work_rate_F, ATP_F, ADP_F, Pi_F, MVO2_F, PCrATP_F, XB_turnover_F, ATPase_F, efficiency_F, SRX_store_F] = cardiovascular_model(Ca_value, HF, plotting, dATP_percent, ka_i, kd_i, k1_i, k_1_i, k2_i, k_2_i, k3_i, krecruit_i, k_on_i, k_off_i, k_coop_i)

tic 
     
%% Flags
% Calcium
Ca_flag = Ca_value; % 0 = ATP, 1 = dATP

% Healthy or HF
HF_protocol = HF; % 0 = Healthy, 1 = Heart failure

% Suppress plotting and calculation outputs
suppress_plotting = plotting; % 0 = do not plot, 1 = plot

% Percent dATP
dATP = dATP_percent;
ATP = 100 - dATP;

for r = 1:length(dATP)
if HF_protocol == 0
    flag_swap_metabolite = 0; % Healthy
else
    flag_swap_metabolite = 1; % Failure
end

%% Set parameters
% XB parameters
ka = 250*(ATP/100) + ka_i*(dATP/100); % Myosin actin associaiton rate (P to A1) (s^-1)
kd = 304.7*(ATP/100) + kd_i*(dATP/100); % Myosin actin dissociation rate (A1 to P) (s^-1) 
k1 = 4*(ATP/100) + k1_i*(dATP/100); % A1 to A2 transition forward rate constant (s^-1)
k_1 = 2*(ATP/100) + k_1_i*(dATP/100); % A2 to A1 transition reverse rate constant (s^-1) 
k2 = 80*(ATP/100) + k2_i*(dATP/100); % A2 to A3 transition forward rate constant (s^-1)
k_2 = 4*(ATP/100) + k_2_i*(dATP/100); % A3 to A2 transition reverse rate constant (s^-1) 
k3 = 25*(ATP/100) + k3_i*(dATP/100); % A3 to P transition forward rate constant (s^-1)
krecruit = 0.4*(ATP/100) + krecruit_i*(dATP/100); % Force dependence of transition to super-relaxed state (N^-1 m^-1)
k_on = k_on_i; % Rate constant of Ca2+ binding to troponin C (uM^-1s^-1)
k_off = k_off_i; % Rate constant of Ca2+ unbinding from troponin C (s^-1)
k_coop = k_coop_i; % Strength of thin filament cooperativity

visc = 0.01; % Viscosity (mmHg s/um)
kstiff1 = 1.8219e+04; % Stiffness constant due to myosin-actin interaction (kPa/um)
kstiff2 = 4.7822e+05; % Stiffness constant due to working stroke of XBs (kPa/um)
kSR = 1.3; % Forward rate constant of force dependent super-relaxed state transition (USR to N or P) (s^-1)
k_SR = 30; % Reverse rate constant of force dependent super-relaxed transition (s^-1)
k_passive = 2000; % Passive stiffness constant (kPa/um)
Lthin = 1.200; % Length of thin filament (nm)
Lthick = 1.670; % Length of thick filament (nm)
Lbare = 0.100; % Bare length of thin filament (nm)
dr = 0.01; % Power-stroke size (um)
SLcollagen = 2.25; % Threshold for collagen activation (um)
PConcollagen = 0.01; % Scale factor for passive force contributed by collagen
PExpcollagen = 70; % um^-1
Lsref = 1.9; % Resting sarcomere length (um)
L_rest_pas = 1; % Length at which passive force = 0 (um)
K_T = 0.4897; % MgATP dissociation constant (mM)
K_D = 0.194; % MgADP dissociation constant (mM)
K_Pi = 4.00; %Pi dissociation constant (mM)
alpha1  = 10.0; % Stretch sensing parameter for k1 and k_1 (1/um)
alpha2  = 9.1; % Stretch sensing parameter for k2 and k_2 (1/um)
alpha3 = 0.1*59.3; % Stretch sensing parameter for k3 (1/um)
s3 = 9.9e-3;  % Strain at which k3 is minimum (um)
scale = 0.12; % Scale factor for scaling myocyte force to ventricular level
kse = 100000; % Series element elastance (mmHg/um) 

%% Specifying which simulation should be run
rat_number = 9; % Mean sham rat 

%% Reading the adujstable scale variables for each individual rat
adjvar_all_rest = xlsread('Adjustable_parameters_table_rest.xlsx','B2:J21'); 
adjvar_all_swap = xlsread('Adjustable_parameters_table_swap.xlsx','B2:J21'); 

adjvar_all = adjvar_all_rest;
if flag_swap_metabolite == 1 % Failure
 adjvar_all = adjvar_all_swap;  
end

rate_of_XB_turnover_mean_sham = 5.0147;

if rat_number<=9
    shamRat = 1;
    delta_p = 0;
else
    shamRat = 0;
end

%%  Adjustable variables 
% adjvar = [Reference area LV , Reference area Septal,  Reference area RV, ksr& krecruit, Blood volume, R_SA, R_TAC, ATP_tune_Coeff]
adjvar = adjvar_all(rat_number,:);

%% Read the experimental data for SHAM and TAC rats from the excel file 
data = xlsread('data1.xlsx','A3:W23');
BW  = data(rat_number, 1); % g
LVW = data(rat_number, 2); % mg
RVW = data(rat_number, 3); % mg
LW = data(rat_number, 5); % mg

HR = 454; % beats/min, based on experimental protocol from Nowakowski et al. 2013 

edLV_target = data(rat_number,13);  % uL
esLV_target = data(rat_number,14); % uL

SV_LV_target = 30*1000/454; %mL, based on experimental protocol from Nowakowski et al. 2013 
CO_target = 30; % mL/min Based on experimental protocol from Nowakowski et al. 2013 
EF_LV_target = 80; % Percent, Based on experimental protocol from Nowakowski et al. 2013 
TAN = data(rat_number,16)/1000; % Total adenine nucleotide pool, mol/L cell
CRtot = data(rat_number,18)/1000; % Total creatine pool, mol/L cell
Ox_capacity = data(rat_number,21)/data(9,21); % Computed relative to mean sham
Ox_capacity_sham = 1; % Healthy
Ox_capacity_TAC = data(20,21)/data(9,21); % Failing

% Average sham
TAN_sham = data(9,16)/1000; % Total adenine nucleotide pool (healthy), mol/L cell
CRtot_sham = data(9,18)/1000; % Total creatine pool (healthy), mol/L cell

% Average TAC
TAN_TAC =  data(20,16)/1000; % Total adenine nucleotide pool (failure), mol/L cell
CRtot_TAC = data(20,18)/1000; % Total creatine pool (failure), mol/L cell

Ao = 10.260e-3; % mol/L cell
Co = 43.007e-3; % mol/L cell
Po = 35.446e-3; % mol/L cell

TEP = Po - (0.283e-3)*(Ao-TAN)/(0.082e-3); % Total exchangeable phosphate pool, mol/L cell
TEP_sham = Po - (0.283e-3)*(Ao-TAN_sham)/(0.082e-3); % Total exchangeable phosphate pool (healthy), mol/L cell
TEP_TAC = Po - (0.283e-3)*(Ao-TAN_TAC)/(0.082e-3); % Total exchangeable phosphate pool (failure), mol/L cell

if flag_swap_metabolite ==1
    if shamRat == 0
    TAN = TAN_sham;
    CRtot = CRtot_sham;
    TEP = TEP_sham;
    Ox_capacity = Ox_capacity_sham;
    else 
    TAN = TAN_TAC;
    CRtot = CRtot_TAC;
    TEP = TEP_TAC;
    Ox_capacity = Ox_capacity_TAC;
    end
end

if shamRat == 0
    preV = data(rat_number , 11); % mm/s
    postV = data(rat_number , 12); % mm/s
    postV = postV/1000; % m/s
    preV = preV/1000; % m/s
    rho_blood = 1060; % kg/m^3
    delta_p = 0.5*(postV^2-preV^2)*rho_blood; % Pa
    delta_p = 0.0075*delta_p; % mmHg
    if delta_p > 60 % Replacing delta_p with Max delta_P for the animals without measured post TAC velocity
        delta_p = 31.48;
    end
end

R_TAC = adjvar(8); % Resistance due to transverse aortic constriction (mmHg*sec/mL)

tune_ATPase_LV =  adjvar(9)* (1/ 0.6801) *1.0e-3; % ATP hydrolysis rate: M/s/L cytosol

Amref_LV  = adjvar(1)*1.0253; % LV midwall reference surface area (cm^2)
Amref_SEP = adjvar(2)*1.0253; % Septal midwall reference surface area (cm^2)
Amref_RV  = adjvar(3)*1.0253; % RV midwall reference surface area (cm^2)

Vw_LV = (LVW*2/3)/1000; % LV wall volume (mL)
Vw_SEP =(LVW/3)/1000; % Septal wall volume (mL)
Vw_RV = RVW/1000; % RV wall volume (mL)

% Slight differences in geometry for model stability
if flag_swap_metabolite == 1 % Failure
Vw_LV = Vw_LV*1.0174; 
Vw_SEP = Vw_SEP*1.0174; 
Vw_RV = Vw_RV*1.0174; 
else % Healthy
Vw_LV = Vw_LV*1.0143; %1.0142
Vw_SEP = Vw_SEP*1.0143; 
Vw_RV = Vw_RV*1.0143; 
end

% Lumped circulatory parameters
C_Ao = 0.0022045;  % Proximal aortic compliance (mL/mmHg)
C_SA = 0.0077157; % Systemic arterial compliance (mL/mmHg)
C_SV = 2.5; % Systemic venous compliance (mL/mmHg) 
C_PV = 0.25; % Pulmonary venous compliance (mL/mmHg)
C_PA = 0.013778; % Pulmonary arterial compliance (mL/mmHg)
R_Ao   = 2.5; % Resistance of aorta (mmHg*s/mL)
R_Ao = R_Ao/3;
R_SA   = adjvar(7); % Systemic vasculature resistance (mmHg*s/mL)
if flag_swap_metabolite == 1
    R_PA   = 12/CO_target*60; % Pulmonary vasculature resistance (mmHg*sec/mL)
else 
    R_PA   = adjvar(7)*12/88; % Pulmonary vasculature resistance (mmHg*sec/mL)
end
R_SV   = 0.25; % Resistance of systemic veins (mmHg*s/mL)
R_PV   = 0.25; % Resistance of pulmonary veins (mmHg*s/mL)
R_vlv  = 0.24; %  Valve resistance (mmHg*s/mL)
R_AV   = R_vlv + R_TAC; % Resistance across aortic valve (mmHg/s*mL)
R_tAo  = 0.5; % Transmural aortic resistance (mmHg*s/mL)
R_tSA  = 4; % Transmural systemic arterial resistance (mmHg*s/mL)

%% Run energetics model to get initial metabolite concentrations
energtics_output  = energetics_model(TAN, CRtot, TEP, Ox_capacity, tune_ATPase_LV);

MgATP_cyto = energtics_output(1); % Cytosolic MgATP concentration (M·(L cytosol water)^-1)
MgADP_cyto = energtics_output(2); % Cytosolic MgADP concentration (M·(L cytosol water)^-1)
fPi_cytoplasm = energtics_output(3); % Cytosolic unchelated Pi concentration (M·(L cytosol water)^-1)
MVO2_tissue = energtics_output(5); % Oxygen consumption rate (uM·min^-1·(g tissue)^-1)
dGrATPase = energtics_output(6); % ATP hydrolysis free energy, kJ/mol

PCrATP =  energtics_output(7); % Creatine phoosphate ATP ratio, unitless
ATP_cyto = energtics_output(8); % Cytosolic total ATP concentration (M·(L cytosol water)^-1)
ADP_cyto = energtics_output(9); % Cytosolic total ADP concentration (M·(L cytosol water)^-1)
Pi_cyto = energtics_output(10)*1000; % Cytosolic total Pi concentration (M·(L cytosol water)^-1)

% Store outputs
ATP_store(1) = MgATP_cyto;
ADP_store(1) = MgADP_cyto*1000;
Pi_store(1) = Pi_cyto;
MVO2_tissue_store(1) = MVO2_tissue;
PCrATP_store(1) = PCrATP;
   
%% Run cardiovascular mechanics model (run initially to steady state without coupled energetics model)
stim_period = 1/(HR/60); % Stimulation period (Hz)
para = [ka, kd, k1, k_1, k2, k_2, k3, krecruit, k_on ,k_off, k_coop, visc, kstiff1, kstiff2, kSR, k_SR, k_passive,... 
        Lthin, Lthick, Lbare, dr, SLcollagen, PConcollagen, PExpcollagen, Lsref, L_rest_pas, K_T, K_D, K_Pi, alpha1,...
        alpha2, alpha3, s3, scale, kse, Vw_LV, Vw_SEP, Vw_RV, Amref_LV, Amref_SEP, Amref_RV, C_Ao, C_SA, C_SV, C_PV,...
        C_PA, R_Ao, R_SA, R_PA, R_SV, R_PV, R_vlv, R_AV, R_tAo, R_tSA, MgATP_cyto, MgADP_cyto, Pi_cyto, MgATP_cyto,... 
        MgADP_cyto, Pi_cyto, MgATP_cyto, MgADP_cyto, Pi_cyto, flag_swap_metabolite, Ca_flag, r, stim_period, dATP];

M = speye(47); % Sparse identity matrix
M(1,1) = 0;
M(2,2) = 0;
M(3,3) = 0;
M(4,4) = 0; 

% Initial conditions for state variables
P1_0_LV = 0; % 0th moment state A1, LV
P1_1_LV = 0; % 1st moment state A1, LV
P1_2_LV = 0; % 2nd moment state A1, LV
P2_0_LV = 0; % 0th moment state A2, LV
P2_1_LV = 0; % 1st moment state A2, LV
P2_2_LV = 0; % 2nd moment state A2, LV
P3_0_LV = 0; % 0th moment state A3, LV
P3_1_LV = 0; % 1st moment state A3, LV
P3_2_LV = 0; % 2nd moment state A3, LV
N_LV = 1; % Non-permissive state, LV
U_NR_LV = 0; % Non-super-relaxed state, LV
P1_0_RV = 0; % 0th moment state A1, RV
P1_1_RV = 0; % 1st moment state A1, RV
P1_2_RV = 0; % 2nd moment state A1, RV
P2_0_RV = 0; % 0th moment state A2, RV
P2_1_RV = 0; % 1st moment state A2, RV
P2_2_RV = 0; % 2nd moment state A2, RV
P3_0_RV = 0; % 0th moment state A3, RV
P3_1_RV = 0; % 1st moment state A3, RV
P3_2_RV = 0; % 2nd moment state A3, RV
N_RV = 1; % Non-permissive state, RV
U_NR_RV = 0; % Non-super-relaxed state, RV
P1_0_SEP = 0; % 0th moment state A1, septum
P1_1_SEP = 0; % 1st moment state A1, septum
P1_2_SEP = 0; % 2nd moment state A1, septum
P2_0_SEP = 0; % 0th moment state A2, septum
P2_1_SEP= 0; % 1st moment state A2, septum
P2_2_SEP = 0; % 2nd moment state A2, septum
P3_0_SEP = 0; % 0th moment state A3, septum
P3_1_SEP = 0; % 1st moment state A3, septum
P3_2_SEP = 0; % 2nd moment state A3, septum
N_SEP = 1; % Non-permissive state, septum
U_NR_SEP = 0; % Non-super-relaxed state, septum

% Assign initial condtion for LV and RV
V_LV  = edLV_target/1000 + 0.1; % Intial LV volume (mL)
V_RV  = edLV_target/1000 + 0.1; % Initial RV volume (mL)

% Initial volumes
V_SA = adjvar(6)* 3.0; % Volume of systemic arteries (mL)
V_SV = adjvar(6)* 4.80; % Volume of systemic veins (mL)
V_PA = adjvar(6)* 0.5; % Volume of pulmonary arteries (mL)
V_PV = adjvar(6)* 1.0; % Volume of pulmonary veins (mL)
V_Ao = adjvar(6)* 1.0; % Volume of proximal aorta (mL)

% Initial conditons for geometry 
xm_LV   = -0.60; % Maximal axial distance from LV midwall surface to origin (cm)
xm_SEP  = 0.40; % Maximal axial distance from septal midwall surface to origin (cm)
xm_RV   = 1.0; % Maximal axial distance from RV midwall surface to origin (cm)
ym    = 0.50; % Radius of midwall junction circle (cm)

SL_LV  = 2.2; % LV sarcomere length (um)
SL_SEP = 2.2; % Septal sarcomere length (um)
SL_RV  = 2.2; % RV sarcomere length (um)

% Initial conditions 
init = [xm_LV ,xm_SEP ,xm_RV ,ym , SL_LV, SL_SEP, SL_RV, V_LV, V_RV, ...
       P1_0_LV, P1_1_LV, P1_2_LV ,P2_0_LV, P2_1_LV, P2_2_LV, P3_0_LV, P3_1_LV, P3_2_LV, N_LV, U_NR_LV,...
       P1_0_SEP,P1_1_SEP,P1_2_SEP,P2_0_SEP,P2_1_SEP,P2_2_SEP,P3_0_SEP,P3_1_SEP,P3_2_SEP,N_SEP,U_NR_SEP,...
       P1_0_RV, P1_1_RV, P1_2_RV, P2_0_RV, P2_1_RV, P2_2_RV, P3_0_RV, P3_1_RV, P3_2_RV, N_RV, U_NR_RV,...
       V_SV, V_PV ,V_SA ,V_PA, V_Ao]';

% Triseg equations to get initial geometry parameters
opts = optimset('MaxFunEvals',10000,'MaxIter',1000);
x_triseg = fsolve(@TrisegEquations,init(1:4),opts,Vw_LV,Vw_SEP,Vw_RV,SL_LV,SL_SEP,SL_RV,V_LV,V_RV,Amref_LV,Amref_SEP,Amref_RV); % Solve triseg equations

% Run cardiovascular mechanics model
init(1:4) = x_triseg;
options = odeset('Mass',M,'MassSingular','yes','RelTol',1e-7,'AbsTol',1e-7,'MaxStep',stim_period/200);
[t,Y] = ode15s(@dXdT_cardiovascular_mechanics,[0 120*stim_period],init,options,para);

% Store output
t_store{1} = t;
Y_store{1} = Y;

% Assigning solution of ODEs to variables
xm_LV  = Y(:,1); % Maximal axial distance from LV midwall surface to origin (cm)
xm_RV  = Y(:,3); % Maximal axial distance from LV midwall surface to origin (cm)
ym     = Y(:,4); % Radius of midwall junction circle (cm)
SL_LV  = Y(:,5); % LV sarcomere length (um)
SL_SEP = Y(:,6); % Septal sarcomere length (um)
SL_RV  = Y(:,7); % RV sarcomere length (um)
V_LV   = Y(:,8); % LV volume (mL)
V_RV   = Y(:,9); % RV volume (mL)

P1_0_LV = Y(:,10); % 0th moment state A1, LV
P1_1_LV = Y(:,11); % 1st moment state A1, LV
P1_2_LV = Y(:,12); % 2nd moment state A1, LV
P2_0_LV = Y(:,13); % 0th moment state A2, LV
P2_1_LV = Y(:,14); % 1st moment state A2, LV
P2_2_LV = Y(:,15); % 2nd moment state A2, LV
P3_0_LV = Y(:,16); % 0th moment state A3, LV
P3_1_LV = Y(:,17); % 1st moment state A3, LV
P3_2_LV = Y(:,18); % 2nd moment state A3, LV
N_LV    = Y(:,19); % Non-permissive state, LV
U_NR_LV = Y(:,20); % Non-super-relaxed state, LV

% Store output
A1_store{1} = P1_0_LV;
A2_store{1} = P2_0_LV;
A3_store{1} = P3_0_LV;
N_store{1} = N_LV;
SRX_store{1} = U_NR_LV;

P3_0_SEP = Y(:,27); % 0th moment state A3, SEP
P3_1_SEP = Y(:,28); % 1st moment state A3, SEP
P3_2_SEP = Y(:,29); % 2nd moment state A3, SEP

V_SV   = Y(:,43); % Volume of systemic veins (mL)
V_PV   = Y(:,44); % Volume of pulmonary veins (mL)
V_SA   = Y(:,45); % Volume of systemic arterys (mL)
V_PA   = Y(:,46); % Volume of pulmonary arterys (mL)
V_Ao   = Y(:,47); % Volume of proximal aorta (mL)

% Store output 
V_SV_store{1} = V_SV;
V_PV_store{1} = V_PV;
V_SA_store{1} = V_SA;
V_PA_store{1} = V_PA;
V_Ao_store{1} = V_Ao;

%% Calculate
% Pulmonary Pressures
P_PV = V_PV/C_PV; % Pulmonary venous pressure (mmHg)
P_SV = V_SV/C_SV; % Systemic venous pressure (mmHg)
P_PA = V_PA/C_PA; % Pulmonary arterial pressure (mmHg)
P_SA = V_SA/C_SA; % Systemic arterial pressure (mmHg)

% Store output
P_PV_store{1} = P_PV;
P_SV_store{1} = P_SV;
P_PA_store{1} = P_PA;
P_SA_store{1} = P_SA;

Am_LV = pi*(xm_LV.^2 + ym.^2); % LV midwall surface area (cm^2)
Am_RV = pi*(xm_RV.^2 + ym.^2); % RV midwall surface area (cm^2)
Cm_LV = 2*xm_LV./(xm_LV.^2 + ym.^2); % Curvature of midwall surface, LV
Cm_RV = 2*xm_RV./(xm_RV.^2 + ym.^2); % Curvature of midwall surface, RV
z_LV = 3*Cm_LV.*Vw_LV./(2*Am_LV);
z_RV = 3*Cm_RV.*Vw_RV./(2*Am_RV);

epsf_LV = (1/2)*log(Am_LV./Amref_LV) - (1/12)*z_LV.^2 - 0.019*z_LV.^4; % Fiber strain
epsf_RV = (1/2)*log(Am_RV./Amref_RV) - (1/12)*z_RV.^2 - 0.019*z_RV.^4; % Fiber strain
SLo_LV = Lsref*exp(epsf_LV); % LV sarcomere length (um)
SLo_RV = Lsref*exp(epsf_RV); % RV sarcomere length (um)

% Total forces
sigmaf_LV = -kse*(SL_LV - SLo_LV); % LV force
sigmaf_RV = -kse*(SL_RV - SLo_RV); % RV force
sigmaf_LV_store{1} = sigmaf_LV; % Store LV force

% Equilibrium of forces at junction circle
Tm_LV = (Vw_LV.*sigmaf_LV./(2*Am_LV)).*(1 + z_LV.^2/3 + z_LV.^4/5); % Tension in LV midwall
Tm_RV = (Vw_RV.*sigmaf_RV./(2*Am_RV)).*(1 + z_RV.^2/3 + z_RV.^4/5); % Tension in RV midwall
sinalpha_LV = 2*xm_LV.*ym./(xm_LV.^2 + ym.^2);
sinalpha_RV = 2*xm_RV.*ym./(xm_RV.^2 + ym.^2);
Tx_LV = Tm_LV.*sinalpha_LV; % Axial component of tension
Tx_RV = Tm_RV.*sinalpha_RV; % Radial component of tension

% Ventricular pressure
ptrans_LV = 2*Tx_LV./ym; % LV transmural pressure (mmHg)
ptrans_RV = 2*Tx_RV./ym; % RV transmural pressure (mmHg)
P_LV = -ptrans_LV; % LV pressure (mmHg)
P_RV = ptrans_RV; % RV pressure (mmHg)

% Ao valve closed equations
P_Ao_closed = (C_SA*R_Ao*R_SA*V_Ao + C_SA*R_Ao*R_tSA*V_Ao + C_SA*R_SA*R_tSA*V_Ao + C_Ao*R_SA*R_tAo*V_SA + C_Ao*C_SA*P_SV*R_tSA*R_tAo)/(C_Ao*C_SA*(R_Ao*R_SA + R_Ao*R_tSA + R_SA*R_tSA + R_SA*R_tAo + R_tSA*R_tAo));

% Ao valve open equations 
P_Ao_open = (C_SA*R_Ao*R_SA*R_AV*V_Ao + C_SA*R_Ao*R_tSA*R_AV*V_Ao + C_SA*R_SA*R_tSA*R_AV*V_Ao + C_Ao*R_SA*R_tAo*R_AV*V_SA + C_Ao*C_SA*P_LV*R_Ao*R_SA*R_tAo + C_Ao*C_SA*P_LV*R_Ao*R_tSA*R_tAo + C_Ao*C_SA*P_LV*R_SA*R_tSA*R_tAo + C_Ao*C_SA*P_SV*R_tSA*R_tAo*R_AV)/(C_Ao*C_SA*(R_Ao*R_SA*R_tAo + R_Ao*R_SA*R_AV + R_Ao*R_tSA*R_tAo + R_Ao*R_tSA*R_AV + R_SA*R_tSA*R_tAo + R_SA*R_tSA*R_AV + R_SA*R_tAo*R_AV + R_tSA*R_tAo*R_AV));

P_Ao = P_Ao_open.*(P_LV>P_Ao_open) + P_Ao_closed.*(P_LV<=P_Ao_open); % Aortic pressure (mmHg)
   
P_Ao_store{1} = P_Ao; % Store output

SV_LV_sim = max(V_LV) - min(V_LV); % Stroke volume (mL/beat)
EF_LV_sim = SV_LV_sim/max(V_LV) * 100; % Ejection fraction (%)
    
edLV_sim =  max(1e3*V_LV); % End-diastolic LV volume (mL)
esLV_sim =  min(1e3*V_LV); % End-systolic LV volume (mL)
edRV_sim =  max(1e3*V_RV); % End-diastolic RV volume (mL)
esRV_sim =  min(1e3*V_RV); % End-systolic RV volume (mL)

g2_LV = (MgATP_cyto/K_T)/(1.0 + MgATP_cyto/K_T + MgADP_cyto/K_D);
k3_LV = k3*g2_LV;
g2_SEP = (MgATP_cyto/K_T)/(1.0 + MgATP_cyto/K_T + MgADP_cyto/K_D);
k3_SEP = k3*g2_SEP;
f_alpha3o_LV  = (P3_0_LV + alpha3*(s3*s3*P3_0_LV + 2.0*s3*P3_1_LV + P3_2_LV)); 
f_alpha3o_SEP = (P3_0_SEP + alpha3*(s3*s3*P3_0_SEP + 2.0*s3*P3_1_SEP + P3_2_SEP)); 

ti = 0:0.00001:stim_period;
MAP = mean(interp1(t,P_Ao,ti)); % Mean arterial pressure (mmHg)

r_LV  = interp1(t,k3_LV*f_alpha3o_LV,ti);
r_SEP = interp1(t,k3_SEP*f_alpha3o_SEP,ti);

% LV XB turnover rate
Vw_LV_W = (2/3)*LVW/1000;
Vw_SEP_W= (1/3)*LVW/1000;
rate_of_XB_turnover_ave = (Vw_LV_W*mean(r_LV) + Vw_SEP_W*mean(r_SEP))/(Vw_LV_W + Vw_SEP_W);

% ATPase rate
ATP_ase_mechanics_Averge_LV_SEP = (1.3/rate_of_XB_turnover_mean_sham)*rate_of_XB_turnover_ave; 
tune_ATPase_LV =  ATP_ase_mechanics_Averge_LV_SEP * (1/ 0.6801) *1.0e-3; % ATP hydrolysis rate: M/s/L cytosol

% Store output
ATPase_store(1) = tune_ATPase_LV;
P_LV_store{1} = P_LV;
V_LV_store{1} = V_LV;
shortening_store{1} = Y(:,5); 


%% Run coupled cardiovascular mechanics/energetics model
% Initialize
p = 2;
beat = 1;

for j = 122:420 % Run for 300 beats
tspan_beat = [(stim_period)*(j-1) (stim_period)*j];

% Run energetics model
if rem(j,3) == 0 % Run every 3 beats (for model stability)
energtics_output  = energetics_model(TAN, CRtot, TEP, Ox_capacity, tune_ATPase_LV);

MgATP_cyto = energtics_output(1); % Cytosolic MgATP concentration (M·(L cytosol water)^-1)
MgADP_cyto = energtics_output(2); % Cytosolic MgADP concentration (M·(L cytosol water)^-1)
fPi_cytoplasm = energtics_output(3); % Cytosolic unchelated Pi concentration (M·(L cytosol water)^-1)
MVO2_tissue = energtics_output(5); % Oxygen consumption rate (uM·min^-1·(g tissue)^-1)
dGrATPase = energtics_output(6); % ATP hydrolysis free energy, kJ/mol

PCrATP =  energtics_output(7); % Creatine phoosphate ATP ratio, unitless
ATP_cyto = energtics_output(8); % Cytosolic total ATP concentration (M·(L cytosol water)^-1)
ADP_cyto = energtics_output(9); % Cytosolic total ADP concentration (M·(L cytosol water)^-1)
Pi_cyto = energtics_output(10)*1000; % Cytosolic total Pi concentration (M·(L cytosol water)^-1)
end

% Store output
ATP_store(p) = MgATP_cyto;
ADP_store(p) = MgADP_cyto*1000;
Pi_store(p) = Pi_cyto;
MVO2_tissue_store(p) = MVO2_tissue;
PCrATP_store(p) = PCrATP;

% Run mechanics model
init = Y(end,:);
para = [ka, kd, k1, k_1, k2, k_2, k3, krecruit, k_on ,k_off, k_coop, visc, kstiff1, kstiff2, kSR, k_SR, k_passive,... 
        Lthin, Lthick, Lbare, dr, SLcollagen, PConcollagen, PExpcollagen, Lsref, L_rest_pas, K_T, K_D, K_Pi, alpha1,...
        alpha2, alpha3, s3, scale, kse, Vw_LV, Vw_SEP, Vw_RV, Amref_LV, Amref_SEP, Amref_RV, C_Ao, C_SA, C_SV, C_PV,...
        C_PA, R_Ao, R_SA, R_PA, R_SV, R_PV, R_vlv, R_AV, R_tAo, R_tSA, MgATP_cyto, MgADP_cyto, Pi_cyto, MgATP_cyto,... 
        MgADP_cyto, Pi_cyto, MgATP_cyto, MgADP_cyto, Pi_cyto, flag_swap_metabolite, Ca_flag, r, stim_period, dATP];
[t,Y] = ode15s(@dXdT_cardiovascular_mechanics, tspan_beat, init, options, para);

% Store output
t_store{p} = t;
Y_store{p} = Y;

xm_LV  = Y(:,1); % Maximal axial distance from LV midwall surface to origin (cm)
xm_RV  = Y(:,3); % Maximal axial distance from LV midwall surface to origin (cm)
ym     = Y(:,4); % Radius of midwall junction circle (cm)
SL_LV  = Y(:,5); % LV sarcomere length (um)
SL_SEP = Y(:,6); % Septal sarcomere length (um)
SL_RV  = Y(:,7); % RV sarcomere length (um)
V_LV   = Y(:,8); % LV volume (mL)
V_RV   = Y(:,9); % RV volume (mL)

P1_0_LV = Y(:,10); % 0th moment state A1, LV
P1_1_LV = Y(:,11); % 1st moment state A1, LV
P1_2_LV = Y(:,12); % 2nd moment state A1, LV
P2_0_LV = Y(:,13); % 0th moment state A2, LV
P2_1_LV = Y(:,14); % 1st moment state A2, LV
P2_2_LV = Y(:,15); % 2nd moment state A2, LV
P3_0_LV = Y(:,16); % 0th moment state A3, LV
P3_1_LV = Y(:,17); % 1st moment state A3, LV
P3_2_LV = Y(:,18); % 2nd moment state A3, LV
N_LV    = Y(:,19); % Non-permissive state, LV
U_NR_LV = Y(:,20); % Non-super-relaxed state, LV

% Store output
A1_store{p} = P1_0_LV;
A2_store{p} = P2_0_LV;
A3_store{p} = P3_0_LV;
N_store{p} = N_LV;
SRX_store{p} = U_NR_LV;

P3_0_SEP = Y(:,27); % 0th moment state A3, SEP
P3_1_SEP = Y(:,28); % 1st moment state A3, SEP
P3_2_SEP = Y(:,29); % 2nd moment state A3, SEP

V_SV   = Y(:,43); % Volume of systemic veins (mL)
V_PV   = Y(:,44); % Volume of pulmonary veins (mL)
V_SA   = Y(:,45); % Volume of systemic arterys (mL)
V_PA   = Y(:,46); % Volume of pulmonary arterys (mL)
V_Ao   = Y(:,47); % Volume of proximal aorta (mL)

V_SV_store{p} = V_SV;
V_PV_store{p} = V_PV;
V_SA_store{p} = V_SA;
V_PA_store{p} = V_PA;
V_Ao_store{p} = V_Ao;

%% Calculate
% Pulmonary Pressures
P_PV = V_PV/C_PV; % Pulmonary venous pressure (mmHg)
P_SV = V_SV/C_SV; % Systemic venous pressure (mmHg)
P_PA = V_PA/C_PA; % Pulmonary arterial pressure (mmHg)
P_SA = V_SA/C_SA; % Systemic arterial pressure (mmHg)

P_SV_store{p} = P_SV;
P_PV_store{p} = P_PV;
P_SA_store{p} = P_SA;
P_PA_store{p} = P_PA;

Am_LV = pi*(xm_LV.^2 + ym.^2); % LV midwall surface area (cm^2)
Am_RV = pi*(xm_RV.^2 + ym.^2); % RV midwall surface area (cm^2)
Cm_LV = 2*xm_LV./(xm_LV.^2 + ym.^2); % Curvature of midwall surface, LV
Cm_RV = 2*xm_RV./(xm_RV.^2 + ym.^2); % Curvature of midwall surface, RV
z_LV = 3*Cm_LV.*Vw_LV./(2*Am_LV);
z_RV = 3*Cm_RV.*Vw_RV./(2*Am_RV);

epsf_LV = (1/2)*log(Am_LV./Amref_LV) - (1/12)*z_LV.^2 - 0.019*z_LV.^4; % Fiber strain
epsf_RV = (1/2)*log(Am_RV./Amref_RV) - (1/12)*z_RV.^2 - 0.019*z_RV.^4; % Fiber strain
SLo_LV = Lsref*exp(epsf_LV); % LV sarcomere length (um)
SLo_RV = Lsref*exp(epsf_RV); % RV sarcomere length (um)

% Total forces
sigmaf_LV = -kse*(SL_LV - SLo_LV); % LV force
sigmaf_RV = -kse*(SL_RV - SLo_RV); % RV force
sigmaf_LV_store{p} = sigmaf_LV; % Store LV force

% Equilibrium of forces at junction circle
Tm_LV = (Vw_LV.*sigmaf_LV./(2*Am_LV)).*(1 + z_LV.^2/3 + z_LV.^4/5); % Tension in LV midwall
Tm_RV = (Vw_RV.*sigmaf_RV./(2*Am_RV)).*(1 + z_RV.^2/3 + z_RV.^4/5); % Tension in RV midwall
sinalpha_LV = 2*xm_LV.*ym./(xm_LV.^2 + ym.^2);
sinalpha_RV = 2*xm_RV.*ym./(xm_RV.^2 + ym.^2);
Tx_LV = Tm_LV.*sinalpha_LV; % Axial component of tension
Tx_RV = Tm_RV.*sinalpha_RV; % Radial component of tension

% Ventricular pressure
ptrans_LV = 2*Tx_LV./ym; % LV transmural pressure (mmHg)
ptrans_RV = 2*Tx_RV./ym; % RV transmural pressure (mmHg)
P_LV = -ptrans_LV; % LV pressure (mmHg)
P_RV = ptrans_RV; % RV pressure (mmHg)

% Ao valve closed equations
P_Ao_closed = (C_SA*R_Ao*R_SA*V_Ao + C_SA*R_Ao*R_tSA*V_Ao + C_SA*R_SA*R_tSA*V_Ao + C_Ao*R_SA*R_tAo*V_SA + C_Ao*C_SA*P_SV*R_tSA*R_tAo)/(C_Ao*C_SA*(R_Ao*R_SA + R_Ao*R_tSA + R_SA*R_tSA + R_SA*R_tAo + R_tSA*R_tAo));

% Ao valve open equations 
P_Ao_open = (C_SA*R_Ao*R_SA*R_AV*V_Ao + C_SA*R_Ao*R_tSA*R_AV*V_Ao + C_SA*R_SA*R_tSA*R_AV*V_Ao + C_Ao*R_SA*R_tAo*R_AV*V_SA + C_Ao*C_SA*P_LV*R_Ao*R_SA*R_tAo + C_Ao*C_SA*P_LV*R_Ao*R_tSA*R_tAo + C_Ao*C_SA*P_LV*R_SA*R_tSA*R_tAo + C_Ao*C_SA*P_SV*R_tSA*R_tAo*R_AV)/(C_Ao*C_SA*(R_Ao*R_SA*R_tAo + R_Ao*R_SA*R_AV + R_Ao*R_tSA*R_tAo + R_Ao*R_tSA*R_AV + R_SA*R_tSA*R_tAo + R_SA*R_tSA*R_AV + R_SA*R_tAo*R_AV + R_tSA*R_tAo*R_AV));

P_Ao = P_Ao_open.*(P_LV>P_Ao_open) + P_Ao_closed.*(P_LV<=P_Ao_open); % Aortic pressure (mmHg)
   
P_Ao_store{p} = P_Ao;

SV_LV_sim = max(V_LV) - min(V_LV); % Stroke volume (mL/beat)
EF_LV_sim = SV_LV_sim/max(V_LV) * 100; % Ejection fraction (%)
    
edLV_sim =  max(1e3*V_LV); % End-diastolic LV volume (mL)
esLV_sim =  min(1e3*V_LV); % End-systolic LV volume (mL)
edRV_sim =  max(1e3*V_RV); % End-diastolic RV volume (mL)
esRV_sim =  min(1e3*V_RV); % End-systolic RV volume (mL)

g2_LV = (MgATP_cyto/K_T)/(1.0 + MgATP_cyto/K_T + MgADP_cyto/K_D);
k3_LV = k3*g2_LV;
g2_SEP = (MgATP_cyto/K_T)/(1.0 + MgATP_cyto/K_T + MgADP_cyto/K_D);
k3_SEP = k3*g2_SEP;
f_alpha3o_LV  = (P3_0_LV + alpha3*(s3*s3*P3_0_LV + 2.0*s3*P3_1_LV + P3_2_LV)); 
f_alpha3o_SEP = (P3_0_SEP + alpha3*(s3*s3*P3_0_SEP + 2.0*s3*P3_1_SEP + P3_2_SEP)); 

ti = stim_period*(j-1):0.00001:stim_period*j;
MAP = mean(interp1(t,P_Ao,ti)); % Mean arterial pressure (mmHg)

r_LV  = interp1(t,k3_LV*f_alpha3o_LV,ti);
r_SEP = interp1(t,k3_SEP*f_alpha3o_SEP,ti);

% LV XB turnover rate
Vw_LV_W = (2/3)*LVW/1000;
Vw_SEP_W= (1/3)*LVW/1000;
rate_of_XB_turnover_ave = (Vw_LV_W*mean(r_LV) + Vw_SEP_W*mean(r_SEP))/(Vw_LV_W + Vw_SEP_W);

% ATPase rate
ATP_ase_mechanics_Averge_LV_SEP = (1.3/rate_of_XB_turnover_mean_sham)*rate_of_XB_turnover_ave; 
tune_ATPase_LV =  ATP_ase_mechanics_Averge_LV_SEP * (1/ 0.6801) *1.0e-3; % ATP hydrolysis rate: M/s/L cytosol

% Store output
ATPase_store(p) = tune_ATPase_LV;
P_LV_store{p} = P_LV;
V_LV_store{p} = V_LV;
shortening_store{p} = Y(:,5); 

beat = p % Output which beat the simulation is on
p = p+1;
end


%% Finalize output
% Initialize
k = 1;
V_LV_f = [];
P_LV_f = [];
t_f = [];
shortening_f = [];
sigmaf_LV_f = [];
A1_f = [];
A2_f = [];
A3_f = [];
N_f = [];
SRX_f = [];
P_SV_f = [];
P_PV_f = [];
P_SA_f = [];
P_PA_f = [];
P_Ao_f = [];
V_SV_f = [];
V_PV_f = [];
V_SA_f = [];
V_PA_f = [];
V_Ao_f = [];

while k < p-1
    V_LV1 = V_LV_store{k};
    V_LV2 = V_LV_store{k+1};
    V_LV_new = vertcat(V_LV1(2:end),V_LV2(2:end));
    V_LV_f  =  vertcat(V_LV_f,V_LV_new);
    
    P_LV1 = P_LV_store{k};
    P_LV2 = P_LV_store{k+1};
    P_LV_new = vertcat(P_LV1(2:end),P_LV2(2:end));
    P_LV_f  =  vertcat(P_LV_f,P_LV_new);
    
    t_1 = t_store{k};
    t_2 = t_store{k+1};
    t_new = vertcat(t_1(2:end),t_2(2:end));
    t_f = vertcat(t_f,t_new);
    
    shortening_1 = shortening_store{k};
    shortening_2 = shortening_store{k+1};
    shortening_new = vertcat(shortening_1(2:end),shortening_2(2:end));
    shortening_f  =  vertcat(shortening_f,shortening_new);
    
    sigmaf_LV_1 = sigmaf_LV_store{k};
    sigmaf_LV_2 = sigmaf_LV_store{k+1};
    sigmaf_LV_new = vertcat(sigmaf_LV_1(2:end),sigmaf_LV_2(2:end));
    sigmaf_LV_f  =  vertcat(sigmaf_LV_f,sigmaf_LV_new);
    
    A1_1 = A1_store{k};
    A1_2 = A1_store{k+1};
    A1_new = vertcat(A1_1(2:end),A1_2(2:end));
    A1_f  =  vertcat(A1_f,A1_new);
    
    A2_1 = A2_store{k};
    A2_2 = A2_store{k+1};
    A2_new = vertcat(A2_1(2:end),A2_2(2:end));
    A2_f  =  vertcat(A2_f,A2_new);
    
    A3_1 = A3_store{k};
    A3_2 = A3_store{k+1};
    A3_new = vertcat(A3_1(2:end),A3_2(2:end));
    A3_f  =  vertcat(A3_f,A3_new);
    
    N_1 = N_store{k};
    N_2 = N_store{k+1};
    N_new = vertcat(N_1(2:end),N_2(2:end));
    N_f  =  vertcat(N_f,N_new);
    
    SRX_1 = SRX_store{k};
    SRX_2 = SRX_store{k+1};
    SRX_new = vertcat(SRX_1(2:end),SRX_2(2:end));
    SRX_f  =  vertcat(SRX_f,SRX_new);
    
    P_SV_1 = P_SV_store{k};
    P_SV_2 = P_SV_store{k+1};
    P_SV_new = vertcat(P_SV_1(2:end),P_SV_2(2:end));
    P_SV_f  =  vertcat(P_SV_f,P_SV_new);
    
    P_PV_1 = P_PV_store{k};
    P_PV_2 = P_PV_store{k+1};
    P_PV_new = vertcat(P_PV_1(2:end),P_PV_2(2:end));
    P_PV_f  =  vertcat(P_PV_f,P_PV_new);
    
    P_SA_1 = P_SA_store{k};
    P_SA_2 = P_SA_store{k+1};
    P_SA_new = vertcat(P_SA_1(2:end),P_SA_2(2:end));
    P_SA_f  =  vertcat(P_SA_f,P_SA_new);
    
    P_PA_1 = P_PA_store{k};
    P_PA_2 = P_PA_store{k+1};
    P_PA_new = vertcat(P_PA_1(2:end),P_PA_2(2:end));
    P_PA_f  =  vertcat(P_PA_f,P_PA_new);
    
    P_Ao_1 = P_Ao_store{k};
    P_Ao_2 = P_Ao_store{k+1};
    P_Ao_new = vertcat(P_Ao_1(2:end),P_Ao_2(2:end));
    P_Ao_f  =  vertcat(P_Ao_f,P_Ao_new);
    
    V_SV_1 = V_SV_store{k};
    V_SV_2 = V_SV_store{k+1};
    V_SV_new = vertcat(V_SV_1(2:end),V_SV_2(2:end));
    V_SV_f  =  vertcat(V_SV_f,V_SV_new);
    
    V_PV_1 = V_PV_store{k};
    V_PV_2 = V_PV_store{k+1};
    V_PV_new = vertcat(V_PV_1(2:end),V_PV_2(2:end));
    V_PV_f  =  vertcat(V_PV_f,V_PV_new);
    
    V_SA_1 = V_SA_store{k};
    V_SA_2 = V_SA_store{k+1};
    V_SA_new = vertcat(V_SA_1(2:end),V_SA_2(2:end));
    V_SA_f  =  vertcat(V_SA_f,V_SA_new);
    
    V_PA_1 = V_PA_store{k};
    V_PA_2 = V_PA_store{k+1};
    V_PA_new = vertcat(V_PA_1(2:end),V_PA_2(2:end));
    V_PA_f  =  vertcat(V_PA_f,V_PA_new);
    
    V_Ao_1 = V_Ao_store{k};
    V_Ao_2 = V_Ao_store{k+1};
    V_Ao_new = vertcat(V_Ao_1(2:end),V_Ao_2(2:end));
    V_Ao_f  =  vertcat(V_Ao_f,V_Ao_new);
      
    k = k+2;
end

beats  = 1:420;
ATPase_orig = zeros(1,120);
for i = 1:length(ATPase_orig)
    ATPase_orig(i) = ATPase_store(1);
end
ATPase_full = horzcat(ATPase_orig,ATPase_store);

MVO2_tissue_orig = zeros(1,120);
for i = 1:length(MVO2_tissue_orig)
    MVO2_tissue_orig(i) = MVO2_tissue_store(1);
end
MVO2_tissue_full = horzcat(MVO2_tissue_orig,MVO2_tissue_store);

PCrATP_orig = zeros(1,120);
for i = 1:length(PCrATP_orig)
    PCrATP_orig(i) = PCrATP_store(1);
end
PCrATP_full = horzcat(PCrATP_orig,PCrATP_store);

ATP_orig = zeros(1,120);
for i = 1:length(ATP_orig)
    ATP_orig(i) = ATP_store(1);
end
ATP_full = horzcat(ATP_orig,ATP_store);

ADP_orig = zeros(1,120);
for i = 1:length(ADP_orig)
    ADP_orig(i) = ADP_store(1);
end
ADP_full = horzcat(ADP_orig,ADP_store);

Pi_orig = zeros(1,120);
for i = 1:length(Pi_orig)
    Pi_orig(i) = Pi_store(1);
end
Pi_full = horzcat(Pi_orig,Pi_store);

% Get just last beat
last_beat = 419; 
last_beat_t = last_beat*stim_period;
for i = 1:length(t_f)
    if t_f(i) > last_beat_t
        idx_last = i;
        break
    end
end


%% Final calculations
t_f_last = t_f(idx_last:end);
P_LV_last = P_LV_f(idx_last:end);
V_LV_last = V_LV_f(idx_last:end);
shortening_last = shortening_f(idx_last:end);
force_last = sigmaf_LV_f(idx_last:end);
A1_last = A1_f(idx_last:end);
A2_last = A2_f(idx_last:end);
A3_last = A3_f(idx_last:end);
N_last = N_f(idx_last:end);
SRX_last = SRX_f(idx_last:end);
P_SV_last = P_SV_f(idx_last:end);
P_PV_last = P_PV_f(idx_last:end);
P_SA_last = P_SA_f(idx_last:end);
P_PA_last = P_PA_f(idx_last:end);
P_Ao_last = P_Ao_f(idx_last:end);
V_SV_last = V_SV_f(idx_last:end);
V_PV_last = V_PV_f(idx_last:end);
V_SA_last = V_SA_f(idx_last:end);
V_PA_last = V_PA_f(idx_last:end);
V_Ao_last = V_Ao_f(idx_last:end);

V_LV_store_F{r} = V_LV_last;
P_LV_store_F{r} = P_LV_last;
t_store_F{r} = t_f_last;
sigmaf_LV_store_F{r} = force_last;
shortening_store_F{r} = shortening_last;
A1_store_F{r} = A1_last;
A2_store_F{r} = A2_last;
A3_store_F{r} = A3_last;
N_store_F{r} = N_last;
P_SV_store_F{r} = P_SV_last;
P_PV_store_F{r} = P_PV_last;
P_SA_store_F{r} = P_SA_last;
P_PA_store_F{r} = P_PA_last;
P_Ao_store_F{r} = P_Ao_last;
V_SV_store_F{r} = V_SV_last;
V_PV_store_F{r} = V_PV_last;
V_SA_store_F{r} = V_SA_last;
V_PA_store_F{r} = V_PA_last;
V_Ao_store_F{r} = V_Ao_last;
SRX_store_F{r} = SRX_last;

max_force{r} = max(force_last); % Max force (kPa)
FS{r} = ((min(shortening_last)-max(shortening_last))/max(shortening_last))*100; % Fractional shortening (%)
CO{r} = (max(V_LV_last)-min(V_LV_last))*HR; % Cardiac output (mL/min)
EF{r} = (max(V_LV_last)-min(V_LV_last))/max(V_LV_last); % Ejection fracion (%)

max_P = max(P_LV_last); % End-systolic pressure (mmHg)
min_P = min(P_LV_last); % End-diastolic pressure (mmHg)
LVDP{r} = max_P-min_P; % LV Developed pressure (mmHg)

% % Calculate rate of pressure development and decline
for i = 1:length(P_LV_last)
    if P_LV_last(i) >= max_P
        idx_max = i;
        break
    end
end

up_slope = P_LV_last(1:idx_max);
t_up = t_f_last(1:idx_max);
down_slope = P_LV_last(idx_max:end);
t_down = t_f_last(idx_max:end);

plus_deriv = gradient(up_slope,0.001);
minus_deriv = gradient(down_slope,0.001);

plus_dPdt_F{r} = max(plus_deriv); % Rate of pressure development (mmHg/mL)
minus_dPdt_F{r} = min(minus_deriv); % Rate of pressure decline (mmHg/mL)

work_rate_F{r} = SV_LV_sim*MAP*HR/60; % Work rate (mmHg*mL/s)
ATP_F{r} = ATP_store(300); % ATP (mM)
ADP_F{r} = ADP_store(300); % ADP (mM)
Pi_F{r} = Pi_store(300); % Pi (mM)
MVO2_F{r} = MVO2_tissue_store(300); % Oxygen consumption (uM O2/min/g tissue)
PCrATP_F{r} = PCrATP_store(300); % Creatine phosphate ATP ratio (unitless)
XB_turnover_F{r} = rate_of_XB_turnover_ave; % XB turnover rate 
ATPase_F{r} =  ATPase_store(300); % ATPase rate (M/ms/L)
work_per_beat = ((SV_LV_sim*MAP*HR/60)/HR)*60; % Work rate (mmHg*mL/s)
efficiency_F{r} = (work_per_beat/ATPase_store(300))*1000; % Efficiency (mL^2*mmHg*ms/M)
end


%% Plot
if suppress_plotting == 1
figure
hold on
plot(V_LV_store_F{1},P_LV_store_F{1}-47,'linewidth',3)
xlabel('Volume (mL)')
ylabel('Pressure (mmHg)')
set(gca, 'fontsize', 14)  
end

end
