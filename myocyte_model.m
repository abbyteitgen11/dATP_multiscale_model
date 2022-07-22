%%%% Main code to run myocyte model
% Adapted by Abigail Teitgen based on code from: 

% Tewari, S. G., Bugenhagen, S. M., Palmer, B. M, Beard,
% D. A. (2016). Dynamics of corss-bridge cycling, ATP hydrolysis, force
% generation, and deformation in cardiac muscle. J. Mol. Cell Cardiol. 96:
% 11-25.

% Tewari,S. G., Bugenhagen, S. M., Vinnakota, K. C., Rice, J. J., Janssen,
% M. L., Beard, D. A. (2016). Influence of metabolic dysfuntion on cardiac
% mechanics in decompensated hypertrophy and heart failure. J. Mol.
% Cell Cardiol. 94: 162-175.

% Lopez, R. Marzban, B., Gao, X. Lauinger, E., Van den Bergh, F.,
% Whitesall, S. E., Converso-Baran, K., Burant, C. F., Michele, D. E.,
% Beard, D. A. (2020). Impaired myocardial energetics causes mechanical dysfunction
% in decompensated failing hearts. Function, 1(2): zqaa018

% Note: code requires parallel computing toolbox for use of parallel for
% loops (alternatively, this can be removed and normal for loops can be
% used)


function [T_final_XB, force_final, idx_XB, Shortening_final, SS_Ftotal_fpca, Ftotal_ktr, t_ktr] = myocyte_model(Ca_value, XB, Ktr, plotting, dATP_percent, ka_i, kd_i, k1_i, k_1_i, k2_i, k_2_i, k3_i, krecruit_i, k_on_i, k_off_i, k_coop_i)

tic 
     
%% Flags
% Calcium
Ca_flag = Ca_value; % 0 = ATP, 1 = dATP

% Protocol
XB_protocol = XB; % 0 = force pCa, 1 = Ktr, 2 = twitch, 3 = all

% pCa
Ktr_protocol = Ktr; % 0 = pCa 4.0, 1 = pCa 4.5, 2 = pCa 5.0, 3 = pCa 5.5, 4 = pCa 6.0, 5 = pCa 6.5, 6 = pCa 7.0

% Suppress plotting and calculation outputs
suppress_plotting = plotting; % 0 = do not plot, 1 = plot

% Percent dATP
dATP = dATP_percent; % dATP fraction
ATP = 100 - dATP; % ATP fraction

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
k_on = 85*(ATP/100) + k_on_i*(dATP/100); % Rate constant of Ca2+ binding to troponin C (uM^-1s^-1)
k_off = 900*(ATP/100) + k_off_i*(dATP/100); % Rate constant of Ca2+ unbinding from troponin C (s^-1)
k_coop = k_coop_i; % Strength of thin filament cooperativity

visc = 0.01; % Viscosity (mmHg*s/um)
kstiff1 = 1.8219e+04; % Stiffness constant due to myosin-actin interaction (kPa/um)
kstiff2 = 4.7822e+05; % Stiffness constant due to working stroke of XBs (kPa/um)
kSR = 1.3; % Forward rate constant of force dependent super-relaxed state transition (USR to N or P) (s^-1)
k_SR = 30; % Reverse rate constant of force dependent super-relaxed transition (s^-1)
k_passive = 0.1; % Passive stiffness constant (kPa/um)
Lthin = 1.200; % Length of thin filament (nm)
Lthick = 1.670; % Length of thick filament (nm)
Lbare = 0.100; % Bare length of thin filament (nm)
dr = 0.01; % Power-stroke size (um)
SLcollagen = 2.25; % Threshold for collagen activation (um)
PConcollagen = 0.01; % Scale factor for passive force contributed by collagen
PExpcollagen = 70; % um^-1
Lsref = 1.9; % Resting sarcomere length (um)
stim_period = 1; %Hz

% Metabolite concentrations based on mean sham rat (Lopez et al. 2020)
MgATP = 7.8873; % mM
MgADP = 0.0501; % mM
Pi = 1.2308; % mM

%% Run model
%% Force pCa
if XB_protocol == 0 || XB_protocol == 3 % Force pCa or all
flag = 0;
kse = 50000; % Series element stiffness (mmHg/um) 
SL0_fpca = 2.25; % (um) based on experimental protocol (Regnier et al. 2004)

% For input into ODE solver
para = [MgATP, MgADP, Pi, kstiff1, kstiff2, k_passive, kse, k_coop, k_on, k_off, kSR, krecruit, k_SR, ka, kd, k1, k_1, k2, k_2, k3, visc, stim_period];

% Defining time vector
tspan_fpca = 0:0.001:0.3;

% Defining Ca range
Ca_fraction_fpca = [0.1:0.1:100]; % pCa 7 to 4

init_fpca  = [zeros(1,10), SL0_fpca,0]; % Initial conditions for the model
init_fpca(10) = 1; % Setting the initial value for nonpermissible state equal to 1
SS_Ftotal_fpca = zeros(1,length(Ca_fraction_fpca)); % Initialize storage

parfor k = 1:length(Ca_fraction_fpca) % Parallel for loop
    options = odeset('RelTol',1e-3,'AbsTol',1e-6,'MaxStep',5e-3);
    [t_fpca, Y_fpca] = ode15s(@dXdT_myocyte_mechanics, tspan_fpca, init_fpca, options, para, Ca_fraction_fpca(k), flag, Ca_flag);
    
    % Get states
    p1_0_fpca = Y_fpca(:,1);
    p2_0_fpca = Y_fpca(:,4);
    p2_1_fpca = Y_fpca(:,5);
    p3_0_fpca = Y_fpca(:,7);
    p3_1_fpca = Y_fpca(:,8);
    SL_fpca = Y_fpca(:,11);
    N_fpca = Y_fpca(:,10);
       
    % Overlap function
    OV_Zaxis_fpca = min(Lthick/2, SL_fpca/2); % Overlap region closest to Z-axis (nm)
    OV_Mline_fpca = max(SL_fpca/2-(SL_fpca-Lthin), Lbare/2); % Overlap region closest to M-line (nm)
    LOV_fpca = OV_Zaxis_fpca - OV_Mline_fpca; % Length of overlap (nm)
    N_overlap_thick_fpca = LOV_fpca*2/(Lthick - Lbare); % Fraction of thick filament overlap

    % Active force
    B_process_fpca = kstiff2 * dr * p3_0_fpca;   % Force due to XB ratcheting
    C_process_fpca = kstiff1 * (p2_1_fpca + p3_1_fpca); % Force due to stretching of XBs
    F_XB_fpca = N_overlap_thick_fpca.*(B_process_fpca + C_process_fpca); % Active force

    % (Linear) passive force model
    sigma_collagen_fpca  = PConcollagen*(exp(PExpcollagen*(SL_fpca - SLcollagen)) - 1).*(SL_fpca > SLcollagen); % Collagen force
    F_passive_fpca  = k_passive*(SL_fpca/2-Lsref/2)   + sigma_collagen_fpca; % Passive Force 
    Ftotal_fpca = F_XB_fpca + F_passive_fpca; % Total force
    SS_Ftotal_fpca(k) = Ftotal_fpca(end)/10; % Store force output
end
end

%% Ktr
if XB_protocol == 1 || XB_protocol == 3 % Ktr or all
flag = 1;
kse = 50000; % Series element stiffness (mmHg/um) 
SL0_ktr = 2.25; % (um) based on experimental protocol (Regnier et al. 2004)

% For input into ODE solver
para = [MgATP, MgADP, Pi, kstiff1, kstiff2, k_passive, kse, k_coop, k_on, k_off, kSR, krecruit, k_SR, ka, kd, k1, k_1, k2, k_2, k3, visc, stim_period];

% Defining Ca range
if Ktr_protocol == 0 % pCa 4.0
    Ca_fraction_ktr = 100;

else if Ktr_protocol == 1 % pCa 4.5
    Ca_fraction_ktr = 31.6228;

else if Ktr_protocol == 2 % pCa 5.0
    Ca_fraction_ktr = 10;

else if Ktr_protocol == 3 % pCa 5.5
    Ca_fraction_ktr = 3.1623;

else if Ktr_protocol == 4 % pCa 6.0
    Ca_fraction_ktr = 1;

else if Ktr_protocol == 5 % pCa 6.5
    Ca_fraction_ktr = 0.3162;

else if Ktr_protocol == 6 % pCa 7.0
    Ca_fraction_ktr = 0.1;
end
end
end
end
end
end
end

% 0 to 4.488 s (steady state force development)
tspan_ktr_1 = 0:0.001:4.488;
init_ktr_1 = [zeros(1,10),SL0_ktr,0]; % Initial conditions for the model
options = odeset('RelTol',1e-3,'AbsTol',1e-6,'MaxStep',5e-3);
[t_ktr_1, Y_ktr_1] = ode15s(@dXdT_myocyte_mechanics, tspan_ktr_1, init_ktr_1, options, para, Ca_fraction_ktr, flag, Ca_flag);
init1 = Y_ktr_1(end,:);

% 4.488 to 4.5469 s (shorten to 54% SL0)
init1(11) = 0.54*SL0_ktr; 
tspan_ktr_2 = 4.488:0.001:4.5469;
options = odeset('RelTol',1e-3,'AbsTol',1e-6,'MaxStep',5e-3);
[t_ktr_2, Y_ktr_2] = ode15s(@dXdT_myocyte_mechanics, tspan_ktr_2, init1, options, para, Ca_fraction_ktr, flag, Ca_flag);
init2 = Y_ktr_2(end,:);

% 4.5469 to 6 s (reset to SL0)
init2(11) = SL0_ktr;
tspan_ktr_3 = 4.5469:0.001:6; 
options = odeset('RelTol',1e-3,'AbsTol',1e-6,'MaxStep',5e-3);
[t_ktr_3, Y_ktr_3] = ode15s(@dXdT_myocyte_mechanics, tspan_ktr_3, init2, options, para, Ca_fraction_ktr, flag, Ca_flag);

% Combine output
t_ktr = [t_ktr_1; t_ktr_2; t_ktr_3]; Y_ktr = [Y_ktr_1; Y_ktr_2; Y_ktr_3];

% Get states
p1_0_ktr = Y_ktr(:,1);
p2_0_ktr = Y_ktr(:,4);
p2_1_ktr = Y_ktr(:,5);
p3_0_ktr = Y_ktr(:,7);
p3_1_ktr = Y_ktr(:,8);
SL_ktr = Y_ktr(:,11);
N_ktr = Y_ktr(:,10);

% Overlap function
OV_Zaxis_ktr = min(Lthick/2, SL_ktr/2); % Overlap region closest to Z-axis (nm)
OV_Mline_ktr = max(SL_ktr/2-(SL_ktr-Lthin), Lbare/2); % Overlap region closest to M-line (nm)
LOV_ktr = OV_Zaxis_ktr - OV_Mline_ktr; % Length of overlap (nm)
N_overlap_thick_ktr = LOV_ktr*2/(Lthick - Lbare); % Fraction of thick filament overlap

% Active force
B_process_ktr = kstiff2 * dr * p3_0_ktr;   % Force due to XB ratcheting
C_process_ktr = kstiff1 * (p2_1_ktr + p3_1_ktr);% Force due to stretching of XBs
F_XB_ktr = N_overlap_thick_ktr.*(B_process_ktr + C_process_ktr); % Active force

% (Linear) passive force model
sigma_collagen_ktr  = PConcollagen*(exp(PExpcollagen*(SL_ktr - SLcollagen)) - 1).*(SL_ktr > SLcollagen); % Collagen force
F_passive_ktr  = k_passive*(SL_ktr/2-Lsref/2)   + sigma_collagen_ktr; % Passive force
Ftotal_ktr = F_XB_ktr + F_passive_ktr; % Total force
Ftotal_ktr = Ftotal_ktr/10; % Store force output

end        

if XB_protocol == 2 || XB_protocol == 3 % Twitch or all
flag = 2;
kse = 4.35; % Series element stiffness (mmHg/um) 
SL0_twitch = 1.84; % (um) from Korte et al. 2011 experimental protocol

% For input into ODE solver
para = [MgATP, MgADP, Pi, kstiff1, kstiff2, k_passive, kse, k_coop, k_on, k_off, kSR, krecruit, k_SR, ka, kd, k1, k_1, k2, k_2, k3, visc, stim_period];

beats = 3; % Number of beats
tspan_twitch = 0:0.00001:beats; % Time vector (each beat is 1000 ms)
last_beat = beats*1000-1000; % Point at which last beat starts (for plotting)
y_XB = [zeros(1,10),SL0_twitch,0]; % Initialize
y_XB(10) = 1; % Setting the initial value for nonpermissible state equal to 1
init_twitch = y_XB; % Initial conditions
options = odeset('RelTol',1e-3,'AbsTol',1e-6,'MaxStep',5e-3);
Ca_i = 0;
[t_twitch, Y_twitch] = ode15s(@dXdT_myocyte_mechanics, tspan_twitch, init_twitch, options, para, Ca_i, flag, Ca_flag);

% Get states
SL_twitch = Y_twitch(:,11);
p1_0_twitch = Y_twitch(:,1);
p2_0_twitch = Y_twitch(:,4);
p2_1_twitch = Y_twitch(:,5);
p3_0_twitch = Y_twitch(:,7);
p3_1_twitch = Y_twitch(:,8);
N_twitch = Y_twitch(:,10);

% Overlap function
OV_Zaxis_twitch = min(Lthick/2,SL_twitch/2); % Overlap region closest to Z-axis (nm)
OV_Mline_twitch = max(SL_twitch/2-(SL_twitch-Lthin),Lbare/2); % Overal region closest to M-line (nm)
LOV_twitch = OV_Zaxis_twitch - OV_Mline_twitch; % Length of overlap (nm)
N_overlap_thick_twitch = LOV_twitch*2/(Lthick - Lbare); % Fraction of thick filament overlap

% Active force
B_process_twitch = kstiff2 * dr * p3_0_twitch;   % Force due to XB ratcheting
C_process_twitch = kstiff1 * (p2_1_twitch + p3_1_twitch );% Force due to stretching of XBs
F_XB_twitch = N_overlap_thick_twitch.*( B_process_twitch + C_process_twitch ); % ctive force

% (Linear) passive force model
sigma_collagen_twitch  = PConcollagen*(exp(PExpcollagen*(SL_twitch - SLcollagen)) - 1).*(SL_twitch > SLcollagen); % Collagen force
F_passive_twitch  = k_passive*(SL_twitch/2-Lsref/2)   + sigma_collagen_twitch; % Passive force
Ftotal_twitch = F_XB_twitch + F_passive_twitch; % Total force

T_final_XB = t_twitch*1000; % Convert to ms
force_final = Ftotal_twitch/10; % Store force output
Shortening_final = SL_twitch; % Store shortening output

%% Get just last beat
for i = 1:length(T_final_XB)
    if T_final_XB(i) > last_beat
        idx_XB = i;
        break
    end
end
end

%% Digitized experimental data (interpolating at specified points)
tspan = 0:0.001:2;
% Force pCa (Regnier et al. 2004)
Ca_fraction_fpca = [0.1:0.1:100];
ATP_points = csvread('ATP_points.csv');
dATP_points = csvread('dATP_points.csv');
ATP_Hill = csvread('ATP_Hill.csv');
ATP_Hill1 = ATP_Hill(:,1);
ATP_Hill2 = ATP_Hill(:,2);
dATP_Hill = csvread('dATP_Hill.csv');
dATP_Hill1 = dATP_Hill(:,1);
dATP_Hill2 = dATP_Hill(:,2);
[x, index] = unique(ATP_Hill1); 
ATP_Hill_interp = interp1(x, ATP_Hill2(index), Ca_fraction_fpca);
ATP_Hill_interp(1:40) = 0.9818;
ATP_Hill_interp(58:end) = 0.0061;
[x, index] = unique(dATP_Hill1); 
dATP_Hill_interp = interp1(x, dATP_Hill2(index), Ca_fraction_fpca);

% Ktr (Regnier et al. 2004)
ATP_data = csvread('pCa_4_ATP.csv');
ATP_data1 = ATP_data(:,1);
ATP_data2 = ATP_data(:,2);
[x, index] = unique(ATP_data1+0.5); 
ATP_pCa_4 = interp1(x, ATP_data2(index), tspan);
ATP_pCa_4(1:106) = 1.009;

dATP_data = csvread('pCa_4_dATP.csv');
dATP_data1 = dATP_data(:,1);
dATP_data2 = dATP_data(:,2);
[x, index] = unique(dATP_data1+0.5); 
dATP_pCa_4 = interp1(x, dATP_data2(index), tspan);
dATP_pCa_4(1:103) = 1.3114;

% Normalized Ktr
ATP_pCa_4_norm = ATP_pCa_4./ATP_pCa_4(1);
dATP_pCa_4_norm = dATP_pCa_4./ATP_pCa_4(1);

%% Plot
if suppress_plotting == 1
if XB_protocol == 0 || XB_protocol == 3 % Force pCa or all     
figure
hold on 
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca./38.5,'linewidth',3) % Normalize to max value
set (gca, 'xdir', 'reverse' )
xlabel('pCa');
ylabel('Relative Force');
plot(ATP_points(:,1),ATP_points(:,2),'o','markersize',8','linewidth',2,'color',[0.2, 0.1333, 0.5333]) % ATP data
plot(dATP_points(:,1),dATP_points(:,2),'o','markersize',8,'linewidth',2,'color',[0.2667, 0.6667, 0.6]) % dATP data
legend('Model','ATP data','dATP data')
xlim([4 7])
ylim([0 1.5])
set(gca,'FontSize',14)


else if XB_protocol == 1 || XB_protocol == 3 % Ktr or all 
tspan = 0:0.001:2;
figure
hold on
plot(t_ktr-4.0580,Ftotal_ktr./38.9814,'linewidth',3) % pCa 4 (normalized to dATP), shifted to start at 0
plot(tspan,ATP_pCa_4_norm,'-.','linewidth',3,'color',[0.2, 0.1333, 0.5333])
plot(tspan,dATP_pCa_4_norm,'-.','linewidth',3,'color',[0.2667, 0.6667, 0.6])
xlabel('Time (s)')
ylabel('Normalized force')
legend('Model','ATP data','dATP data')
xlim([0 2])
set(gca,'FontSize',14)

else if XB_protocol == 2 || XB_protocol == 3
% Force
figure
hold on
plot(T_final_XB(idx_XB:end)/1000-last_beat/1000,force_final(idx_XB:end),'linewidth',3)
xlabel('Time (s)')
ylabel('Force (kPa)')
set(gca,'FontSize',14)

% Shortening
figure
hold on
plot(T_final_XB(idx_XB:end)/1000-last_beat/1000,Shortening_final(idx_XB:end)./max(Shortening_final(idx_XB:end)),'linewidth',3)
xlabel('Time (s)')
ylabel('Relative Shortening')
set(gca,'FontSize',14)

% Sarcomere length
figure
hold on
plot(T_final_XB(idx_XB:end)/1000-last_beat/1000,Shortening_final(idx_XB:end),'linewidth',3)
xlabel('Time (s)')
ylabel('Sarcomere length')
set(gca,'FontSize',14)

end
end
end
end

%% Calculate
if XB_protocol == 0 || XB_protocol == 3 % Force pCa or all
[hill, ec50_n] = pCa_calculate((Ca_fraction_fpca')*10^(-6),(SS_Ftotal_fpca./max(SS_Ftotal_fpca))')
Ca50 = -log10(ec50_n)
peak_SS_force = max(SS_Ftotal_fpca)

else if XB_protocol == 1 || XB_protocol == 3 % Ktr or all
for i = 1:length(Ftotal_ktr)
    if t_ktr(i) > 4.5469
        idx_sl = i;
        break
    end
end

Force_Ktr_norm = Ftotal_ktr./Ftotal_ktr(idx_sl-1000);
Fnew = Force_Ktr_norm(idx_sl+2:end);
tnew = t_ktr(idx_sl+2:end);
[min_force,idx_min] = min(Fnew);
[max_force,idx_max] = max(Fnew);
force_t_1_2 = (max_force-min_force)*0.5 + min_force;
for i = idx_min:length(Fnew)
    if Fnew(i) > force_t_1_2 % max_force/2
        t_1_2a = tnew(i);
        t_1_2 = t_1_2a - t_ktr(idx_sl+2)
        break
    end
end
Ktr = 1/(1.264*t_1_2)
        
else if XB_protocol == 2 || XB_protocol == 3 % Twitch or all
force_final_last = force_final(idx_XB:end);
T_XB_final_last = T_final_XB(idx_XB:end);
Shortening_final_last = Shortening_final(idx_XB:end);

% Twitch
[min_force,idx_min] = min(force_final_last)
[max_force,idx_max] = max(force_final_last)
TTP_twitch = T_XB_final_last(idx_max) - last_beat
RT50 = (max(force_final_last)-min(force_final_last))*0.5 + min(force_final_last);
time_RT50_twitch = 0;
for i = idx_max:length(force_final_last)
    if force_final_last(i) < RT50
        time_RT50_twitch = ((T_XB_final_last(i)+T_XB_final_last(i-1))/2);
        time_RT50_twitch = time_RT50_twitch - last_beat - TTP_twitch
        break
    end
end

RT90 = (max(force_final_last)-min(force_final_last))*0.1 + min(force_final_last);
time_RT90_twitch = 0;
for i = idx_max:length(force_final_last)
    if force_final_last(i) < RT90
        time_RT90_twitch = ((T_XB_final_last(i)+T_XB_final_last(i-1))/2);
        time_RT90_twitch = time_RT90_twitch - last_beat - TTP_twitch
        break
    end
end

% Shortening
peak_shortening = min(Shortening_final_last./Shortening_final_last(1));
FS = -((min(Shortening_final_last)-max(Shortening_final_last))/max(Shortening_final_last))*100
[min_SL,idx_min] = min(Shortening_final_last);
TTP_shortening = T_XB_final_last(idx_min) - last_beat

RT50 = max(Shortening_final_last) - (max(Shortening_final_last)-min(Shortening_final_last))*0.5;
time_RT50_shortening = 0;
for i = idx_min:length(Shortening_final_last)
    if Shortening_final_last(i) > RT50
        time_RT50_shortening = ((T_XB_final_last(i)+T_XB_final_last(i-1))/2);
        time_RT50_shortening = time_RT50_shortening - last_beat - TTP_shortening
        break
    end
end

RT90 = max(Shortening_final_last) - (max(Shortening_final_last)-min(Shortening_final_last))*0.1;
time_RT90_shortening = 0;
for i = idx_min:length(Shortening_final_last)
    if Shortening_final_last(i) > RT90
        time_RT90_shortening = ((T_XB_final_last(i)+T_XB_final_last(i-1))/2);
        time_RT90_shortening = time_RT90_shortening - last_beat - TTP_shortening
        break
    end
end

end
end
end

toc

% For output
if XB_protocol == 0 % force pCa
Ftotal_ktr = 0;
t_ktr = 0;
T_final_XB = 0;
idx_XB = 0;
Shortening_final = 0;
force_final = 0;  

else if XB_protocol == 1 % Ktr
SS_Ftotal_fpca = 0;
T_final_XB = 0;
idx_XB = 0;
Shortening_final = 0;
force_final = 0; 

else if XB_protocol == 2 % Twitch    
SS_Ftotal_fpca = 0;
Ftotal_ktr = 0;
t_ktr = 0;
end
end
end


end