%% Digitized experimental data (interpolating at specified points)
% ***Run this section before any of other figure sections
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

%% Figure 3
% A-D
[~, ~, ~, ~, SS_Ftotal_fpca_ATP, ~, ~] = myocyte_model(0, 0, 0, 0, 0, 250, 304.7, 4, 2, 80, 4, 25, 0.4, 85, 900, 4); % ATP force pCa
[~, ~, ~, ~, SS_Ftotal_fpca_dATP_ka, ~, ~] = myocyte_model(0, 0, 0, 0, 100, 538, 304.7, 4, 2, 80, 4, 25, 0.4, 85, 900, 4); % dATP ka force pCa
[~, ~, ~, ~, SS_Ftotal_fpca_dATP_ka_kd_k1, ~, ~] = myocyte_model(0, 0, 0, 0, 100, 538, 540, 5.7, 2, 80, 4, 25, 0.4, 85, 900, 4); % dATP XB (ka, kd, k1) force pCa
[~, ~, ~, ~, ~, Ftotal_ktr_ATP, t_ktr_ATP] = myocyte_model(0, 1, 0, 0, 0, 250, 304.7, 4, 2, 80, 4, 25, 0.4, 85, 900, 4); % ATP Ktr
[~, ~, ~, ~, ~, Ftotal_ktr_dATP_ka, t_ktr_dATP_ka] = myocyte_model(0, 1, 0, 0, 100, 538, 304.7, 4, 2, 80, 4, 25, 0.4, 85, 900, 4); % dATP ka Ktr
[~, ~, ~, ~, ~, Ftotal_ktr_dATP_ka_kd_k1, t_ktr_dATP_ka_kd_k1] = myocyte_model(0, 1, 0, 0, 100, 538, 540, 5.7, 2, 80, 4, 25, 0.4, 85, 900, 4); % dATP XB (ka, kd, k1) Ktr

figure
subplot(2,2,1)
hold on
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_ATP/max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.2, 0.1333, 0.5333])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_dATP_ka./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.2667, 0.6667, 0.6])
plot(ATP_points(:,1),ATP_points(:,2),'o','markersize',8','linewidth',2,'color',[0.2, 0.1333, 0.5333]) 
plot(dATP_points(:,1),dATP_points(:,2),'o','markersize',8,'linewidth',2,'color',[0.2667, 0.6667, 0.6]) 
xlim([4 7])
ylim([0 1.5])
set(gca,'FontSize',14)
set (gca, 'xdir', 'reverse' )
xlabel('pCa');
ylabel('Relative force');
legend('ATP model','dATP model','ATP data','dATP data')

subplot(2,2,2)
hold on
tspan = 0:0.001:2;
plot(t_ktr_ATP-4.0580,Ftotal_ktr_ATP./max(Ftotal_ktr_ATP),'linewidth',3,'color',[0.2, 0.1333, 0.5333])
plot(t_ktr_dATP_ka-4.0580,Ftotal_ktr_dATP_ka./max(Ftotal_ktr_ATP),'linewidth',3,'color',[0.2667, 0.6667, 0.6]) 
plot(tspan,ATP_pCa_4_norm,'-.','linewidth',3,'color',[0.2, 0.1333, 0.5333])
plot(tspan,dATP_pCa_4_norm,'-.','linewidth',3,'color',[0.2667, 0.6667, 0.6])
xlabel('Time (s)')
ylabel('Relative force')
legend('ATP model','dATP model','ATP data','dATP data')
xlim([0 2])
ylim([0 1.5])
set(gca,'FontSize',14)

subplot(2,2,3)
hold on
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_ATP/max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.2, 0.1333, 0.5333])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_dATP_ka_kd_k1./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.2667, 0.6667, 0.6])
plot(ATP_points(:,1),ATP_points(:,2),'o','markersize',8','linewidth',2,'color',[0.2, 0.1333, 0.5333]) 
plot(dATP_points(:,1),dATP_points(:,2),'o','markersize',8,'linewidth',2,'color',[0.2667, 0.6667, 0.6]) 
xlim([4 7])
ylim([0 1.5])
set(gca,'FontSize',14)
set (gca, 'xdir', 'reverse' )
xlabel('pCa');
ylabel('Relative force');
legend('ATP model','dATP model','ATP data','dATP data')

subplot(2,2,4)
hold on
tspan = 0:0.001:2;
plot(t_ktr_ATP-4.0580,Ftotal_ktr_ATP./max(Ftotal_ktr_ATP),'linewidth',3,'color',[0.2, 0.1333, 0.5333]) 
plot(t_ktr_dATP_ka_kd_k1-4.0580,Ftotal_ktr_dATP_ka_kd_k1./max(Ftotal_ktr_ATP),'linewidth',3,'color',[0.2667, 0.6667, 0.6]) 
plot(tspan,ATP_pCa_4_norm,'-.','linewidth',3,'color',[0.2, 0.1333, 0.5333])
plot(tspan,dATP_pCa_4_norm,'-.','linewidth',3,'color',[0.2667, 0.6667, 0.6])
xlabel('Time (s)')
ylabel('Relative force')
legend('ATP model','dATP model','ATP data','dATP data')
xlim([0 2])
ylim([0 1.5])
set(gca,'FontSize',14)

%E
[T_final_XB_ATP, force_final_ATP, idx_XB_ATP, Shortening_final_ATP, ~, ~, ~] = myocyte_model(0, 2, 0, 0, 0, 250, 304.7, 4, 2, 80, 4, 25, 0.4, 85, 900, 4); % ATP shortening
[T_final_XB_dATP_XB_SRX_Ca, force_final_dATP_XB_SRX_Ca, idx_XB_dATP_XB_SRX_Ca, Shortening_final_dATP_XB_SRX_Ca, ~, ~, ~] = myocyte_model(1, 2, 0, 0, 2, 538, 540, 5.7, 2, 80, 4, 25, 350, 85, 900, 4); % dATP XB + SRX + Ca shortening

last_beat = 2000;

figure
hold on
plot(T_final_XB_ATP(idx_XB_ATP:end)/1000-last_beat/1000,Shortening_final_ATP(idx_XB_ATP:end)./max(Shortening_final_ATP(idx_XB_ATP:end)),'linewidth',3,'color',[0.2, 0.1333, 0.5333])
plot(T_final_XB_dATP_XB_SRX_Ca(idx_XB_dATP_XB_SRX_Ca:end)/1000-last_beat/1000,Shortening_final_dATP_XB_SRX_Ca(idx_XB_dATP_XB_SRX_Ca:end)./max(Shortening_final_dATP_XB_SRX_Ca(idx_XB_dATP_XB_SRX_Ca:end)),'linewidth',3,'color',[0.2667, 0.6667, 0.6])
xlabel('Time (s)')
ylabel('Relative shortening')
legend('ATP model','dATP model')
set(gca,'FontSize',14)

%F
[T_final_XB_ATP, force_final_ATP, idx_XB_ATP, Shortening_final_ATP, ~, ~, ~] = myocyte_model(0, 2, 0, 0, 0, 250, 304.7, 4, 2, 80, 4, 25, 0.4, 85, 900, 4); % ATP shortening
force_final_last = force_final_ATP(idx_XB_ATP:end);
T_XB_final_last = T_final_XB_ATP(idx_XB_ATP:end);
Shortening_final_last = Shortening_final_ATP(idx_XB_ATP:end);
peak_shortening = min(Shortening_final_last./Shortening_final_last(1));
FS = -((min(Shortening_final_last)-max(Shortening_final_last))/max(Shortening_final_last))*100;
[min_SL,idx_min] = min(Shortening_final_last);
TTP_shortening = T_XB_final_last(idx_min) - last_beat;

RT50 = max(Shortening_final_last) - (max(Shortening_final_last)-min(Shortening_final_last))*0.5;
time_RT50_shortening = 0;
for i = idx_min:length(Shortening_final_last)
    if Shortening_final_last(i) > RT50
        time_RT50_shortening = ((T_XB_final_last(i)+T_XB_final_last(i-1))/2);
        time_RT50_shortening = time_RT50_shortening - last_beat - TTP_shortening;
        break
    end
end

RT90 = max(Shortening_final_last) - (max(Shortening_final_last)-min(Shortening_final_last))*0.1;
time_RT90_shortening = 0;
for i = idx_min:length(Shortening_final_last)
    if Shortening_final_last(i) > RT90
        time_RT90_shortening = ((T_XB_final_last(i)+T_XB_final_last(i-1))/2);
        time_RT90_shortening = time_RT90_shortening - last_beat - TTP_shortening;
        break
    end
end

FS_store(1) = FS;
TTP_store(1) = TTP_shortening;
RT50_store(1) = time_RT50_shortening;
RT90_store(1) = time_RT90_shortening;

[T_final_XB_dATP_XB, force_final_dATP_XB, idx_XB_dATP_XB, Shortening_final_dATP_XB, ~, ~, ~] = myocyte_model(0, 2, 0, 0, 2, 538, 540, 5.7, 2, 80, 4, 25, 0.4, 85, 900, 4); % dATP XB shortening
force_final_last = force_final_dATP_XB(idx_XB_dATP_XB:end);
T_XB_final_last = T_final_XB_dATP_XB(idx_XB_dATP_XB:end);
Shortening_final_last = Shortening_final_dATP_XB(idx_XB_dATP_XB:end);
peak_shortening = min(Shortening_final_last./Shortening_final_last(1));
FS = -((min(Shortening_final_last)-max(Shortening_final_last))/max(Shortening_final_last))*100;
[min_SL,idx_min] = min(Shortening_final_last);
TTP_shortening = T_XB_final_last(idx_min) - last_beat;

RT50 = max(Shortening_final_last) - (max(Shortening_final_last)-min(Shortening_final_last))*0.5;
time_RT50_shortening = 0;
for i = idx_min:length(Shortening_final_last)
    if Shortening_final_last(i) > RT50
        time_RT50_shortening = ((T_XB_final_last(i)+T_XB_final_last(i-1))/2);
        time_RT50_shortening = time_RT50_shortening - last_beat - TTP_shortening;
        break
    end
end

RT90 = max(Shortening_final_last) - (max(Shortening_final_last)-min(Shortening_final_last))*0.1;
time_RT90_shortening = 0;
for i = idx_min:length(Shortening_final_last)
    if Shortening_final_last(i) > RT90
        time_RT90_shortening = ((T_XB_final_last(i)+T_XB_final_last(i-1))/2);
        time_RT90_shortening = time_RT90_shortening - last_beat - TTP_shortening;
        break
    end
end

FS_store(2) = FS;
TTP_store(2) = TTP_shortening;
RT50_store(2) = time_RT50_shortening;
RT90_store(2) = time_RT90_shortening;

[T_final_XB_dATP_SRX, force_final_dATP_SRX, idx_XB_dATP_SRX, Shortening_final_dATP_SRX, ~, ~, ~] = myocyte_model(0, 2, 0, 0, 2, 250, 304.7, 4, 2, 80, 4, 25, 350, 85, 900, 4); % dATP SRX shortening
force_final_last = force_final_dATP_SRX(idx_XB_dATP_SRX:end);
T_XB_final_last = T_final_XB_dATP_SRX(idx_XB_dATP_SRX:end);
Shortening_final_last = Shortening_final_dATP_SRX(idx_XB_dATP_SRX:end);
peak_shortening = min(Shortening_final_last./Shortening_final_last(1));
FS = -((min(Shortening_final_last)-max(Shortening_final_last))/max(Shortening_final_last))*100;
[min_SL,idx_min] = min(Shortening_final_last);
TTP_shortening = T_XB_final_last(idx_min) - last_beat;

RT50 = max(Shortening_final_last) - (max(Shortening_final_last)-min(Shortening_final_last))*0.5;
time_RT50_shortening = 0;
for i = idx_min:length(Shortening_final_last)
    if Shortening_final_last(i) > RT50
        time_RT50_shortening = ((T_XB_final_last(i)+T_XB_final_last(i-1))/2);
        time_RT50_shortening = time_RT50_shortening - last_beat - TTP_shortening;
        break
    end
end

RT90 = max(Shortening_final_last) - (max(Shortening_final_last)-min(Shortening_final_last))*0.1;
time_RT90_shortening = 0;
for i = idx_min:length(Shortening_final_last)
    if Shortening_final_last(i) > RT90
        time_RT90_shortening = ((T_XB_final_last(i)+T_XB_final_last(i-1))/2);
        time_RT90_shortening = time_RT90_shortening - last_beat - TTP_shortening;
        break
    end
end

FS_store(3) = FS;
TTP_store(3) = TTP_shortening;
RT50_store(3) = time_RT50_shortening;
RT90_store(3) = time_RT90_shortening;

[T_final_XB_dATP_Ca, force_final_dATP_Ca, idx_XB_dATP_Ca, Shortening_final_dATP_Ca, ~, ~, ~] = myocyte_model(1, 2, 0, 0, 2, 250, 304.7, 4, 2, 80, 4, 25, 0.4, 85, 900, 4); % dATP Ca shortening
force_final_last = force_final_dATP_Ca(idx_XB_dATP_Ca:end);
T_XB_final_last = T_final_XB_dATP_Ca(idx_XB_dATP_Ca:end);
Shortening_final_last = Shortening_final_dATP_Ca(idx_XB_dATP_Ca:end);
peak_shortening = min(Shortening_final_last./Shortening_final_last(1));
FS = -((min(Shortening_final_last)-max(Shortening_final_last))/max(Shortening_final_last))*100;
[min_SL,idx_min] = min(Shortening_final_last);
TTP_shortening = T_XB_final_last(idx_min) - last_beat;

RT50 = max(Shortening_final_last) - (max(Shortening_final_last)-min(Shortening_final_last))*0.5;
time_RT50_shortening = 0;
for i = idx_min:length(Shortening_final_last)
    if Shortening_final_last(i) > RT50
        time_RT50_shortening = ((T_XB_final_last(i)+T_XB_final_last(i-1))/2);
        time_RT50_shortening = time_RT50_shortening - last_beat - TTP_shortening;
        break
    end
end

RT90 = max(Shortening_final_last) - (max(Shortening_final_last)-min(Shortening_final_last))*0.1;
time_RT90_shortening = 0;
for i = idx_min:length(Shortening_final_last)
    if Shortening_final_last(i) > RT90
        time_RT90_shortening = ((T_XB_final_last(i)+T_XB_final_last(i-1))/2);
        time_RT90_shortening = time_RT90_shortening - last_beat - TTP_shortening;
        break
    end
end

FS_store(4) = FS;
TTP_store(4) = TTP_shortening;
RT50_store(4) = time_RT50_shortening;
RT90_store(4) = time_RT90_shortening;

[T_final_XB_dATP_XB_SRX_Ca, force_final_dATP_XB_SRX_Ca, idx_XB_dATP_XB_SRX_Ca, Shortening_final_dATP_XB_SRX_Ca, ~, ~, ~] = myocyte_model(1, 2, 0, 0, 2, 538, 540, 5.7, 2, 80, 4, 25, 350, 85, 900, 4); % dATP XB + SRX + Ca shortening
force_final_last = force_final_dATP_XB_SRX_Ca(idx_XB_dATP_XB_SRX_Ca:end);
T_XB_final_last = T_final_XB_dATP_XB_SRX_Ca(idx_XB_dATP_XB_SRX_Ca:end);
Shortening_final_last = Shortening_final_dATP_XB_SRX_Ca(idx_XB_dATP_XB_SRX_Ca:end);
peak_shortening = min(Shortening_final_last./Shortening_final_last(1));
FS = -((min(Shortening_final_last)-max(Shortening_final_last))/max(Shortening_final_last))*100;
[min_SL,idx_min] = min(Shortening_final_last);
TTP_shortening = T_XB_final_last(idx_min) - last_beat;

RT50 = max(Shortening_final_last) - (max(Shortening_final_last)-min(Shortening_final_last))*0.5;
time_RT50_shortening = 0;
for i = idx_min:length(Shortening_final_last)
    if Shortening_final_last(i) > RT50
        time_RT50_shortening = ((T_XB_final_last(i)+T_XB_final_last(i-1))/2);
        time_RT50_shortening = time_RT50_shortening - last_beat - TTP_shortening;
        break
    end
end

RT90 = max(Shortening_final_last) - (max(Shortening_final_last)-min(Shortening_final_last))*0.1;
time_RT90_shortening = 0;
for i = idx_min:length(Shortening_final_last)
    if Shortening_final_last(i) > RT90
        time_RT90_shortening = ((T_XB_final_last(i)+T_XB_final_last(i-1))/2);
        time_RT90_shortening = time_RT90_shortening - last_beat - TTP_shortening;
        break
    end
end

FS_store(5) = FS;
TTP_store(5) = TTP_shortening;
RT50_store(5) = time_RT50_shortening;
RT90_store(5) = time_RT90_shortening;

% From experimental data (Korte et al. 2011) (% change)
FS_store(6) = 0.55;
TTP_store(6) = -0.03;
RT50_store(6) = -0.29;
RT90_store(6) = -0.43;

x_axis = [0 1 2 3];

figure
subplot(4,1,1)
hold on
plot(x_axis,(FS_store(2:5)-FS_store(1))./FS_store(1),'x','color',[0.5, 0.5, 0.5],'linewidth',3,'markersize',12)
yline((FS_store(1)-FS_store(1))./FS_store(1),':','color',[0.2, 0.1333, 0.5333],'linewidth',3);
yline(FS_store(6),':','color',[0.2667, 0.6667, 0.6],'linewidth',3);
ylabel('FS')
xticks(x_axis)
set(gca,'FontSize',14,'XTick',[],'XTickLabel',{'','','',''})

subplot(4,1,2)
hold on
plot(x_axis,(TTP_store(2:5)-TTP_store(1))./TTP_store(1),'x','color',[0.5, 0.5, 0.5],'linewidth',3,'markersize',12)
yline((TTP_store(1)-TTP_store(1))./TTP_store(1),':','color',[0.2, 0.1333, 0.5333],'linewidth',3);
yline(TTP_store(6),':','color',[0.2667, 0.6667, 0.6],'linewidth',3);
ylabel('TTP')
xticks(x_axis)
set(gca,'FontSize',14,'XTick',[],'XTickLabel',{'','','',''})

subplot(4,1,3)
hold on
plot(x_axis,(RT50_store(2:5)-RT50_store(1))./RT50_store(1),'x','color',[0.5, 0.5, 0.5],'linewidth',3,'markersize',12)
yline((RT50_store(1)-RT50_store(1))./RT50_store(1),':','color',[0.2, 0.1333, 0.5333],'linewidth',3);
yline(RT50_store(6),':','color',[0.2667, 0.6667, 0.6],'linewidth',3);
ylabel('RT50')
xticks(x_axis)
set(gca,'FontSize',14,'XTick',[],'XTickLabel',{'','','',''})

subplot(4,1,4)
hold on
plot(x_axis,(RT90_store(2:5)-RT90_store(1))./RT90_store(1),'x','color',[0.5, 0.5, 0.5],'linewidth',3,'markersize',12)
yline((RT90_store(1)-RT90_store(1))./RT90_store(1),':','color',[0.2, 0.1333, 0.5333],'linewidth',3);
yline(RT90_store(6),':','color',[0.2667, 0.6667, 0.6],'linewidth',3);
ylabel('RT90')
xticks(x_axis)
set(gca,'FontSize',14,'XTickLabel',...
    {'XB','SRX','Ca','XB + Ca + SRX'})

%% Figure S3
[~, ~, ~, ~, SS_Ftotal_fpca_ATP, ~, ~] = myocyte_model(0, 0, 0, 0, 0, 250, 304.7, 4, 2, 80, 4, 25, 0.4, 85, 900, 4); % ATP force pCa
[~, ~, ~, ~, SS_Ftotal_fpca_dATP, ~, ~] = myocyte_model(0, 0, 0, 0, 100, 538, 540, 5.7, 2, 80, 4, 25, 0.4, 85, 900, 4); % dATP force pCa
[~, ~, ~, ~, SS_Ftotal_fpca_ATP_shifted, ~, ~] = myocyte_model(0, 0, 0, 0, 100, 250, 304.7, 4, 2, 80, 4, 25, 0.4, 280, 800, 4); % ATP force pCa shifted
[~, ~, ~, ~, SS_Ftotal_fpca_dATP_shifted, ~, ~] = myocyte_model(0, 0, 0, 0, 100, 538, 540, 5.7, 2, 80, 4, 25, 0.4, 280, 800, 4); % dATP force pCa shifted

figure
hold on
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_ATP/max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.2, 0.1333, 0.5333])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_dATP./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.2667, 0.6667, 0.6])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_ATP_shifted/max(SS_Ftotal_fpca_ATP_shifted),'--','linewidth',3,'color',[0.2, 0.1333, 0.5333])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_dATP_shifted./max(SS_Ftotal_fpca_ATP_shifted),'--','linewidth',3,'color',[0.2667, 0.6667, 0.6])
plot(ATP_points(:,1),ATP_points(:,2),'o','markersize',8','linewidth',2,'color',[0.2, 0.1333, 0.5333]) 
plot(dATP_points(:,1),dATP_points(:,2),'o','markersize',8,'linewidth',2,'color',[0.2667, 0.6667, 0.6]) 
xlim([4 7])
ylim([0 1.6])
set(gca,'FontSize',14)
set (gca, 'xdir', 'reverse' )
xlabel('pCa');
ylabel('Relative Force');
legend('ATP model','dATP model','ATP data','dATP data')

%% Figure S4
% ATP
[~, ~, ~, ~, SS_Ftotal_fpca_ATP, ~, ~] = myocyte_model(0, 0, 0, 0, 0, 250, 304.7, 4, 2, 80, 4, 25, 0.4, 85, 900, 4); % ATP force pCa

% ka
[~, ~, ~, ~, SS_Ftotal_fpca_ka_10, ~, ~] = myocyte_model(0, 0, 0, 0, 100, 250*0.1, 304.7, 4, 2, 80, 4, 25, 0.4, 85, 900, 4); % force pCa scaling ka 10%
[~, ~, ~, ~, SS_Ftotal_fpca_ka_50, ~, ~] = myocyte_model(0, 0, 0, 0, 100, 250*0.5, 304.7, 4, 2, 80, 4, 25, 0.4, 85, 900, 4); % force pCa scaling ka 50%
[~, ~, ~, ~, SS_Ftotal_fpca_ka_150, ~, ~] = myocyte_model(0, 0, 0, 0, 100, 250*1.5, 304.7, 4, 2, 80, 4, 25, 0.4, 85, 900, 4); % force pCa scaling ka 150%
[~, ~, ~, ~, SS_Ftotal_fpca_ka_200, ~, ~] = myocyte_model(0, 0, 0, 0, 100, 250*2, 304.7, 4, 2, 80, 4, 25, 0.4, 85, 900, 4); % force pCa scaling ka 200%

% kd
[~, ~, ~, ~, SS_Ftotal_fpca_kd_10, ~, ~] = myocyte_model(0, 0, 0, 0, 100, 250, 304.7*0.1, 4, 2, 80, 4, 25, 0.4, 85, 900, 4); % force pCa scaling kd 10%
[~, ~, ~, ~, SS_Ftotal_fpca_kd_50, ~, ~] = myocyte_model(0, 0, 0, 0, 100, 250, 304.7*0.5, 4, 2, 80, 4, 25, 0.4, 85, 900, 4); % force pCa scaling kd 50%
[~, ~, ~, ~, SS_Ftotal_fpca_kd_150, ~, ~] = myocyte_model(0, 0, 0, 0, 100, 250, 304.7*1.5, 4, 2, 80, 4, 25, 0.4, 85, 900, 4); % force pCa scaling kd 150%
[~, ~, ~, ~, SS_Ftotal_fpca_kd_200, ~, ~] = myocyte_model(0, 0, 0, 0, 100, 250, 304.7*2, 4, 2, 80, 4, 25, 0.4, 85, 900, 4); % force pCa scaling kd 200%

% k1
[~, ~, ~, ~, SS_Ftotal_fpca_k1_10, ~, ~] = myocyte_model(0, 0, 0, 0, 100, 250, 304.7, 4*0.1, 2, 80, 4, 25, 0.4, 85, 900, 4); % force pCa scaling k1 10%
[~, ~, ~, ~, SS_Ftotal_fpca_k1_50, ~, ~] = myocyte_model(0, 0, 0, 0, 100, 250, 304.7, 4*0.5, 2, 80, 4, 25, 0.4, 85, 900, 4); % force pCa scaling k1 50%
[~, ~, ~, ~, SS_Ftotal_fpca_k1_150, ~, ~] = myocyte_model(0, 0, 0, 0, 100, 250, 304.7, 4*1.5, 2, 80, 4, 25, 0.4, 85, 900, 4); % force pCa scaling k1 150%
[~, ~, ~, ~, SS_Ftotal_fpca_k1_200, ~, ~] = myocyte_model(0, 0, 0, 0, 100, 250, 304.7, 4*2, 2, 80, 4, 25, 0.4, 85, 900, 4); % force pCa scaling k1 200%

% k_1
[~, ~, ~, ~, SS_Ftotal_fpca_k_1_10, ~, ~] = myocyte_model(0, 0, 0, 0, 100, 250, 304.7, 4, 2*0.1, 80, 4, 25, 0.4, 85, 900, 4); % force pCa scaling k_1 10%
[~, ~, ~, ~, SS_Ftotal_fpca_k_1_50, ~, ~] = myocyte_model(0, 0, 0, 0, 100, 250, 304.7, 4, 2*0.5, 80, 4, 25, 0.4, 85, 900, 4); % force pCa scaling k_1 50%
[~, ~, ~, ~, SS_Ftotal_fpca_k_1_150, ~, ~] = myocyte_model(0, 0, 0, 0, 100, 250, 304.7, 4, 2*1.5, 80, 4, 25, 0.4, 85, 900, 4); % force pCa scaling k_1 150%
[~, ~, ~, ~, SS_Ftotal_fpca_k_1_200, ~, ~] = myocyte_model(0, 0, 0, 0, 100, 250, 304.7, 4, 2*2, 80, 4, 25, 0.4, 85, 900, 4); % force pCa scaling k_1 200%

% k2
[~, ~, ~, ~, SS_Ftotal_fpca_k2_10, ~, ~] = myocyte_model(0, 0, 0, 0, 100, 250, 304.7, 4, 2, 80*0.1, 4, 25, 0.4, 85, 900, 4); % force pCa scaling k2 10%
[~, ~, ~, ~, SS_Ftotal_fpca_k2_50, ~, ~] = myocyte_model(0, 0, 0, 0, 100, 250, 304.7, 4, 2, 80*0.5, 4, 25, 0.4, 85, 900, 4); % force pCa scaling k2 50%
[~, ~, ~, ~, SS_Ftotal_fpca_k2_150, ~, ~] = myocyte_model(0, 0, 0, 0, 100, 250, 304.7, 4, 2, 80*1.5, 4, 25, 0.4, 85, 900, 4); % force pCa scaling k2 150%
[~, ~, ~, ~, SS_Ftotal_fpca_k2_200, ~, ~] = myocyte_model(0, 0, 0, 0, 100, 250, 304.7, 4, 2, 80*2, 4, 25, 0.4, 85, 900, 4); % force pCa scaling k2 200%

% k_2
[~, ~, ~, ~, SS_Ftotal_fpca_k_2_10, ~, ~] = myocyte_model(0, 0, 0, 0, 100, 250, 304.7, 4, 2, 80, 4*0.1, 25, 0.4, 85, 900, 4); % force pCa scaling k_2 10%
[~, ~, ~, ~, SS_Ftotal_fpca_k_2_50, ~, ~] = myocyte_model(0, 0, 0, 0, 100, 250, 304.7, 4, 2, 80, 4*0.5, 25, 0.4, 85, 900, 4); % force pCa scaling k_2 50%
[~, ~, ~, ~, SS_Ftotal_fpca_k_2_150, ~, ~] = myocyte_model(0, 0, 0, 0, 100, 250, 304.7, 4, 2, 80, 4*1.5, 25, 0.4, 85, 900, 4); % force pCa scaling k_2 150%
[~, ~, ~, ~, SS_Ftotal_fpca_k_2_200, ~, ~] = myocyte_model(0, 0, 0, 0, 100, 250, 304.7, 4, 2, 80, 4*2, 25, 0.4, 85, 900, 4); % force pCa scaling k_2 200%

% k3
[~, ~, ~, ~, SS_Ftotal_fpca_k3_10, ~, ~] = myocyte_model(0, 0, 0, 0, 100, 250, 304.7, 4, 2, 80, 4, 25*0.1, 0.4, 85, 900, 4); % force pCa scaling k3 10%
[~, ~, ~, ~, SS_Ftotal_fpca_k3_50, ~, ~] = myocyte_model(0, 0, 0, 0, 100, 250, 304.7, 4, 2, 80, 4, 25*0.5, 0.4, 85, 900, 4); % force pCa scaling k3 50%
[~, ~, ~, ~, SS_Ftotal_fpca_k3_150, ~, ~] = myocyte_model(0, 0, 0, 0, 100, 250, 304.7, 4, 2, 80, 4, 25*1.5, 0.4, 85, 900, 4); % force pCa scaling k3 150%
[~, ~, ~, ~, SS_Ftotal_fpca_k3_200, ~, ~] = myocyte_model(0, 0, 0, 0, 100, 250, 304.7, 4, 2, 80, 4, 25*2, 0.4, 85, 900, 4); % force pCa scaling k3 200%


figure
subplot(2,4,1)
hold on
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_ATP/max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0, 0, 0])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_ka_10./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.20,0.13,0.53])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_ka_50./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.05,0.42,0.22])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_ka_150./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.87,0.80,0.47])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_ka_200./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.67,0.27,0.60])
xlim([4 7])
ylim([0 2])
set(gca,'FontSize',14)
set (gca, 'xdir', 'reverse' )
xlabel('pCa');
ylabel('Relative Force');
legend('ATP','10%','50%','150%','200%')
title('k_a')

subplot(2,4,2)
hold on
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_ATP/max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0, 0, 0])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_kd_10./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.20,0.13,0.53])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_kd_50./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.05,0.42,0.22])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_kd_150./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.87,0.80,0.47])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_kd_200./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.67,0.27,0.60])
xlim([4 7])
ylim([0 2])
set(gca,'FontSize',14)
set (gca, 'xdir', 'reverse' )
xlabel('pCa');
ylabel('Relative Force');
legend('ATP','10%','50%','150%','200%')
title('k_d')

subplot(2,4,3)
hold on
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_ATP/max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0, 0, 0])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_k1_10./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.20,0.13,0.53])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_k1_50./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.05,0.42,0.22])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_k1_150./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.87,0.80,0.47])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_k1_200./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.67,0.27,0.60])
xlim([4 7])
ylim([0 2])
set(gca,'FontSize',14)
set (gca, 'xdir', 'reverse' )
xlabel('pCa');
ylabel('Relative Force');
legend('ATP','10%','50%','150%','200%')
title('k_1')

subplot(2,4,4)
hold on
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_ATP/max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0, 0, 0])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_k_1_10./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.20,0.13,0.53])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_k_1_50./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.05,0.42,0.22])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_k_1_150./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.87,0.80,0.47])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_k_1_200./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.67,0.27,0.60])
xlim([4 7])
ylim([0 2])
set(gca,'FontSize',14)
set (gca, 'xdir', 'reverse' )
xlabel('pCa');
ylabel('Relative Force');
legend('ATP','10%','50%','150%','200%')
title('k_{-1}')

subplot(2,4,5)
hold on
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_ATP/max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0, 0, 0])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_k2_10./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.20,0.13,0.53])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_k2_50./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.05,0.42,0.22])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_k2_150./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.87,0.80,0.47])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_k2_200./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.67,0.27,0.60])
xlim([4 7])
ylim([0 2])
set(gca,'FontSize',14)
set (gca, 'xdir', 'reverse' )
xlabel('pCa');
ylabel('Relative Force');
legend('ATP','10%','50%','150%','200%')
title('k_2')

subplot(2,4,6)
hold on
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_ATP/max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0, 0, 0])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_k_2_10./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.20,0.13,0.53])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_k_2_50./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.05,0.42,0.22])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_k_2_150./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.87,0.80,0.47])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_k_2_200./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.67,0.27,0.60])
xlim([4 7])
ylim([0 2])
set(gca,'FontSize',14)
set (gca, 'xdir', 'reverse' )
xlabel('pCa');
ylabel('Relative Force');
legend('ATP','10%','50%','150%','200%')
title('k_{-2}')

subplot(2,4,7)
hold on
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_ATP/max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0, 0, 0])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_k3_10./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.20,0.13,0.53])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_k3_50./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.05,0.42,0.22])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_k3_150./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.87,0.80,0.47])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_k3_200./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.67,0.27,0.60])
xlim([4 7])
ylim([0 2])
set(gca,'FontSize',14)
set (gca, 'xdir', 'reverse' )
xlabel('pCa');
ylabel('Relative Force');
legend('ATP','10%','50%','150%','200%')
title('k_3')


%% Figure S5
% ATP
[~, ~, ~, ~, SS_Ftotal_fpca_ATP, ~, ~] = myocyte_model(0, 0, 0, 0, 0, 250, 304.7, 4, 2, 80, 4, 25, 0.4, 85, 900, 4); % ATP force pCa

% ka
[~, ~, ~, ~, SS_Ftotal_fpca_ka, ~, ~] = myocyte_model(0, 0, 0, 0, 100, 538, 304.7, 4, 2, 80, 4, 25, 0.4, 85, 900, 4); % force pCa scaling ka 10%

% ka + kd
[~, ~, ~, ~, SS_Ftotal_fpca_ka_kd, ~, ~] = myocyte_model(0, 0, 0, 0, 100, 538, 367, 4, 2, 80, 4, 25, 0.4, 85, 900, 4); % force pCa scaling ka 10%

% ka + k1
[~, ~, ~, ~, SS_Ftotal_fpca_ka_k1, ~, ~] = myocyte_model(0, 0, 0, 0, 100, 538, 304.7, 5, 2, 80, 4, 25, 0.4, 85, 900, 4); % force pCa scaling ka 10%

% ka + k3
[~, ~, ~, ~, SS_Ftotal_fpca_ka_k3, ~, ~] = myocyte_model(0, 0, 0, 0, 100, 538, 304.7, 4, 2, 80, 4, 30, 0.4, 85, 900, 4); % force pCa scaling ka 10%

% ka + kd + k1
[~, ~, ~, ~, SS_Ftotal_fpca_ka_kd_k1, ~, ~] = myocyte_model(0, 0, 0, 0, 100, 538, 540, 5.7, 2, 80, 4, 25, 0.4, 85, 900, 4); % force pCa scaling ka 10%

% ka + k1 + k3
[~, ~, ~, ~, SS_Ftotal_fpca_ka_k1_k3, ~, ~] = myocyte_model(0, 0, 0, 0, 100, 538, 304.7, 6, 2, 80, 4, 32, 0.4, 85, 900, 4); % force pCa scaling ka 10%

% ka + kd + k3
[~, ~, ~, ~, SS_Ftotal_fpca_ka_kd_k3, ~, ~] = myocyte_model(0, 0, 0, 0, 100, 538, 360, 4, 2, 80, 4, 26, 0.4, 85, 900, 4); % force pCa scaling ka 10%

% ka + kd + k1 + k3
[~, ~, ~, ~, SS_Ftotal_fpca_ka_kd_k1_k3, ~, ~] = myocyte_model(0, 0, 0, 0, 100, 538, 540, 5.7, 2, 80, 4, 25.5, 0.4, 85, 900, 4); % force pCa scaling ka 10%

figure
hold on
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_ATP/max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.2, 0.1333, 0.5333])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_ka/max(SS_Ftotal_fpca_ATP),':','linewidth',3,'color',[0.5, 0.5, 0.5])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_ka_kd/max(SS_Ftotal_fpca_ATP),'-.','linewidth',3,'color',[0.5, 0.5, 0.5])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_ka_k1/max(SS_Ftotal_fpca_ATP),'--','linewidth',3,'color',[0.5, 0.5, 0.5])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_ka_k3/max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.5, 0.5, 0.5])
plot(ATP_points(:,1),ATP_points(:,2),'o','markersize',8','linewidth',2,'color',[0.2, 0.1333, 0.5333]) 
plot(dATP_points(:,1),dATP_points(:,2),'o','markersize',8,'linewidth',2,'color',[0.2667, 0.6667, 0.6]) 
xlim([4 7])
ylim([0 1.6])
set(gca,'FontSize',14)
set (gca, 'xdir', 'reverse' )
xlabel('pCa');
ylabel('Relative Force');
legend('ATP','dATP ka','dATP ka + kd','dATP ka + k1','dATP ka + k3','ATP data','dATP data')

figure
hold on
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_ATP/max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.2, 0.1333, 0.5333])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_ka/max(SS_Ftotal_fpca_ATP),':','linewidth',3,'color',[0.5, 0.5, 0.5])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_ka_kd_k1/max(SS_Ftotal_fpca_ATP),'-.','linewidth',3,'color',[0.5, 0.5, 0.5])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_ka_k1_k3/max(SS_Ftotal_fpca_ATP),'--','linewidth',3,'color',[0.5, 0.5, 0.5])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_ka_kd_k3/max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.5, 0.5, 0.5])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_ka_kd_k1_k3/max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.8, 0.8, 0.8])
plot(ATP_points(:,1),ATP_points(:,2),'o','markersize',8','linewidth',2,'color',[0.2, 0.1333, 0.5333]) 
plot(dATP_points(:,1),dATP_points(:,2),'o','markersize',8,'linewidth',2,'color',[0.2667, 0.6667, 0.6]) 
xlim([4 7])
ylim([0 1.6])
set(gca,'FontSize',14)
set (gca, 'xdir', 'reverse' )
xlabel('pCa');
ylabel('Relative Force');
legend('ATP','dATP ka','dATP ka + kd + k1','dATP ka + k1 + k3','dATP ka + kd + k3','dATP ka + kd + k1 + k3','ATP data','dATP data')

%% Figure S6
% force pCa
[~, ~, ~, ~, SS_Ftotal_fpca_ATP, ~, ~] = myocyte_model(0, 0, 0, 0, 0, 250, 304.7, 4, 2, 80, 4, 25, 0.4, 85, 900, 4); % ATP force pCa
[~, ~, ~, ~, SS_Ftotal_fpca_dATP_ka_kd_k1, ~, ~] = myocyte_model(0, 0, 0, 0, 100, 538, 540, 5.7, 2, 80, 4, 25, 0.4, 85, 900, 4); % dATP ka kd k1 force pCa
[~, ~, ~, ~, SS_Ftotal_fpca_dATP_ka_kd_k1_k2, ~, ~] = myocyte_model(0, 0, 0, 0, 100, 538, 540, 5.7, 2, 1e6, 4, 25, 0.4, 85, 900, 4); % dATP ka kd k1 k2 force pCa
[~, ~, ~, ~, SS_Ftotal_fpca_dATP_ka_kd_k1_k3, ~, ~] = myocyte_model(0, 0, 0, 0, 100, 538, 540, 5.7, 2, 80, 4, 0.01, 0.4, 85, 900, 4); % dATP ka kd k1 k3 force pCa
[~, ~, ~, ~, SS_Ftotal_fpca_dATP_ka_kd_k1_krecruit, ~, ~] = myocyte_model(0, 0, 0, 0, 100, 538, 540, 5.7, 2, 80, 4, 25, 120, 85, 900, 4); % dATP ka kd k1 krecruit force pCa
[~, ~, ~, ~, SS_Ftotal_fpca_dATP_ka_kd_k1_krecruit_kcoop, ~, ~] = myocyte_model(0, 0, 0, 0, 100, 538, 540, 5.7, 2, 80, 4, 25, 120, 85, 900, 9); % dATP ka kd k1 krecruit kcoop force pCa

% Ktr
[~, ~, ~, ~, ~, Ftotal_ktr_ATP, t_ktr_ATP] = myocyte_model(0, 1, 0, 0, 0, 250, 304.7, 4, 2, 80, 4, 25, 0.4, 85, 900, 4); % ATP force pCa
[~, ~, ~, ~, ~, Ftotal_ktr_dATP_ka_kd_k1, t_ktr_dATP_ka_kd_k1] = myocyte_model(0, 1, 0, 0, 100, 538, 540, 5.7, 2, 80, 4, 25, 0.4, 85, 900, 4); % dATP ka kd k1 Ktr
[~, ~, ~, ~, ~, Ftotal_ktr_dATP_ka_kd_k1_k2, t_ktr_dATP_ka_kd_k1_k2] = myocyte_model(0, 1, 0, 0, 100, 538, 540, 5.7, 2, 1e6, 4, 25, 0.4, 85, 900, 4); % dATP ka kd k1 k2 Ktr
[~, ~, ~, ~, ~, Ftotal_ktr_dATP_ka_kd_k1_k3, t_ktr_dATP_ka_kd_k1_k3] = myocyte_model(0, 1, 0, 0, 100, 538, 540, 5.7, 2, 80, 4, 0.01, 0.4, 85, 900, 4); % dATP ka kd k1 k3 Ktr
[~, ~, ~, ~, ~, Ftotal_ktr_dATP_ka_kd_k1_krecruit, t_ktr_dATP_ka_kd_k1_krecruit] = myocyte_model(0, 1, 0, 0, 100, 538, 540, 5.7, 2, 80, 4, 25, 120, 85, 900, 4); % dATP ka kd k1 krecruit Ktr
[~, ~, ~, ~, ~, Ftotal_ktr_dATP_ka_kd_k1_krecruit_kcoop, t_ktr_dATP_ka_kd_k1_krecruit_kcoop] = myocyte_model(0, 1, 0, 0, 100, 538, 540, 5.7, 2, 80, 4, 25, 120, 85, 900, 9); % dATP ka kd k1 krecruit kcoop Ktr

% Shortening
[T_final_XB_ATP, force_final_ATP, idx_XB_ATP, Shortening_final_ATP, ~, ~, ~] = myocyte_model(0, 2, 0, 0, 0, 250, 304.7, 4, 2, 80, 4, 25, 0.4, 85, 900, 4); % ATP force pCa
[T_final_XB_dATP_ka_kd_k1, force_final_dATP_ka_kd_k1, idx_XB_dATP_ka_kd_k1, Shortening_final_dATP_ka_kd_k1, ~, ~, ~] = myocyte_model(0, 2, 0, 0, 2, 538, 540, 5.7, 2, 80, 4, 25, 0.4, 85, 900, 4); % dATP ka kd k1 shortening
[T_final_XB_dATP_ka_kd_k1_k2, force_final_dATP_ka_kd_k1_k2, idx_XB_dATP_ka_kd_k1_k2, Shortening_final_dATP_ka_kd_k1_k2, ~, ~, ~] = myocyte_model(0, 2, 0, 0, 2, 538, 540, 5.7, 2, 1e6, 4, 25, 0.4, 85, 900, 4); % dATP ka kd k1 k2 shortening
[T_final_XB_dATP_ka_kd_k1_k3, force_final_dATP_ka_kd_k1_k3, idx_XB_dATP_ka_kd_k1_k3, Shortening_final_dATP_ka_kd_k1_k3, ~, ~, ~] = myocyte_model(0, 2, 0, 0, 2, 538, 540, 5.7, 2, 80, 4, 0.01, 0.4, 85, 900, 4); % dATP ka kd k1 k3 shortening
[T_final_XB_dATP_ka_kd_k1_krecruit, force_final_dATP_ka_kd_k1_krecruit, idx_XB_dATP_ka_kd_k1_krecruit, Shortening_final_dATP_ka_kd_k1_krecruit, ~, ~, ~] = myocyte_model(0, 2, 0, 0, 2, 538, 540, 5.7, 2, 80, 4, 25, 120, 85, 900, 4); % dATP ka kd k1 krecruit shortening
[T_final_XB_dATP_ka_kd_k1_krecruit_kcoop, force_final_dATP_ka_kd_k1_krecruit_kcoop, idx_XB_dATP_ka_kd_k1_krecruit_kcoop, Shortening_final_dATP_ka_kd_k1_krecruit_kcoop, ~, ~, ~] = myocyte_model(0, 2, 0, 0, 2, 538, 540, 5.7, 2, 80, 4, 25, 120, 85, 900, 9); % dATP ka kd k1 krecruit kcoop shortening

last_beat = 2;

figure
subplot(5,3,1)
hold on
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_ATP/max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.2, 0.1333, 0.5333])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_dATP_ka_kd_k1./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.2667, 0.6667, 0.6])
plot(ATP_points(:,1),ATP_points(:,2),'o','markersize',8','linewidth',2,'color',[0.2, 0.1333, 0.5333]) 
plot(dATP_points(:,1),dATP_points(:,2),'o','markersize',8,'linewidth',2,'color',[0.2667, 0.6667, 0.6]) 
xlim([4 7])
ylim([0 1.6])
set(gca,'FontSize',14)
set (gca, 'xdir', 'reverse' )
xlabel('pCa');
ylabel('Relative Force');
legend('ATP model','dATP model','ATP data','dATP data')

subplot(5,3,2)
hold on
tspan = 0:0.001:2;
plot(t_ktr_ATP-4.0580,Ftotal_ktr_ATP./max(Ftotal_ktr_ATP),'linewidth',3,'color',[0.2, 0.1333, 0.5333]) 
plot(t_ktr_dATP_ka_kd_k1-4.0580,Ftotal_ktr_dATP_ka_kd_k1./max(Ftotal_ktr_ATP),'linewidth',3,'color',[0.2667, 0.6667, 0.6])
plot(tspan,ATP_pCa_4_norm,'-.','linewidth',3,'color',[0.2, 0.1333, 0.5333])
plot(tspan,dATP_pCa_4_norm,'-.','linewidth',3,'color',[0.2667, 0.6667, 0.6])
xlabel('Time (s)')
ylabel('Normalized force')
legend('ATP model','dATP model','ATP data','dATP data')
xlim([0 2])
ylim([0 1.6])
set(gca,'FontSize',14)

subplot(5,3,3)
hold on
plot(T_final_XB_ATP(idx_XB_ATP:end)/1000-last_beat,Shortening_final_ATP(idx_XB_ATP:end)./max(Shortening_final_ATP(idx_XB_ATP:end)),'linewidth',3,'color',[0.2, 0.1333, 0.5333])
plot(T_final_XB_dATP_ka_kd_k1(idx_XB_dATP_ka_kd_k1:end)/1000-last_beat,Shortening_final_dATP_ka_kd_k1(idx_XB_dATP_ka_kd_k1:end)./max(Shortening_final_dATP_ka_kd_k1(idx_XB_dATP_ka_kd_k1:end)),'linewidth',3,'color',[0.2667, 0.6667, 0.6])
xlabel('Time (s)')
ylabel('Relative Shortening')
legend('ATP model','dATP model')
ylim([0.87 1])
set(gca,'FontSize',14)

subplot(5,3,4)
hold on
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_ATP/max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.2, 0.1333, 0.5333])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_dATP_ka_kd_k1_k2./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.2667, 0.6667, 0.6])
plot(ATP_points(:,1),ATP_points(:,2),'o','markersize',8','linewidth',2,'color',[0.2, 0.1333, 0.5333]) 
plot(dATP_points(:,1),dATP_points(:,2),'o','markersize',8,'linewidth',2,'color',[0.2667, 0.6667, 0.6]) 
xlim([4 7])
ylim([0 1.6])
set(gca,'FontSize',14)
set (gca, 'xdir', 'reverse' )
xlabel('pCa');
ylabel('Relative Force');
legend('ATP model','dATP model','ATP data','dATP data')

subplot(5,3,5)
hold on
tspan = 0:0.001:2;
plot(t_ktr_ATP-4.0580,Ftotal_ktr_ATP./max(Ftotal_ktr_ATP),'linewidth',3,'color',[0 0 0])
plot(t_ktr_dATP_ka_kd_k1_k2-4.0580,Ftotal_ktr_dATP_ka_kd_k1_k2./max(Ftotal_ktr_ATP),'linewidth',3,'color',[0.2667, 0.6667, 0.6]) 
plot(tspan,ATP_pCa_4_norm,'-.','linewidth',3,'color',[0.2, 0.1333, 0.5333])
plot(tspan,dATP_pCa_4_norm,'-.','linewidth',3,'color',[0.2667, 0.6667, 0.6])
xlabel('Time (s)')
ylabel('Normalized force')
legend('ATP model','dATP model','ATP data','dATP data')
xlim([0 2])
ylim([0 1.6])
set(gca,'FontSize',14)

subplot(5,3,6)
hold on
plot(T_final_XB_ATP(idx_XB_ATP:end)/1000-last_beat,Shortening_final_ATP(idx_XB_ATP:end)./max(Shortening_final_ATP(idx_XB_ATP:end)),'linewidth',3,'color',[0.2, 0.1333, 0.5333])
plot(T_final_XB_dATP_ka_kd_k1_k2(idx_XB_dATP_ka_kd_k1_k2:end)/1000-last_beat,Shortening_final_dATP_ka_kd_k1_k2(idx_XB_dATP_ka_kd_k1_k2:end)./max(Shortening_final_dATP_ka_kd_k1_k2(idx_XB_dATP_ka_kd_k1_k2:end)),'linewidth',3,'color',[0.2667, 0.6667, 0.6])
xlabel('Time (s)')
ylabel('Relative Shortening')
legend('ATP model','dATP model')
ylim([0.87 1])
set(gca,'FontSize',14)

subplot(5,3,7)
hold on
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_ATP/max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.2, 0.1333, 0.5333])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_dATP_ka_kd_k1_k3./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.2667, 0.6667, 0.6])
plot(ATP_points(:,1),ATP_points(:,2),'o','markersize',8','linewidth',2,'color',[0.2, 0.1333, 0.5333]) 
plot(dATP_points(:,1),dATP_points(:,2),'o','markersize',8,'linewidth',2,'color',[0.2667, 0.6667, 0.6]) 
xlim([4 7])
ylim([0 1.6])
set(gca,'FontSize',14)
set (gca, 'xdir', 'reverse' )
xlabel('pCa');
ylabel('Relative Force');
legend('ATP model','dATP model','ATP data','dATP data')

subplot(5,3,8)
hold on
tspan = 0:0.001:2;
plot(t_ktr_ATP-4.0580,Ftotal_ktr_ATP./max(Ftotal_ktr_ATP),'linewidth',3,'color',[0.2, 0.1333, 0.5333]) 
plot(t_ktr_dATP_ka_kd_k1_k3-4.0580,Ftotal_ktr_dATP_ka_kd_k1_k3./max(Ftotal_ktr_ATP),'linewidth',3,'color',[0.2667, 0.6667, 0.6]) 
plot(tspan,ATP_pCa_4_norm,'-.','linewidth',3,'color',[0.2, 0.1333, 0.5333])
plot(tspan,dATP_pCa_4_norm,'-.','linewidth',3,'color',[0.2667, 0.6667, 0.6])
xlabel('Time (s)')
ylabel('Normalized force')
legend('ATP model','dATP model','ATP data','dATP data')
xlim([0 2])
ylim([0 1.6])
set(gca,'FontSize',14)

subplot(5,3,9)
hold on
plot(T_final_XB_ATP(idx_XB_ATP:end)/1000-last_beat,Shortening_final_ATP(idx_XB_ATP:end)./max(Shortening_final_ATP(idx_XB_ATP:end)),'linewidth',3,'color',[0.2, 0.1333, 0.5333])
plot(T_final_XB_dATP_ka_kd_k1_k3(idx_XB_dATP_ka_kd_k1_k3:end)/1000-last_beat,Shortening_final_dATP_ka_kd_k1_k3(idx_XB_dATP_ka_kd_k1_k3:end)./max(Shortening_final_dATP_ka_kd_k1_k3(idx_XB_dATP_ka_kd_k1_k3:end)),'linewidth',3,'color',[0.2667, 0.6667, 0.6])
xlabel('Time (s)')
ylabel('Relative Shortening')
legend('ATP model','dATP model')
ylim([0.87 1])
set(gca,'FontSize',14)

subplot(5,3,10)
hold on
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_ATP/max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.2, 0.1333, 0.5333])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_dATP_ka_kd_k1_krecruit./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.2667, 0.6667, 0.6])
plot(ATP_points(:,1),ATP_points(:,2),'o','markersize',8','linewidth',2,'color',[0.2, 0.1333, 0.5333]) 
plot(dATP_points(:,1),dATP_points(:,2),'o','markersize',8,'linewidth',2,'color',[0.2667, 0.6667, 0.6]) 
xlim([4 7])
ylim([0 1.6])
set(gca,'FontSize',14)
set (gca, 'xdir', 'reverse' )
xlabel('pCa');
ylabel('Relative Force');
legend('ATP model','dATP model','ATP data','dATP data')

subplot(5,3,11)
hold on
tspan = 0:0.001:2;
plot(t_ktr_ATP-4.0580,Ftotal_ktr_ATP./max(Ftotal_ktr_ATP),'linewidth',3,'color',[0.2, 0.1333, 0.5333]) 
plot(t_ktr_dATP_ka_kd_k1_krecruit-4.0580,Ftotal_ktr_dATP_ka_kd_k1_krecruit./max(Ftotal_ktr_ATP),'linewidth',3,'color',[0.2667, 0.6667, 0.6]) 
plot(tspan,ATP_pCa_4_norm,'-.','linewidth',3,'color',[0.2, 0.1333, 0.5333])
plot(tspan,dATP_pCa_4_norm,'-.','linewidth',3,'color',[0.2667, 0.6667, 0.6])
xlabel('Time (s)')
ylabel('Normalized force')
legend('ATP model','dATP model','ATP data','dATP data')
xlim([0 2])
ylim([0 1.6])
set(gca,'FontSize',14)

subplot(5,3,12)
hold on
plot(T_final_XB_ATP(idx_XB_ATP:end)/1000-last_beat,Shortening_final_ATP(idx_XB_ATP:end)./max(Shortening_final_ATP(idx_XB_ATP:end)),'linewidth',3,'color',[0.2, 0.1333, 0.5333])
plot(T_final_XB_dATP_ka_kd_k1_krecruit(idx_XB_dATP_ka_kd_k1_krecruit:end)/1000-last_beat,Shortening_final_dATP_ka_kd_k1_krecruit(idx_XB_dATP_ka_kd_k1_krecruit:end)./max(Shortening_final_dATP_ka_kd_k1_krecruit(idx_XB_dATP_ka_kd_k1_krecruit:end)),'linewidth',3,'color',[0.2667, 0.6667, 0.6])
xlabel('Time (s)')
ylabel('Relative Shortening')
legend('ATP model','dATP model')
ylim([0.87 1])
set(gca,'FontSize',14)

subplot(5,3,13)
hold on
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_ATP/max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.2, 0.1333, 0.5333])
plot(-log10(Ca_fraction_fpca*(10^(-6))),SS_Ftotal_fpca_dATP_ka_kd_k1_krecruit_kcoop./max(SS_Ftotal_fpca_ATP),'linewidth',3,'color',[0.2667, 0.6667, 0.6])
plot(ATP_points(:,1),ATP_points(:,2),'o','markersize',8','linewidth',2,'color',[0.2, 0.1333, 0.5333]) 
plot(dATP_points(:,1),dATP_points(:,2),'o','markersize',8,'linewidth',2,'color',[0.2667, 0.6667, 0.6]) 
xlim([4 7])
ylim([0 1.6])
set(gca,'FontSize',14)
set (gca, 'xdir', 'reverse' )
xlabel('pCa');
ylabel('Relative Force');
legend('ATP model','dATP model','ATP data','dATP data')

subplot(5,3,14)
hold on
tspan = 0:0.001:2;
plot(t_ktr_ATP-4.0580,Ftotal_ktr_ATP./max(Ftotal_ktr_ATP),'linewidth',3,'color',[0.2, 0.1333, 0.5333]) 
plot(t_ktr_dATP_ka_kd_k1_krecruit_kcoop-4.0580,Ftotal_ktr_dATP_ka_kd_k1_krecruit_kcoop./max(Ftotal_ktr_ATP),'linewidth',3,'color',[0.2667, 0.6667, 0.6]) 
plot(tspan,ATP_pCa_4_norm,'-.','linewidth',3,'color',[0.2, 0.1333, 0.5333])
plot(tspan,dATP_pCa_4_norm,'-.','linewidth',3,'color',[0.2667, 0.6667, 0.6])
xlabel('Time (s)')
ylabel('Normalized force')
legend('ATP model','dATP model','ATP data','dATP data')
xlim([0 2])
ylim([0 1.6])
set(gca,'FontSize',14)

subplot(5,3,15)
hold on
plot(T_final_XB_ATP(idx_XB_ATP:end)/1000-last_beat,Shortening_final_ATP(idx_XB_ATP:end)./max(Shortening_final_ATP(idx_XB_ATP:end)),'linewidth',3,'color',[0.2, 0.1333, 0.5333])
plot(T_final_XB_dATP_ka_kd_k1_krecruit_kcoop(idx_XB_dATP_ka_kd_k1_krecruit_kcoop:end)/1000-last_beat,Shortening_final_dATP_ka_kd_k1_krecruit_kcoop(idx_XB_dATP_ka_kd_k1_krecruit_kcoop:end)./max(Shortening_final_dATP_ka_kd_k1_krecruit_kcoop(idx_XB_dATP_ka_kd_k1_krecruit_kcoop:end)),'linewidth',3,'color',[0.2667, 0.6667, 0.6])
xlabel('Time (s)')
ylabel('Relative Shortening')
legend('ATP model','dATP model')
ylim([0.87 1])
set(gca,'FontSize',14)

%% Figure S7
dATP_sim = [0 2 4 6 8 10 20 30 40 50 60 70 80 90 100];
last_beat = 2000;

% ATP
for i = 1:length(dATP_sim)
dATP = dATP_sim(i);
[T_final_XB_ATP_s, force_final_ATP_s, idx_XB_ATP_s, Shortening_final_ATP_s, ~, ~, ~] = myocyte_model(0, 2, 0, 0, dATP, 250, 304.7, 4, 2, 80, 4, 25, 0.4, 85, 900, 4); % dATP ka shortening
Force_final_ATPa{i} = force_final_ATP_s;
T_final_XB_APTa{i} = T_final_XB_ATP_s;
idx_XB_ATPa{i} = idx_XB_ATP_s;
Shortening_final_ATPa{i} = Shortening_final_ATP_s;
end

force_final_last = Force_final_ATPa{2}(idx_XB_ATPa{2}:end);
T_XB_final_last = T_final_XB_APTa{2}(idx_XB_ATPa{2}:end);
Shortening_final_last = Shortening_final_ATPa{2}(idx_XB_ATPa{2}:end);
peak_shortening = min(Shortening_final_last./Shortening_final_last(1));
FS = -((min(Shortening_final_last)-max(Shortening_final_last))/max(Shortening_final_last))*100;
[min_SL,idx_min] = min(Shortening_final_last);
TTP_shortening = T_XB_final_last(idx_min) - last_beat;

RT50 = max(Shortening_final_last) - (max(Shortening_final_last)-min(Shortening_final_last))*0.5;
time_RT50_shortening = 0;
for i = idx_min:length(Shortening_final_last)
    if Shortening_final_last(i) > RT50
        time_RT50_shortening = ((T_XB_final_last(i)+T_XB_final_last(i-1))/2);
        time_RT50_shortening = time_RT50_shortening - last_beat - TTP_shortening;
        break
    end
end

RT90 = max(Shortening_final_last) - (max(Shortening_final_last)-min(Shortening_final_last))*0.1;
time_RT90_shortening = 0;
for i = idx_min:length(Shortening_final_last)
    if Shortening_final_last(i) > RT90
        time_RT90_shortening = ((T_XB_final_last(i)+T_XB_final_last(i-1))/2);
        time_RT90_shortening = time_RT90_shortening - last_beat - TTP_shortening;
        break
    end
end

FS_store(1) = FS;
TTP_store(1) = TTP_shortening;
RT50_store(1) = time_RT50_shortening;
RT90_store(1) = time_RT90_shortening;

% ka
for i = 1:length(dATP_sim)
dATP = dATP_sim(i);
[T_final_XB_dATP_ka_s, force_final_dATP_ka_s, idx_XB_dATP_ka_s, Shortening_final_dATP_ka_s, ~, ~, ~] = myocyte_model(0, 2, 0, 0, dATP, 538, 304.7, 4, 2, 80, 4, 25, 0.4, 85, 900, 4); % dATP ka shortening
Force_final_dATP_kaa{i} = force_final_dATP_ka_s;
T_final_XB_dATP_kaa{i} = T_final_XB_dATP_ka_s;
idx_XB_dATP_kaa{i} = idx_XB_dATP_ka_s;
Shortening_final_dATP_kaa{i} = Shortening_final_dATP_ka_s;
end

force_final_last = Force_final_dATP_kaa{2}(idx_XB_dATP_kaa{2}:end);
T_XB_final_last = T_final_XB_dATP_kaa{2}(idx_XB_dATP_kaa{2}:end);
Shortening_final_last = Shortening_final_dATP_kaa{2}(idx_XB_dATP_kaa{2}:end);
peak_shortening = min(Shortening_final_last./Shortening_final_last(1));
FS = -((min(Shortening_final_last)-max(Shortening_final_last))/max(Shortening_final_last))*100;
[min_SL,idx_min] = min(Shortening_final_last);
TTP_shortening = T_XB_final_last(idx_min) - last_beat;

RT50 = max(Shortening_final_last) - (max(Shortening_final_last)-min(Shortening_final_last))*0.5;
time_RT50_shortening = 0;
for i = idx_min:length(Shortening_final_last)
    if Shortening_final_last(i) > RT50
        time_RT50_shortening = ((T_XB_final_last(i)+T_XB_final_last(i-1))/2);
        time_RT50_shortening = time_RT50_shortening - last_beat - TTP_shortening;
        break
    end
end

RT90 = max(Shortening_final_last) - (max(Shortening_final_last)-min(Shortening_final_last))*0.1;
time_RT90_shortening = 0;
for i = idx_min:length(Shortening_final_last)
    if Shortening_final_last(i) > RT90
        time_RT90_shortening = ((T_XB_final_last(i)+T_XB_final_last(i-1))/2);
        time_RT90_shortening = time_RT90_shortening - last_beat - TTP_shortening;
        break
    end
end

FS_store(2) = FS;
TTP_store(2) = TTP_shortening;
RT50_store(2) = time_RT50_shortening;
RT90_store(2) = time_RT90_shortening;

% XB
for i = 1:length(dATP_sim)
dATP = dATP_sim(i);
[T_final_XB_dATP_XB_s, force_final_dATP_XB_s, idx_XB_dATP_XB_s, Shortening_final_dATP_XB_s, ~, ~, ~] = myocyte_model(0, 2, 0, 0, dATP, 538, 540, 5.7, 2, 80, 4, 25, 0.4, 85, 900, 4); % dATP XB (ka kd k1) shortening
Force_final_dATP_XBa{i} = force_final_dATP_XB_s;
T_final_XB_dATP_XBa{i} = T_final_XB_dATP_XB_s;
idx_XB_dATP_XBa{i} = idx_XB_dATP_XB_s;
Shortening_final_dATP_XBa{i} = Shortening_final_dATP_XB_s;
end

force_final_last = Force_final_dATP_XBa{2}(idx_XB_dATP_XBa{2}:end);
T_XB_final_last = T_final_XB_dATP_XBa{2}(idx_XB_dATP_XBa{2}:end);
Shortening_final_last = Shortening_final_dATP_XBa{2}(idx_XB_dATP_XBa{2}:end);
peak_shortening = min(Shortening_final_last./Shortening_final_last(1));
FS = -((min(Shortening_final_last)-max(Shortening_final_last))/max(Shortening_final_last))*100;
[min_SL,idx_min] = min(Shortening_final_last);
TTP_shortening = T_XB_final_last(idx_min) - last_beat;

RT50 = max(Shortening_final_last) - (max(Shortening_final_last)-min(Shortening_final_last))*0.5;
time_RT50_shortening = 0;
for i = idx_min:length(Shortening_final_last)
    if Shortening_final_last(i) > RT50
        time_RT50_shortening = ((T_XB_final_last(i)+T_XB_final_last(i-1))/2);
        time_RT50_shortening = time_RT50_shortening - last_beat - TTP_shortening;
        break
    end
end

RT90 = max(Shortening_final_last) - (max(Shortening_final_last)-min(Shortening_final_last))*0.1;
time_RT90_shortening = 0;
for i = idx_min:length(Shortening_final_last)
    if Shortening_final_last(i) > RT90
        time_RT90_shortening = ((T_XB_final_last(i)+T_XB_final_last(i-1))/2);
        time_RT90_shortening = time_RT90_shortening - last_beat - TTP_shortening;
        break
    end
end

FS_store(3) = FS;
TTP_store(3) = TTP_shortening;
RT50_store(3) = time_RT50_shortening;
RT90_store(3) = time_RT90_shortening;

% SRX
for i = 1:length(dATP_sim)
dATP = dATP_sim(i);
[T_final_XB_dATP_SRX_s, force_final_dATP_SRX_s, idx_XB_dATP_SRX_s, Shortening_final_dATP_SRX_s, ~, ~, ~] = myocyte_model(0, 2, 0, 0, dATP, 250, 304.7, 4, 2, 80, 4, 25, 350, 85, 900, 4); % dATP SRX (krecruit) shortening
Force_final_dATP_SRXa{i} = force_final_dATP_SRX_s;
T_final_XB_dATP_SRXa{i} = T_final_XB_dATP_SRX_s;
idx_XB_dATP_SRXa{i} = idx_XB_dATP_SRX_s;
Shortening_final_dATP_SRXa{i} = Shortening_final_dATP_SRX_s;
end

force_final_last = Force_final_dATP_SRXa{2}(idx_XB_dATP_SRXa{2}:end);
T_XB_final_last = T_final_XB_dATP_SRXa{2}(idx_XB_dATP_SRXa{2}:end);
Shortening_final_last = Shortening_final_dATP_SRXa{2}(idx_XB_dATP_SRXa{2}:end);
peak_shortening = min(Shortening_final_last./Shortening_final_last(1));
FS = -((min(Shortening_final_last)-max(Shortening_final_last))/max(Shortening_final_last))*100;
[min_SL,idx_min] = min(Shortening_final_last);
TTP_shortening = T_XB_final_last(idx_min) - last_beat;

RT50 = max(Shortening_final_last) - (max(Shortening_final_last)-min(Shortening_final_last))*0.5;
time_RT50_shortening = 0;
for i = idx_min:length(Shortening_final_last)
    if Shortening_final_last(i) > RT50
        time_RT50_shortening = ((T_XB_final_last(i)+T_XB_final_last(i-1))/2);
        time_RT50_shortening = time_RT50_shortening - last_beat - TTP_shortening;
        break
    end
end

RT90 = max(Shortening_final_last) - (max(Shortening_final_last)-min(Shortening_final_last))*0.1;
time_RT90_shortening = 0;
for i = idx_min:length(Shortening_final_last)
    if Shortening_final_last(i) > RT90
        time_RT90_shortening = ((T_XB_final_last(i)+T_XB_final_last(i-1))/2);
        time_RT90_shortening = time_RT90_shortening - last_beat - TTP_shortening;
        break
    end
end

FS_store(4) = FS;
TTP_store(4) = TTP_shortening;
RT50_store(4) = time_RT50_shortening;
RT90_store(4) = time_RT90_shortening;

% Ca
for i = 1:length(dATP_sim)
dATP = dATP_sim(i);
[T_final_XB_dATP_Ca_s, force_final_dATP_Ca_s, idx_XB_dATP_Ca_s, Shortening_final_dATP_Ca_s, ~, ~, ~] = myocyte_model(1, 2, 0, 0, dATP, 250, 304.7, 4, 2, 80, 4, 25, 0.4, 85, 900, 4); % dATP Ca shortening
Force_final_dATP_Caa{i} = force_final_dATP_Ca_s;
T_final_XB_dATP_Caa{i} = T_final_XB_dATP_Ca_s;
idx_XB_dATP_Caa{i} = idx_XB_dATP_Ca_s;
Shortening_final_dATP_Caa{i} = Shortening_final_dATP_Ca_s;
end

force_final_last = Force_final_dATP_Caa{2}(idx_XB_dATP_Caa{2}:end);
T_XB_final_last = T_final_XB_dATP_Caa{2}(idx_XB_dATP_Caa{2}:end);
Shortening_final_last = Shortening_final_dATP_Caa{2}(idx_XB_dATP_Caa{2}:end);
peak_shortening = min(Shortening_final_last./Shortening_final_last(1));
FS = -((min(Shortening_final_last)-max(Shortening_final_last))/max(Shortening_final_last))*100;
[min_SL,idx_min] = min(Shortening_final_last);
TTP_shortening = T_XB_final_last(idx_min) - last_beat;

RT50 = max(Shortening_final_last) - (max(Shortening_final_last)-min(Shortening_final_last))*0.5;
time_RT50_shortening = 0;
for i = idx_min:length(Shortening_final_last)
    if Shortening_final_last(i) > RT50
        time_RT50_shortening = ((T_XB_final_last(i)+T_XB_final_last(i-1))/2);
        time_RT50_shortening = time_RT50_shortening - last_beat - TTP_shortening;
        break
    end
end

RT90 = max(Shortening_final_last) - (max(Shortening_final_last)-min(Shortening_final_last))*0.1;
time_RT90_shortening = 0;
for i = idx_min:length(Shortening_final_last)
    if Shortening_final_last(i) > RT90
        time_RT90_shortening = ((T_XB_final_last(i)+T_XB_final_last(i-1))/2);
        time_RT90_shortening = time_RT90_shortening - last_beat - TTP_shortening;
        break
    end
end

FS_store(5) = FS;
TTP_store(5) = TTP_shortening;
RT50_store(5) = time_RT50_shortening;
RT90_store(5) = time_RT90_shortening;


% XB + SRX
for i = 1:length(dATP_sim)
dATP = dATP_sim(i);
[T_final_XB_dATP_XB_SRX_s, force_final_dATP_XB_SRX_s, idx_XB_dATP_XB_SRX_s, Shortening_final_dATP_XB_SRX_s, ~, ~, ~] = myocyte_model(0, 2, 0, 0, dATP, 538, 540, 5.7, 2, 80, 4, 25, 350, 85, 900, 4); % dATP XB + SRX shortening
Force_final_dATP_SRXa{i} = force_final_dATP_XB_SRX_s;
T_final_XB_dATP_XB_SRXa{i} = T_final_XB_dATP_XB_SRX_s;
idx_XB_dATP_XB_SRXa{i} = idx_XB_dATP_XB_SRX_s;
Shortening_final_dATP_XB_SRXa{i} = Shortening_final_dATP_XB_SRX_s;
end

force_final_last = Force_final_dATP_SRXa{2}(idx_XB_dATP_XB_SRXa{2}:end);
T_XB_final_last = T_final_XB_dATP_XB_SRXa{2}(idx_XB_dATP_XB_SRXa{2}:end);
Shortening_final_last = Shortening_final_dATP_XB_SRXa{2}(idx_XB_dATP_XB_SRXa{2}:end);
peak_shortening = min(Shortening_final_last./Shortening_final_last(1));
FS = -((min(Shortening_final_last)-max(Shortening_final_last))/max(Shortening_final_last))*100;
[min_SL,idx_min] = min(Shortening_final_last);
TTP_shortening = T_XB_final_last(idx_min) - last_beat;

RT50 = max(Shortening_final_last) - (max(Shortening_final_last)-min(Shortening_final_last))*0.5;
time_RT50_shortening = 0;
for i = idx_min:length(Shortening_final_last)
    if Shortening_final_last(i) > RT50
        time_RT50_shortening = ((T_XB_final_last(i)+T_XB_final_last(i-1))/2);
        time_RT50_shortening = time_RT50_shortening - last_beat - TTP_shortening;
        break
    end
end

RT90 = max(Shortening_final_last) - (max(Shortening_final_last)-min(Shortening_final_last))*0.1;
time_RT90_shortening = 0;
for i = idx_min:length(Shortening_final_last)
    if Shortening_final_last(i) > RT90
        time_RT90_shortening = ((T_XB_final_last(i)+T_XB_final_last(i-1))/2);
        time_RT90_shortening = time_RT90_shortening - last_beat - TTP_shortening;
        break
    end
end

FS_store(6) = FS;
TTP_store(6) = TTP_shortening;
RT50_store(6) = time_RT50_shortening;
RT90_store(6) = time_RT90_shortening;

% XB + Ca
for i = 1:length(dATP_sim)
dATP = dATP_sim(i);
[T_final_XB_dATP_XB_Ca_s, force_final_dATP_XB_Ca_s, idx_XB_dATP_XB_Ca_s, Shortening_final_dATP_XB_Ca_s, ~, ~, ~] = myocyte_model(1, 2, 0, 0, dATP, 538, 540, 5.7, 2, 80, 4, 25, 0.4, 85, 900, 4); % dATP XB + Ca shortening
Force_final_dATP_XB_Caa{i} = force_final_dATP_XB_Ca_s;
T_final_XB_dATP_XB_Caa{i} = T_final_XB_dATP_XB_Ca_s;
idx_XB_dATP_XB_Caa{i} = idx_XB_dATP_XB_Ca_s;
Shortening_final_dATP_XB_Caa{i} = Shortening_final_dATP_XB_Ca_s;
end

force_final_last = Force_final_dATP_XB_Caa{2}(idx_XB_dATP_XB_Caa{2}:end);
T_XB_final_last = T_final_XB_dATP_XB_Caa{2}(idx_XB_dATP_XB_Caa{2}:end);
Shortening_final_last = Shortening_final_dATP_XB_Caa{2}(idx_XB_dATP_XB_Caa{2}:end);
peak_shortening = min(Shortening_final_last./Shortening_final_last(1));
FS = -((min(Shortening_final_last)-max(Shortening_final_last))/max(Shortening_final_last))*100;
[min_SL,idx_min] = min(Shortening_final_last);
TTP_shortening = T_XB_final_last(idx_min) - last_beat;

RT50 = max(Shortening_final_last) - (max(Shortening_final_last)-min(Shortening_final_last))*0.5;
time_RT50_shortening = 0;
for i = idx_min:length(Shortening_final_last)
    if Shortening_final_last(i) > RT50
        time_RT50_shortening = ((T_XB_final_last(i)+T_XB_final_last(i-1))/2);
        time_RT50_shortening = time_RT50_shortening - last_beat - TTP_shortening;
        break
    end
end

RT90 = max(Shortening_final_last) - (max(Shortening_final_last)-min(Shortening_final_last))*0.1;
time_RT90_shortening = 0;
for i = idx_min:length(Shortening_final_last)
    if Shortening_final_last(i) > RT90
        time_RT90_shortening = ((T_XB_final_last(i)+T_XB_final_last(i-1))/2);
        time_RT90_shortening = time_RT90_shortening - last_beat - TTP_shortening;
        break
    end
end

FS_store(7) = FS;
TTP_store(7) = TTP_shortening;
RT50_store(7) = time_RT50_shortening;
RT90_store(7) = time_RT90_shortening;

% Ca + SRX
for i = 1:length(dATP_sim)
dATP = dATP_sim(i);
[T_final_XB_dATP_Ca_SRX_s, force_final_dATP_Ca_SRX_s, idx_XB_dATP_Ca_SRX_s, Shortening_final_dATP_Ca_SRX_s, ~, ~, ~] = myocyte_model(1, 2, 0, 0, dATP, 250, 304.7, 4, 2, 80, 4, 25, 350, 85, 900, 4); % dATP Ca + SRX shortening
Force_final_dATP_Ca_SRXa{i} = force_final_dATP_Ca_SRX_s;
T_final_XB_dATP_Ca_SRXa{i} = T_final_XB_dATP_Ca_SRX_s;
idx_XB_dATP_Ca_SRXa{i} = idx_XB_dATP_Ca_SRX_s;
Shortening_final_dATP_Ca_SRXa{i} = Shortening_final_dATP_Ca_SRX_s;
end

force_final_last = Force_final_dATP_Ca_SRXa{2}(idx_XB_dATP_Ca_SRXa{2}:end);
T_XB_final_last = T_final_XB_dATP_Ca_SRXa{2}(idx_XB_dATP_Ca_SRXa{2}:end);
Shortening_final_last = Shortening_final_dATP_Ca_SRXa{2}(idx_XB_dATP_Ca_SRXa{2}:end);
peak_shortening = min(Shortening_final_last./Shortening_final_last(1));
FS = -((min(Shortening_final_last)-max(Shortening_final_last))/max(Shortening_final_last))*100;
[min_SL,idx_min] = min(Shortening_final_last);
TTP_shortening = T_XB_final_last(idx_min) - last_beat;

RT50 = max(Shortening_final_last) - (max(Shortening_final_last)-min(Shortening_final_last))*0.5;
time_RT50_shortening = 0;
for i = idx_min:length(Shortening_final_last)
    if Shortening_final_last(i) > RT50
        time_RT50_shortening = ((T_XB_final_last(i)+T_XB_final_last(i-1))/2);
        time_RT50_shortening = time_RT50_shortening - last_beat - TTP_shortening;
        break
    end
end

RT90 = max(Shortening_final_last) - (max(Shortening_final_last)-min(Shortening_final_last))*0.1;
time_RT90_shortening = 0;
for i = idx_min:length(Shortening_final_last)
    if Shortening_final_last(i) > RT90
        time_RT90_shortening = ((T_XB_final_last(i)+T_XB_final_last(i-1))/2);
        time_RT90_shortening = time_RT90_shortening - last_beat - TTP_shortening;
        break
    end
end

FS_store(8) = FS;
TTP_store(8) = TTP_shortening;
RT50_store(8) = time_RT50_shortening;
RT90_store(8) = time_RT90_shortening;

% XB + Ca + SRX
for i = 1:length(dATP_sim)
dATP = dATP_sim(i);
[T_final_XB_dATP_XB_Ca_SRX_s, force_final_dATP_XB_Ca_SRX_s, idx_XB_dATP_XB_Ca_SRX_s, Shortening_final_dATP_XB_Ca_SRX_s, ~, ~, ~] = myocyte_model(1, 2, 0, 0, dATP, 538, 540, 5.7, 2, 80, 4, 25, 350, 85, 900, 4); % dATP XB + Ca + SRX shortening
Force_final_dATP_XB_Ca_SRXa{i} = force_final_dATP_XB_Ca_SRX_s;
T_final_XB_dATP_XB_Ca_SRXa{i} = T_final_XB_dATP_XB_Ca_SRX_s;
idx_XB_dATP_XB_Ca_SRXa{i} = idx_XB_dATP_XB_Ca_SRX_s;
Shortening_final_dATP_XB_Ca_SRXa{i} = Shortening_final_dATP_XB_Ca_SRX_s;
end

force_final_last = Force_final_dATP_XB_Ca_SRXa{2}(idx_XB_dATP_XB_Ca_SRXa{2}:end);
T_XB_final_last = T_final_XB_dATP_XB_Ca_SRXa{2}(idx_XB_dATP_XB_Ca_SRXa{2}:end);
Shortening_final_last = Shortening_final_dATP_XB_Ca_SRXa{2}(idx_XB_dATP_XB_Ca_SRXa{2}:end);
peak_shortening = min(Shortening_final_last./Shortening_final_last(1));
FS = -((min(Shortening_final_last)-max(Shortening_final_last))/max(Shortening_final_last))*100;
[min_SL,idx_min] = min(Shortening_final_last);
TTP_shortening = T_XB_final_last(idx_min) - last_beat;

RT50 = max(Shortening_final_last) - (max(Shortening_final_last)-min(Shortening_final_last))*0.5;
time_RT50_shortening = 0;
for i = idx_min:length(Shortening_final_last)
    if Shortening_final_last(i) > RT50
        time_RT50_shortening = ((T_XB_final_last(i)+T_XB_final_last(i-1))/2);
        time_RT50_shortening = time_RT50_shortening - last_beat - TTP_shortening;
        break
    end
end

RT90 = max(Shortening_final_last) - (max(Shortening_final_last)-min(Shortening_final_last))*0.1;
time_RT90_shortening = 0;
for i = idx_min:length(Shortening_final_last)
    if Shortening_final_last(i) > RT90
        time_RT90_shortening = ((T_XB_final_last(i)+T_XB_final_last(i-1))/2);
        time_RT90_shortening = time_RT90_shortening - last_beat - TTP_shortening;
        break
    end
end

FS_store(9) = FS;
TTP_store(9) = TTP_shortening;
RT50_store(9) = time_RT50_shortening;
RT90_store(9) = time_RT90_shortening;

% XB + Ca + SRX, coop = 0
for i = 1:length(dATP_sim)
dATP = dATP_sim(i);
[T_final_XB_dATP_XB_Ca_SRX_no_coop_s, force_final_dATP_XB_Ca_SRX_no_coop_s, idx_XB_dATP_XB_Ca_SRX_no_coop_s, Shortening_final_dATP_XB_Ca_SRX_no_coop_s, ~, ~, ~] = myocyte_model(1, 2, 0, 0, dATP, 538, 540, 5.7, 2, 80, 4, 25, 350, 85, 900, 0); % dATP XB + Ca + SRX, kcoop = 0 shortening
Force_final_dATP_XB_Ca_SRXa_no_coop{i} = force_final_dATP_XB_Ca_SRX_no_coop_s;
T_final_XB_dATP_XB_Ca_SRXa_no_coop{i} = T_final_XB_dATP_XB_Ca_SRX_no_coop_s;
idx_XB_dATP_XB_Ca_SRXa_no_coop{i} = idx_XB_dATP_XB_Ca_SRX_no_coop_s;
Shortening_final_dATP_XB_Ca_SRXa_no_coop{i} = Shortening_final_dATP_XB_Ca_SRX_no_coop_s;
end

force_final_last = Force_final_dATP_XB_Ca_SRXa_no_coop{2}(idx_XB_dATP_XB_Ca_SRXa_no_coop{2}:end);
T_XB_final_last = T_final_XB_dATP_XB_Ca_SRXa_no_coop{2}(idx_XB_dATP_XB_Ca_SRXa_no_coop{2}:end);
Shortening_final_last = Shortening_final_dATP_XB_Ca_SRXa_no_coop{2}(idx_XB_dATP_XB_Ca_SRXa_no_coop{2}:end);
peak_shortening = min(Shortening_final_last./Shortening_final_last(1));
FS = -((min(Shortening_final_last)-max(Shortening_final_last))/max(Shortening_final_last))*100;
[min_SL,idx_min] = min(Shortening_final_last);
TTP_shortening = T_XB_final_last(idx_min) - last_beat;

RT50 = max(Shortening_final_last) - (max(Shortening_final_last)-min(Shortening_final_last))*0.5;
time_RT50_shortening = 0;
for i = idx_min:length(Shortening_final_last)
    if Shortening_final_last(i) > RT50
        time_RT50_shortening = ((T_XB_final_last(i)+T_XB_final_last(i-1))/2);
        time_RT50_shortening = time_RT50_shortening - last_beat - TTP_shortening;
        break
    end
end

RT90 = max(Shortening_final_last) - (max(Shortening_final_last)-min(Shortening_final_last))*0.1;
time_RT90_shortening = 0;
for i = idx_min:length(Shortening_final_last)
    if Shortening_final_last(i) > RT90
        time_RT90_shortening = ((T_XB_final_last(i)+T_XB_final_last(i-1))/2);
        time_RT90_shortening = time_RT90_shortening - last_beat - TTP_shortening;
        break
    end
end

FS_store(10) = FS;
TTP_store(10) = TTP_shortening;
RT50_store(10) = time_RT50_shortening;
RT90_store(10) = time_RT90_shortening;

% A-I
last_beat = 2;

figure
subplot(3,3,1)
hold on
plot(T_final_XB_dATP_kaa{1}(idx_XB_dATP_kaa{1}:end)/1000-last_beat,Shortening_final_dATP_kaa{1}(idx_XB_dATP_kaa{1}:end)./max(Shortening_final_dATP_kaa{1}(idx_XB_dATP_kaa{1}:end)),'linewidth',3,'color',[0.2+0.0074*0, 0.1333+0.0593*0, 0.5333+0.0074*0])
plot(T_final_XB_dATP_kaa{2}(idx_XB_dATP_kaa{2}:end)/1000-last_beat,Shortening_final_dATP_kaa{2}(idx_XB_dATP_kaa{1}:end)./max(Shortening_final_dATP_kaa{2}(idx_XB_dATP_kaa{2}:end)),'linewidth',3,'color',[0.2+0.0074*1, 0.1333+0.0593*1, 0.5333+0.0074*1])
plot(T_final_XB_dATP_kaa{3}(idx_XB_dATP_kaa{3}:end)/1000-last_beat,Shortening_final_dATP_kaa{3}(idx_XB_dATP_kaa{3}:end)./max(Shortening_final_dATP_kaa{3}(idx_XB_dATP_kaa{3}:end)),'linewidth',3,'color',[0.2+0.0074*2, 0.1333+0.0593*2, 0.5333+0.0074*2])
plot(T_final_XB_dATP_kaa{4}(idx_XB_dATP_kaa{4}:end)/1000-last_beat,Shortening_final_dATP_kaa{4}(idx_XB_dATP_kaa{4}:end)./max(Shortening_final_dATP_kaa{4}(idx_XB_dATP_kaa{4}:end)),'linewidth',3,'color',[0.2+0.0074*3, 0.1333+0.0593*3, 0.5333+0.0074*3])
plot(T_final_XB_dATP_kaa{5}(idx_XB_dATP_kaa{5}:end)/1000-last_beat,Shortening_final_dATP_kaa{5}(idx_XB_dATP_kaa{5}:end)./max(Shortening_final_dATP_kaa{5}(idx_XB_dATP_kaa{5}:end)),'linewidth',3,'color',[0.2+0.0074*4, 0.1333+0.0593*4, 0.5333+0.0074*4])
plot(T_final_XB_dATP_kaa{6}(idx_XB_dATP_kaa{6}:end)/1000-last_beat,Shortening_final_dATP_kaa{6}(idx_XB_dATP_kaa{6}:end)./max(Shortening_final_dATP_kaa{6}(idx_XB_dATP_kaa{6}:end)),'linewidth',3,'color',[0.2+0.0074*5, 0.1333+0.0593*5, 0.5333+0.0074*5])
plot(T_final_XB_dATP_kaa{7}(idx_XB_dATP_kaa{7}:end)/1000-last_beat,Shortening_final_dATP_kaa{7}(idx_XB_dATP_kaa{7}:end)./max(Shortening_final_dATP_kaa{7}(idx_XB_dATP_kaa{7}:end)),'linewidth',3,'color',[0.2+0.0074*6, 0.1333+0.0593*6, 0.5333+0.0074*6])
plot(T_final_XB_dATP_kaa{8}(idx_XB_dATP_kaa{8}:end)/1000-last_beat,Shortening_final_dATP_kaa{8}(idx_XB_dATP_kaa{8}:end)./max(Shortening_final_dATP_kaa{8}(idx_XB_dATP_kaa{8}:end)),'linewidth',3,'color',[0.2+0.0074*7, 0.1333+0.0593*7, 0.5333+0.0074*7])
plot(T_final_XB_dATP_kaa{9}(idx_XB_dATP_kaa{9}:end)/1000-last_beat,Shortening_final_dATP_kaa{9}(idx_XB_dATP_kaa{9}:end)./max(Shortening_final_dATP_kaa{9}(idx_XB_dATP_kaa{9}:end)),'linewidth',3,'color',[0.2+0.0074*8, 0.1333+0.0593*8, 0.5333+0.0074*8])
plot(T_final_XB_dATP_kaa{10}(idx_XB_dATP_kaa{10}:end)/1000-last_beat,Shortening_final_dATP_kaa{10}(idx_XB_dATP_kaa{10}:end)./max(Shortening_final_dATP_kaa{10}(idx_XB_dATP_kaa{10}:end)),'linewidth',3,'color',[0.2+0.0074*9, 0.1333+0.0593*9, 0.5333+0.0074*9])
legend('ATP','2% dATP','4% dATP','6% dATP','8% dATP','10% dATP','20% dATP','50% dATP','70% dATP','100% dATP')
set(gca, 'fontsize', 14)
xlabel('Time (s)')
ylabel('Relative shortening')
title('k_a')
ylim([0.7, 1])

subplot(3,3,2)
hold on
plot(T_final_XB_dATP_XBa{1}(idx_XB_dATP_XBa{1}:end)/1000-last_beat,Shortening_final_dATP_XBa{1}(idx_XB_dATP_XBa{1}:end)./max(Shortening_final_dATP_XBa{1}(idx_XB_dATP_XBa{1}:end)),'linewidth',3,'color',[0.2+0.0074*0, 0.1333+0.0593*0, 0.5333+0.0074*0])
plot(T_final_XB_dATP_XBa{2}(idx_XB_dATP_XBa{2}:end)/1000-last_beat,Shortening_final_dATP_XBa{2}(idx_XB_dATP_XBa{1}:end)./max(Shortening_final_dATP_XBa{2}(idx_XB_dATP_XBa{2}:end)),'linewidth',3,'color',[0.2+0.0074*1, 0.1333+0.0593*1, 0.5333+0.0074*1])
plot(T_final_XB_dATP_XBa{3}(idx_XB_dATP_XBa{3}:end)/1000-last_beat,Shortening_final_dATP_XBa{3}(idx_XB_dATP_XBa{3}:end)./max(Shortening_final_dATP_XBa{3}(idx_XB_dATP_XBa{3}:end)),'linewidth',3,'color',[0.2+0.0074*2, 0.1333+0.0593*2, 0.5333+0.0074*2])
plot(T_final_XB_dATP_XBa{4}(idx_XB_dATP_XBa{4}:end)/1000-last_beat,Shortening_final_dATP_XBa{4}(idx_XB_dATP_XBa{4}:end)./max(Shortening_final_dATP_XBa{4}(idx_XB_dATP_XBa{4}:end)),'linewidth',3,'color',[0.2+0.0074*3, 0.1333+0.0593*3, 0.5333+0.0074*3])
plot(T_final_XB_dATP_XBa{5}(idx_XB_dATP_XBa{5}:end)/1000-last_beat,Shortening_final_dATP_XBa{5}(idx_XB_dATP_XBa{5}:end)./max(Shortening_final_dATP_XBa{5}(idx_XB_dATP_XBa{5}:end)),'linewidth',3,'color',[0.2+0.0074*4, 0.1333+0.0593*4, 0.5333+0.0074*4])
plot(T_final_XB_dATP_XBa{6}(idx_XB_dATP_XBa{6}:end)/1000-last_beat,Shortening_final_dATP_XBa{6}(idx_XB_dATP_XBa{6}:end)./max(Shortening_final_dATP_XBa{6}(idx_XB_dATP_XBa{6}:end)),'linewidth',3,'color',[0.2+0.0074*5, 0.1333+0.0593*5, 0.5333+0.0074*5])
plot(T_final_XB_dATP_XBa{7}(idx_XB_dATP_XBa{7}:end)/1000-last_beat,Shortening_final_dATP_XBa{7}(idx_XB_dATP_XBa{7}:end)./max(Shortening_final_dATP_XBa{7}(idx_XB_dATP_XBa{7}:end)),'linewidth',3,'color',[0.2+0.0074*6, 0.1333+0.0593*6, 0.5333+0.0074*6])
plot(T_final_XB_dATP_XBa{8}(idx_XB_dATP_XBa{8}:end)/1000-last_beat,Shortening_final_dATP_XBa{8}(idx_XB_dATP_XBa{8}:end)./max(Shortening_final_dATP_XBa{8}(idx_XB_dATP_XBa{8}:end)),'linewidth',3,'color',[0.2+0.0074*7, 0.1333+0.0593*7, 0.5333+0.0074*7])
plot(T_final_XB_dATP_XBa{9}(idx_XB_dATP_XBa{9}:end)/1000-last_beat,Shortening_final_dATP_XBa{9}(idx_XB_dATP_XBa{9}:end)./max(Shortening_final_dATP_XBa{9}(idx_XB_dATP_XBa{9}:end)),'linewidth',3,'color',[0.2+0.0074*8, 0.1333+0.0593*8, 0.5333+0.0074*8])
plot(T_final_XB_dATP_XBa{10}(idx_XB_dATP_XBa{10}:end)/1000-last_beat,Shortening_final_dATP_XBa{10}(idx_XB_dATP_XBa{10}:end)./max(Shortening_final_dATP_XBa{10}(idx_XB_dATP_XBa{10}:end)),'linewidth',3,'color',[0.2+0.0074*9, 0.1333+0.0593*9, 0.5333+0.0074*9])
legend('ATP','2% dATP','4% dATP','6% dATP','8% dATP','10% dATP','20% dATP','50% dATP','70% dATP','100% dATP')
set(gca, 'fontsize', 14)
xlabel('Time (s)')
ylabel('Relative shortening')
title('XB (k_a, k_d, k_1)')
ylim([0.7, 1])

subplot(3,3,3)
hold on
plot(T_final_XB_dATP_SRXa{1}(idx_XB_dATP_SRXa{1}:end)/1000-last_beat,Shortening_final_dATP_SRXa{1}(idx_XB_dATP_SRXa{1}:end)./max(Shortening_final_dATP_SRXa{1}(idx_XB_dATP_SRXa{1}:end)),'linewidth',3,'color',[0.2+0.0074*0, 0.1333+0.0593*0, 0.5333+0.0074*0])
plot(T_final_XB_dATP_SRXa{2}(idx_XB_dATP_SRXa{2}:end)/1000-last_beat,Shortening_final_dATP_SRXa{2}(idx_XB_dATP_SRXa{1}:end)./max(Shortening_final_dATP_SRXa{2}(idx_XB_dATP_SRXa{2}:end)),'linewidth',3,'color',[0.2+0.0074*1, 0.1333+0.0593*1, 0.5333+0.0074*1])
plot(T_final_XB_dATP_SRXa{3}(idx_XB_dATP_SRXa{3}:end)/1000-last_beat,Shortening_final_dATP_SRXa{3}(idx_XB_dATP_SRXa{3}:end)./max(Shortening_final_dATP_SRXa{3}(idx_XB_dATP_SRXa{3}:end)),'linewidth',3,'color',[0.2+0.0074*2, 0.1333+0.0593*2, 0.5333+0.0074*2])
plot(T_final_XB_dATP_SRXa{4}(idx_XB_dATP_SRXa{4}:end)/1000-last_beat,Shortening_final_dATP_SRXa{4}(idx_XB_dATP_SRXa{4}:end)./max(Shortening_final_dATP_SRXa{4}(idx_XB_dATP_SRXa{4}:end)),'linewidth',3,'color',[0.2+0.0074*3, 0.1333+0.0593*3, 0.5333+0.0074*3])
plot(T_final_XB_dATP_SRXa{5}(idx_XB_dATP_SRXa{5}:end)/1000-last_beat,Shortening_final_dATP_SRXa{5}(idx_XB_dATP_SRXa{5}:end)./max(Shortening_final_dATP_SRXa{5}(idx_XB_dATP_SRXa{5}:end)),'linewidth',3,'color',[0.2+0.0074*4, 0.1333+0.0593*4, 0.5333+0.0074*4])
plot(T_final_XB_dATP_SRXa{6}(idx_XB_dATP_SRXa{6}:end)/1000-last_beat,Shortening_final_dATP_SRXa{6}(idx_XB_dATP_SRXa{6}:end)./max(Shortening_final_dATP_SRXa{6}(idx_XB_dATP_SRXa{6}:end)),'linewidth',3,'color',[0.2+0.0074*5, 0.1333+0.0593*5, 0.5333+0.0074*5])
plot(T_final_XB_dATP_SRXa{7}(idx_XB_dATP_SRXa{7}:end)/1000-last_beat,Shortening_final_dATP_SRXa{7}(idx_XB_dATP_SRXa{7}:end)./max(Shortening_final_dATP_SRXa{7}(idx_XB_dATP_SRXa{7}:end)),'linewidth',3,'color',[0.2+0.0074*6, 0.1333+0.0593*6, 0.5333+0.0074*6])
plot(T_final_XB_dATP_SRXa{8}(idx_XB_dATP_SRXa{8}:end)/1000-last_beat,Shortening_final_dATP_SRXa{8}(idx_XB_dATP_SRXa{8}:end)./max(Shortening_final_dATP_SRXa{8}(idx_XB_dATP_SRXa{8}:end)),'linewidth',3,'color',[0.2+0.0074*7, 0.1333+0.0593*7, 0.5333+0.0074*7])
plot(T_final_XB_dATP_SRXa{9}(idx_XB_dATP_SRXa{9}:end)/1000-last_beat,Shortening_final_dATP_SRXa{9}(idx_XB_dATP_SRXa{9}:end)./max(Shortening_final_dATP_SRXa{9}(idx_XB_dATP_SRXa{9}:end)),'linewidth',3,'color',[0.2+0.0074*8, 0.1333+0.0593*8, 0.5333+0.0074*8])
plot(T_final_XB_dATP_SRXa{10}(idx_XB_dATP_SRXa{10}:end)/1000-last_beat,Shortening_final_dATP_SRXa{10}(idx_XB_dATP_SRXa{10}:end)./max(Shortening_final_dATP_SRXa{10}(idx_XB_dATP_SRXa{10}:end)),'linewidth',3,'color',[0.2+0.0074*9, 0.1333+0.0593*9, 0.5333+0.0074*9])
legend('ATP','2% dATP','4% dATP','6% dATP','8% dATP','10% dATP','20% dATP','50% dATP','70% dATP','100% dATP')
set(gca, 'fontsize', 14)
xlabel('Time (s)')
ylabel('Relative shortening')
title('SRX (k_{recruit})')
ylim([0.7, 1])

subplot(3,3,4)
hold on
plot(T_final_XB_dATP_Caa{1}(idx_XB_dATP_Caa{1}:end)/1000-last_beat,Shortening_final_dATP_Caa{1}(idx_XB_dATP_Caa{1}:end)./max(Shortening_final_dATP_Caa{1}(idx_XB_dATP_Caa{1}:end)),'linewidth',3,'color',[0.2+0.0074*0, 0.1333+0.0593*0, 0.5333+0.0074*0])
plot(T_final_XB_dATP_Caa{2}(idx_XB_dATP_Caa{2}:end)/1000-last_beat,Shortening_final_dATP_Caa{2}(idx_XB_dATP_Caa{1}:end)./max(Shortening_final_dATP_Caa{2}(idx_XB_dATP_Caa{2}:end)),'linewidth',3,'color',[0.2+0.0074*1, 0.1333+0.0593*1, 0.5333+0.0074*1])
plot(T_final_XB_dATP_Caa{3}(idx_XB_dATP_Caa{3}:end)/1000-last_beat,Shortening_final_dATP_Caa{3}(idx_XB_dATP_Caa{3}:end)./max(Shortening_final_dATP_Caa{3}(idx_XB_dATP_Caa{3}:end)),'linewidth',3,'color',[0.2+0.0074*2, 0.1333+0.0593*2, 0.5333+0.0074*2])
plot(T_final_XB_dATP_Caa{4}(idx_XB_dATP_Caa{4}:end)/1000-last_beat,Shortening_final_dATP_Caa{4}(idx_XB_dATP_Caa{4}:end)./max(Shortening_final_dATP_Caa{4}(idx_XB_dATP_Caa{4}:end)),'linewidth',3,'color',[0.2+0.0074*3, 0.1333+0.0593*3, 0.5333+0.0074*3])
plot(T_final_XB_dATP_Caa{5}(idx_XB_dATP_Caa{5}:end)/1000-last_beat,Shortening_final_dATP_Caa{5}(idx_XB_dATP_Caa{5}:end)./max(Shortening_final_dATP_Caa{5}(idx_XB_dATP_Caa{5}:end)),'linewidth',3,'color',[0.2+0.0074*4, 0.1333+0.0593*4, 0.5333+0.0074*4])
plot(T_final_XB_dATP_Caa{6}(idx_XB_dATP_Caa{6}:end)/1000-last_beat,Shortening_final_dATP_Caa{6}(idx_XB_dATP_Caa{6}:end)./max(Shortening_final_dATP_Caa{6}(idx_XB_dATP_Caa{6}:end)),'linewidth',3,'color',[0.2+0.0074*5, 0.1333+0.0593*5, 0.5333+0.0074*5])
plot(T_final_XB_dATP_Caa{7}(idx_XB_dATP_Caa{7}:end)/1000-last_beat,Shortening_final_dATP_Caa{7}(idx_XB_dATP_Caa{7}:end)./max(Shortening_final_dATP_Caa{7}(idx_XB_dATP_Caa{7}:end)),'linewidth',3,'color',[0.2+0.0074*6, 0.1333+0.0593*6, 0.5333+0.0074*6])
plot(T_final_XB_dATP_Caa{8}(idx_XB_dATP_Caa{8}:end)/1000-last_beat,Shortening_final_dATP_Caa{8}(idx_XB_dATP_Caa{8}:end)./max(Shortening_final_dATP_Caa{8}(idx_XB_dATP_Caa{8}:end)),'linewidth',3,'color',[0.2+0.0074*7, 0.1333+0.0593*7, 0.5333+0.0074*7])
plot(T_final_XB_dATP_Caa{9}(idx_XB_dATP_Caa{9}:end)/1000-last_beat,Shortening_final_dATP_Caa{9}(idx_XB_dATP_Caa{9}:end)./max(Shortening_final_dATP_Caa{9}(idx_XB_dATP_Caa{9}:end)),'linewidth',3,'color',[0.2+0.0074*8, 0.1333+0.0593*8, 0.5333+0.0074*8])
plot(T_final_XB_dATP_Caa{10}(idx_XB_dATP_Caa{10}:end)/1000-last_beat,Shortening_final_dATP_Caa{10}(idx_XB_dATP_Caa{10}:end)./max(Shortening_final_dATP_Caa{10}(idx_XB_dATP_Caa{10}:end)),'linewidth',3,'color',[0.2+0.0074*9, 0.1333+0.0593*9, 0.5333+0.0074*9])
legend('ATP','2% dATP','4% dATP','6% dATP','8% dATP','10% dATP','20% dATP','50% dATP','70% dATP','100% dATP')
set(gca, 'fontsize', 14)
xlabel('Time (s)')
ylabel('Relative shortening')
title('Ca')
ylim([0.7, 1])

subplot(3,3,5)
hold on
plot(T_final_XB_dATP_XB_SRXa{1}(idx_XB_dATP_XB_SRXa{1}:end)/1000-last_beat,Shortening_final_dATP_XB_SRXa{1}(idx_XB_dATP_XB_SRXa{1}:end)./max(Shortening_final_dATP_XB_SRXa{1}(idx_XB_dATP_XB_SRXa{1}:end)),'linewidth',3,'color',[0.2+0.0074*0, 0.1333+0.0593*0, 0.5333+0.0074*0])
plot(T_final_XB_dATP_XB_SRXa{2}(idx_XB_dATP_XB_SRXa{2}:end)/1000-last_beat,Shortening_final_dATP_XB_SRXa{2}(idx_XB_dATP_XB_SRXa{1}:end)./max(Shortening_final_dATP_XB_SRXa{2}(idx_XB_dATP_XB_SRXa{2}:end)),'linewidth',3,'color',[0.2+0.0074*1, 0.1333+0.0593*1, 0.5333+0.0074*1])
plot(T_final_XB_dATP_XB_SRXa{3}(idx_XB_dATP_XB_SRXa{3}:end)/1000-last_beat,Shortening_final_dATP_XB_SRXa{3}(idx_XB_dATP_XB_SRXa{3}:end)./max(Shortening_final_dATP_XB_SRXa{3}(idx_XB_dATP_XB_SRXa{3}:end)),'linewidth',3,'color',[0.2+0.0074*2, 0.1333+0.0593*2, 0.5333+0.0074*2])
plot(T_final_XB_dATP_XB_SRXa{4}(idx_XB_dATP_XB_SRXa{4}:end)/1000-last_beat,Shortening_final_dATP_XB_SRXa{4}(idx_XB_dATP_XB_SRXa{4}:end)./max(Shortening_final_dATP_XB_SRXa{4}(idx_XB_dATP_XB_SRXa{4}:end)),'linewidth',3,'color',[0.2+0.0074*3, 0.1333+0.0593*3, 0.5333+0.0074*3])
plot(T_final_XB_dATP_XB_SRXa{5}(idx_XB_dATP_XB_SRXa{5}:end)/1000-last_beat,Shortening_final_dATP_XB_SRXa{5}(idx_XB_dATP_XB_SRXa{5}:end)./max(Shortening_final_dATP_XB_SRXa{5}(idx_XB_dATP_XB_SRXa{5}:end)),'linewidth',3,'color',[0.2+0.0074*4, 0.1333+0.0593*4, 0.5333+0.0074*4])
plot(T_final_XB_dATP_XB_SRXa{6}(idx_XB_dATP_XB_SRXa{6}:end)/1000-last_beat,Shortening_final_dATP_XB_SRXa{6}(idx_XB_dATP_XB_SRXa{6}:end)./max(Shortening_final_dATP_XB_SRXa{6}(idx_XB_dATP_XB_SRXa{6}:end)),'linewidth',3,'color',[0.2+0.0074*5, 0.1333+0.0593*5, 0.5333+0.0074*5])
plot(T_final_XB_dATP_XB_SRXa{7}(idx_XB_dATP_XB_SRXa{7}:end)/1000-last_beat,Shortening_final_dATP_XB_SRXa{7}(idx_XB_dATP_XB_SRXa{7}:end)./max(Shortening_final_dATP_XB_SRXa{7}(idx_XB_dATP_XB_SRXa{7}:end)),'linewidth',3,'color',[0.2+0.0074*6, 0.1333+0.0593*6, 0.5333+0.0074*6])
plot(T_final_XB_dATP_XB_SRXa{8}(idx_XB_dATP_XB_SRXa{8}:end)/1000-last_beat,Shortening_final_dATP_XB_SRXa{8}(idx_XB_dATP_XB_SRXa{8}:end)./max(Shortening_final_dATP_XB_SRXa{8}(idx_XB_dATP_XB_SRXa{8}:end)),'linewidth',3,'color',[0.2+0.0074*7, 0.1333+0.0593*7, 0.5333+0.0074*7])
plot(T_final_XB_dATP_XB_SRXa{9}(idx_XB_dATP_XB_SRXa{9}:end)/1000-last_beat,Shortening_final_dATP_XB_SRXa{9}(idx_XB_dATP_XB_SRXa{9}:end)./max(Shortening_final_dATP_XB_SRXa{9}(idx_XB_dATP_XB_SRXa{9}:end)),'linewidth',3,'color',[0.2+0.0074*8, 0.1333+0.0593*8, 0.5333+0.0074*8])
plot(T_final_XB_dATP_XB_SRXa{10}(idx_XB_dATP_XB_SRXa{10}:end)/1000-last_beat,Shortening_final_dATP_XB_SRXa{10}(idx_XB_dATP_XB_SRXa{10}:end)./max(Shortening_final_dATP_XB_SRXa{10}(idx_XB_dATP_XB_SRXa{10}:end)),'linewidth',3,'color',[0.2+0.0074*9, 0.1333+0.0593*9, 0.5333+0.0074*9])
legend('ATP','2% dATP','4% dATP','6% dATP','8% dATP','10% dATP','20% dATP','50% dATP','70% dATP','100% dATP')
set(gca, 'fontsize', 14)
xlabel('Time (s)')
ylabel('Relative shortening')
title('XB + SRX')
ylim([0.7, 1])

subplot(3,3,6)
hold on
plot(T_final_XB_dATP_XB_Caa{1}(idx_XB_dATP_XB_Caa{1}:end)/1000-last_beat,Shortening_final_dATP_XB_Caa{1}(idx_XB_dATP_XB_Caa{1}:end)./max(Shortening_final_dATP_XB_Caa{1}(idx_XB_dATP_XB_Caa{1}:end)),'linewidth',3,'color',[0.2+0.0074*0, 0.1333+0.0593*0, 0.5333+0.0074*0])
plot(T_final_XB_dATP_XB_Caa{2}(idx_XB_dATP_XB_Caa{2}:end)/1000-last_beat,Shortening_final_dATP_XB_Caa{2}(idx_XB_dATP_XB_Caa{1}:end)./max(Shortening_final_dATP_XB_Caa{2}(idx_XB_dATP_XB_Caa{2}:end)),'linewidth',3,'color',[0.2+0.0074*1, 0.1333+0.0593*1, 0.5333+0.0074*1])
plot(T_final_XB_dATP_XB_Caa{3}(idx_XB_dATP_XB_Caa{3}:end)/1000-last_beat,Shortening_final_dATP_XB_Caa{3}(idx_XB_dATP_XB_Caa{3}:end)./max(Shortening_final_dATP_XB_Caa{3}(idx_XB_dATP_XB_Caa{3}:end)),'linewidth',3,'color',[0.2+0.0074*2, 0.1333+0.0593*2, 0.5333+0.0074*2])
plot(T_final_XB_dATP_XB_Caa{4}(idx_XB_dATP_XB_Caa{4}:end)/1000-last_beat,Shortening_final_dATP_XB_Caa{4}(idx_XB_dATP_XB_Caa{4}:end)./max(Shortening_final_dATP_XB_Caa{4}(idx_XB_dATP_XB_Caa{4}:end)),'linewidth',3,'color',[0.2+0.0074*3, 0.1333+0.0593*3, 0.5333+0.0074*3])
plot(T_final_XB_dATP_XB_Caa{5}(idx_XB_dATP_XB_Caa{5}:end)/1000-last_beat,Shortening_final_dATP_XB_Caa{5}(idx_XB_dATP_XB_Caa{5}:end)./max(Shortening_final_dATP_XB_Caa{5}(idx_XB_dATP_XB_Caa{5}:end)),'linewidth',3,'color',[0.2+0.0074*4, 0.1333+0.0593*4, 0.5333+0.0074*4])
plot(T_final_XB_dATP_XB_Caa{6}(idx_XB_dATP_XB_Caa{6}:end)/1000-last_beat,Shortening_final_dATP_XB_Caa{6}(idx_XB_dATP_XB_Caa{6}:end)./max(Shortening_final_dATP_XB_Caa{6}(idx_XB_dATP_XB_Caa{6}:end)),'linewidth',3,'color',[0.2+0.0074*5, 0.1333+0.0593*5, 0.5333+0.0074*5])
plot(T_final_XB_dATP_XB_Caa{7}(idx_XB_dATP_XB_Caa{7}:end)/1000-last_beat,Shortening_final_dATP_XB_Caa{7}(idx_XB_dATP_XB_Caa{7}:end)./max(Shortening_final_dATP_XB_Caa{7}(idx_XB_dATP_XB_Caa{7}:end)),'linewidth',3,'color',[0.2+0.0074*6, 0.1333+0.0593*6, 0.5333+0.0074*6])
plot(T_final_XB_dATP_XB_Caa{8}(idx_XB_dATP_XB_Caa{8}:end)/1000-last_beat,Shortening_final_dATP_XB_Caa{8}(idx_XB_dATP_XB_Caa{8}:end)./max(Shortening_final_dATP_XB_Caa{8}(idx_XB_dATP_XB_Caa{8}:end)),'linewidth',3,'color',[0.2+0.0074*7, 0.1333+0.0593*7, 0.5333+0.0074*7])
plot(T_final_XB_dATP_XB_Caa{9}(idx_XB_dATP_XB_Caa{9}:end)/1000-last_beat,Shortening_final_dATP_XB_Caa{9}(idx_XB_dATP_XB_Caa{9}:end)./max(Shortening_final_dATP_XB_Caa{9}(idx_XB_dATP_XB_Caa{9}:end)),'linewidth',3,'color',[0.2+0.0074*8, 0.1333+0.0593*8, 0.5333+0.0074*8])
plot(T_final_XB_dATP_XB_Caa{10}(idx_XB_dATP_XB_Caa{10}:end)/1000-last_beat,Shortening_final_dATP_XB_Caa{10}(idx_XB_dATP_XB_Caa{10}:end)./max(Shortening_final_dATP_XB_Caa{10}(idx_XB_dATP_XB_Caa{10}:end)),'linewidth',3,'color',[0.2+0.0074*9, 0.1333+0.0593*9, 0.5333+0.0074*9])
legend('ATP','2% dATP','4% dATP','6% dATP','8% dATP','10% dATP','20% dATP','50% dATP','70% dATP','100% dATP')
set(gca, 'fontsize', 14)
xlabel('Time (s)')
ylabel('Relative shortening')
title('XB + Ca')
ylim([0.7, 1])

subplot(3,3,7)
hold on
plot(T_final_XB_dATP_Ca_SRXa{1}(idx_XB_dATP_Ca_SRXa{1}:end)/1000-last_beat,Shortening_final_dATP_Ca_SRXa{1}(idx_XB_dATP_Ca_SRXa{1}:end)./max(Shortening_final_dATP_Ca_SRXa{1}(idx_XB_dATP_Ca_SRXa{1}:end)),'linewidth',3,'color',[0.2+0.0074*0, 0.1333+0.0593*0, 0.5333+0.0074*0])
plot(T_final_XB_dATP_Ca_SRXa{2}(idx_XB_dATP_Ca_SRXa{2}:end)/1000-last_beat,Shortening_final_dATP_Ca_SRXa{2}(idx_XB_dATP_Ca_SRXa{1}:end)./max(Shortening_final_dATP_Ca_SRXa{2}(idx_XB_dATP_Ca_SRXa{2}:end)),'linewidth',3,'color',[0.2+0.0074*1, 0.1333+0.0593*1, 0.5333+0.0074*1])
plot(T_final_XB_dATP_Ca_SRXa{3}(idx_XB_dATP_Ca_SRXa{3}:end)/1000-last_beat,Shortening_final_dATP_Ca_SRXa{3}(idx_XB_dATP_Ca_SRXa{3}:end)./max(Shortening_final_dATP_Ca_SRXa{3}(idx_XB_dATP_Ca_SRXa{3}:end)),'linewidth',3,'color',[0.2+0.0074*2, 0.1333+0.0593*2, 0.5333+0.0074*2])
plot(T_final_XB_dATP_Ca_SRXa{4}(idx_XB_dATP_Ca_SRXa{4}:end)/1000-last_beat,Shortening_final_dATP_Ca_SRXa{4}(idx_XB_dATP_Ca_SRXa{4}:end)./max(Shortening_final_dATP_Ca_SRXa{4}(idx_XB_dATP_Ca_SRXa{4}:end)),'linewidth',3,'color',[0.2+0.0074*3, 0.1333+0.0593*3, 0.5333+0.0074*3])
plot(T_final_XB_dATP_Ca_SRXa{5}(idx_XB_dATP_Ca_SRXa{5}:end)/1000-last_beat,Shortening_final_dATP_Ca_SRXa{5}(idx_XB_dATP_Ca_SRXa{5}:end)./max(Shortening_final_dATP_Ca_SRXa{5}(idx_XB_dATP_Ca_SRXa{5}:end)),'linewidth',3,'color',[0.2+0.0074*4, 0.1333+0.0593*4, 0.5333+0.0074*4])
plot(T_final_XB_dATP_Ca_SRXa{6}(idx_XB_dATP_Ca_SRXa{6}:end)/1000-last_beat,Shortening_final_dATP_Ca_SRXa{6}(idx_XB_dATP_Ca_SRXa{6}:end)./max(Shortening_final_dATP_Ca_SRXa{6}(idx_XB_dATP_Ca_SRXa{6}:end)),'linewidth',3,'color',[0.2+0.0074*5, 0.1333+0.0593*5, 0.5333+0.0074*5])
plot(T_final_XB_dATP_Ca_SRXa{7}(idx_XB_dATP_Ca_SRXa{7}:end)/1000-last_beat,Shortening_final_dATP_Ca_SRXa{7}(idx_XB_dATP_Ca_SRXa{7}:end)./max(Shortening_final_dATP_Ca_SRXa{7}(idx_XB_dATP_Ca_SRXa{7}:end)),'linewidth',3,'color',[0.2+0.0074*6, 0.1333+0.0593*6, 0.5333+0.0074*6])
plot(T_final_XB_dATP_Ca_SRXa{8}(idx_XB_dATP_Ca_SRXa{8}:end)/1000-last_beat,Shortening_final_dATP_Ca_SRXa{8}(idx_XB_dATP_Ca_SRXa{8}:end)./max(Shortening_final_dATP_Ca_SRXa{8}(idx_XB_dATP_Ca_SRXa{8}:end)),'linewidth',3,'color',[0.2+0.0074*7, 0.1333+0.0593*7, 0.5333+0.0074*7])
plot(T_final_XB_dATP_Ca_SRXa{9}(idx_XB_dATP_Ca_SRXa{9}:end)/1000-last_beat,Shortening_final_dATP_Ca_SRXa{9}(idx_XB_dATP_Ca_SRXa{9}:end)./max(Shortening_final_dATP_Ca_SRXa{9}(idx_XB_dATP_Ca_SRXa{9}:end)),'linewidth',3,'color',[0.2+0.0074*8, 0.1333+0.0593*8, 0.5333+0.0074*8])
plot(T_final_XB_dATP_Ca_SRXa{10}(idx_XB_dATP_Ca_SRXa{10}:end)/1000-last_beat,Shortening_final_dATP_Ca_SRXa{10}(idx_XB_dATP_Ca_SRXa{10}:end)./max(Shortening_final_dATP_Ca_SRXa{10}(idx_XB_dATP_Ca_SRXa{10}:end)),'linewidth',3,'color',[0.2+0.0074*9, 0.1333+0.0593*9, 0.5333+0.0074*9])
legend('ATP','2% dATP','4% dATP','6% dATP','8% dATP','10% dATP','20% dATP','50% dATP','70% dATP','100% dATP')
set(gca, 'fontsize', 14)
xlabel('Time (s)')
ylabel('Relative shortening')
title('Ca + SRX')
ylim([0.7, 1])

subplot(3,3,8)
hold on
plot(T_final_XB_dATP_XB_Ca_SRXa{1}(idx_XB_dATP_XB_Ca_SRXa{1}:end)/1000-last_beat,Shortening_final_dATP_XB_Ca_SRXa{1}(idx_XB_dATP_XB_Ca_SRXa{1}:end)./max(Shortening_final_dATP_XB_Ca_SRXa{1}(idx_XB_dATP_XB_Ca_SRXa{1}:end)),'linewidth',3,'color',[0.2+0.0074*0, 0.1333+0.0593*0, 0.5333+0.0074*0])
plot(T_final_XB_dATP_XB_Ca_SRXa{2}(idx_XB_dATP_XB_Ca_SRXa{2}:end)/1000-last_beat,Shortening_final_dATP_XB_Ca_SRXa{2}(idx_XB_dATP_XB_Ca_SRXa{1}:end)./max(Shortening_final_dATP_XB_Ca_SRXa{2}(idx_XB_dATP_XB_Ca_SRXa{2}:end)),'linewidth',3,'color',[0.2+0.0074*1, 0.1333+0.0593*1, 0.5333+0.0074*1])
plot(T_final_XB_dATP_XB_Ca_SRXa{3}(idx_XB_dATP_XB_Ca_SRXa{3}:end)/1000-last_beat,Shortening_final_dATP_XB_Ca_SRXa{3}(idx_XB_dATP_XB_Ca_SRXa{3}:end)./max(Shortening_final_dATP_XB_Ca_SRXa{3}(idx_XB_dATP_XB_Ca_SRXa{3}:end)),'linewidth',3,'color',[0.2+0.0074*2, 0.1333+0.0593*2, 0.5333+0.0074*2])
plot(T_final_XB_dATP_XB_Ca_SRXa{4}(idx_XB_dATP_XB_Ca_SRXa{4}:end)/1000-last_beat,Shortening_final_dATP_XB_Ca_SRXa{4}(idx_XB_dATP_XB_Ca_SRXa{4}:end)./max(Shortening_final_dATP_XB_Ca_SRXa{4}(idx_XB_dATP_XB_Ca_SRXa{4}:end)),'linewidth',3,'color',[0.2+0.0074*3, 0.1333+0.0593*3, 0.5333+0.0074*3])
plot(T_final_XB_dATP_XB_Ca_SRXa{5}(idx_XB_dATP_XB_Ca_SRXa{5}:end)/1000-last_beat,Shortening_final_dATP_XB_Ca_SRXa{5}(idx_XB_dATP_XB_Ca_SRXa{5}:end)./max(Shortening_final_dATP_XB_Ca_SRXa{5}(idx_XB_dATP_XB_Ca_SRXa{5}:end)),'linewidth',3,'color',[0.2+0.0074*4, 0.1333+0.0593*4, 0.5333+0.0074*4])
plot(T_final_XB_dATP_XB_Ca_SRXa{6}(idx_XB_dATP_XB_Ca_SRXa{6}:end)/1000-last_beat,Shortening_final_dATP_XB_Ca_SRXa{6}(idx_XB_dATP_XB_Ca_SRXa{6}:end)./max(Shortening_final_dATP_XB_Ca_SRXa{6}(idx_XB_dATP_XB_Ca_SRXa{6}:end)),'linewidth',3,'color',[0.2+0.0074*5, 0.1333+0.0593*5, 0.5333+0.0074*5])
plot(T_final_XB_dATP_XB_Ca_SRXa{7}(idx_XB_dATP_XB_Ca_SRXa{7}:end)/1000-last_beat,Shortening_final_dATP_XB_Ca_SRXa{7}(idx_XB_dATP_XB_Ca_SRXa{7}:end)./max(Shortening_final_dATP_XB_Ca_SRXa{7}(idx_XB_dATP_XB_Ca_SRXa{7}:end)),'linewidth',3,'color',[0.2+0.0074*6, 0.1333+0.0593*6, 0.5333+0.0074*6])
plot(T_final_XB_dATP_XB_Ca_SRXa{8}(idx_XB_dATP_XB_Ca_SRXa{8}:end)/1000-last_beat,Shortening_final_dATP_XB_Ca_SRXa{8}(idx_XB_dATP_XB_Ca_SRXa{8}:end)./max(Shortening_final_dATP_XB_Ca_SRXa{8}(idx_XB_dATP_XB_Ca_SRXa{8}:end)),'linewidth',3,'color',[0.2+0.0074*7, 0.1333+0.0593*7, 0.5333+0.0074*7])
plot(T_final_XB_dATP_XB_Ca_SRXa{9}(idx_XB_dATP_XB_Ca_SRXa{9}:end)/1000-last_beat,Shortening_final_dATP_XB_Ca_SRXa{9}(idx_XB_dATP_XB_Ca_SRXa{9}:end)./max(Shortening_final_dATP_XB_Ca_SRXa{9}(idx_XB_dATP_XB_Ca_SRXa{9}:end)),'linewidth',3,'color',[0.2+0.0074*8, 0.1333+0.0593*8, 0.5333+0.0074*8])
plot(T_final_XB_dATP_XB_Ca_SRXa{10}(idx_XB_dATP_XB_Ca_SRXa{10}:end)/1000-last_beat,Shortening_final_dATP_XB_Ca_SRXa{10}(idx_XB_dATP_XB_Ca_SRXa{10}:end)./max(Shortening_final_dATP_XB_Ca_SRXa{10}(idx_XB_dATP_XB_Ca_SRXa{10}:end)),'linewidth',3,'color',[0.2+0.0074*9, 0.1333+0.0593*9, 0.5333+0.0074*9])
legend('ATP','2% dATP','4% dATP','6% dATP','8% dATP','10% dATP','20% dATP','50% dATP','70% dATP','100% dATP')
set(gca, 'fontsize', 14)
xlabel('Time (s)')
ylabel('Relative shortening')
title('XB + Ca + SRX')
ylim([0.7, 1])

subplot(3,3,9)
hold on
plot(T_final_XB_dATP_XB_Ca_SRXa_no_coop{1}(idx_XB_dATP_XB_Ca_SRXa_no_coop{1}:end)/1000-last_beat,Shortening_final_dATP_XB_Ca_SRXa_no_coop{1}(idx_XB_dATP_XB_Ca_SRXa_no_coop{1}:end)./max(Shortening_final_dATP_XB_Ca_SRXa_no_coop{1}(idx_XB_dATP_XB_Ca_SRXa_no_coop{1}:end)),'linewidth',3,'color',[0.2+0.0074*1, 0.1333+0.0593*1, 0.5333+0.0074*1])
plot(T_final_XB_dATP_XB_Ca_SRXa_no_coop{2}(idx_XB_dATP_XB_Ca_SRXa_no_coop{2}:end)/1000-last_beat,Shortening_final_dATP_XB_Ca_SRXa_no_coop{2}(idx_XB_dATP_XB_Ca_SRXa_no_coop{1}:end)./max(Shortening_final_dATP_XB_Ca_SRXa_no_coop{2}(idx_XB_dATP_XB_Ca_SRXa_no_coop{2}:end)),'linewidth',3,'color',[0.2+0.0074*1, 0.1333+0.0593*1, 0.5333+0.0074*1])
plot(T_final_XB_dATP_XB_Ca_SRXa_no_coop{3}(idx_XB_dATP_XB_Ca_SRXa_no_coop{3}:end)/1000-last_beat,Shortening_final_dATP_XB_Ca_SRXa_no_coop{3}(idx_XB_dATP_XB_Ca_SRXa_no_coop{3}:end)./max(Shortening_final_dATP_XB_Ca_SRXa_no_coop{3}(idx_XB_dATP_XB_Ca_SRXa_no_coop{3}:end)),'linewidth',3,'color',[0.2+0.0074*2, 0.1333+0.0593*2, 0.5333+0.0074*2])
plot(T_final_XB_dATP_XB_Ca_SRXa_no_coop{4}(idx_XB_dATP_XB_Ca_SRXa_no_coop{4}:end)/1000-last_beat,Shortening_final_dATP_XB_Ca_SRXa_no_coop{4}(idx_XB_dATP_XB_Ca_SRXa_no_coop{4}:end)./max(Shortening_final_dATP_XB_Ca_SRXa_no_coop{4}(idx_XB_dATP_XB_Ca_SRXa_no_coop{4}:end)),'linewidth',3,'color',[0.2+0.0074*3, 0.1333+0.0593*3, 0.5333+0.0074*3])
plot(T_final_XB_dATP_XB_Ca_SRXa_no_coop{5}(idx_XB_dATP_XB_Ca_SRXa_no_coop{5}:end)/1000-last_beat,Shortening_final_dATP_XB_Ca_SRXa_no_coop{5}(idx_XB_dATP_XB_Ca_SRXa_no_coop{5}:end)./max(Shortening_final_dATP_XB_Ca_SRXa_no_coop{5}(idx_XB_dATP_XB_Ca_SRXa_no_coop{5}:end)),'linewidth',3,'color',[0.2+0.0074*4, 0.1333+0.0593*4, 0.5333+0.0074*4])
plot(T_final_XB_dATP_XB_Ca_SRXa_no_coop{6}(idx_XB_dATP_XB_Ca_SRXa_no_coop{6}:end)/1000-last_beat,Shortening_final_dATP_XB_Ca_SRXa_no_coop{6}(idx_XB_dATP_XB_Ca_SRXa_no_coop{6}:end)./max(Shortening_final_dATP_XB_Ca_SRXa_no_coop{6}(idx_XB_dATP_XB_Ca_SRXa_no_coop{6}:end)),'linewidth',3,'color',[0.2+0.0074*5, 0.1333+0.0593*5, 0.5333+0.0074*5])
plot(T_final_XB_dATP_XB_Ca_SRXa_no_coop{7}(idx_XB_dATP_XB_Ca_SRXa_no_coop{7}:end)/1000-last_beat,Shortening_final_dATP_XB_Ca_SRXa_no_coop{7}(idx_XB_dATP_XB_Ca_SRXa_no_coop{7}:end)./max(Shortening_final_dATP_XB_Ca_SRXa_no_coop{7}(idx_XB_dATP_XB_Ca_SRXa_no_coop{7}:end)),'linewidth',3,'color',[0.2+0.0074*6, 0.1333+0.0593*6, 0.5333+0.0074*6])
plot(T_final_XB_dATP_XB_Ca_SRXa_no_coop{8}(idx_XB_dATP_XB_Ca_SRXa_no_coop{8}:end)/1000-last_beat,Shortening_final_dATP_XB_Ca_SRXa_no_coop{8}(idx_XB_dATP_XB_Ca_SRXa_no_coop{8}:end)./max(Shortening_final_dATP_XB_Ca_SRXa_no_coop{8}(idx_XB_dATP_XB_Ca_SRXa_no_coop{8}:end)),'linewidth',3,'color',[0.2+0.0074*7, 0.1333+0.0593*7, 0.5333+0.0074*7])
plot(T_final_XB_dATP_XB_Ca_SRXa_no_coop{9}(idx_XB_dATP_XB_Ca_SRXa_no_coop{9}:end)/1000-last_beat,Shortening_final_dATP_XB_Ca_SRXa_no_coop{9}(idx_XB_dATP_XB_Ca_SRXa_no_coop{9}:end)./max(Shortening_final_dATP_XB_Ca_SRXa_no_coop{9}(idx_XB_dATP_XB_Ca_SRXa_no_coop{9}:end)),'linewidth',3,'color',[0.2+0.0074*8, 0.1333+0.0593*8, 0.5333+0.0074*8])
plot(T_final_XB_dATP_XB_Ca_SRXa_no_coop{10}(idx_XB_dATP_XB_Ca_SRXa_no_coop{10}:end)/1000-last_beat,Shortening_final_dATP_XB_Ca_SRXa_no_coop{10}(idx_XB_dATP_XB_Ca_SRXa_no_coop{10}:end)./max(Shortening_final_dATP_XB_Ca_SRXa_no_coop{10}(idx_XB_dATP_XB_Ca_SRXa_no_coop{10}:end)),'linewidth',3,'color',[0.2+0.0074*9, 0.1333+0.0593*9, 0.5333+0.0074*9])
legend('ATP','2% dATP','4% dATP','6% dATP','8% dATP','10% dATP','20% dATP','50% dATP','70% dATP','100% dATP')
set(gca, 'fontsize', 14)
xlabel('Time (s)')
ylabel('Relative shortening')
title('XB + Ca + SRX, no cooperativity (k_{coop} = 0)')
ylim([0.7, 1])

% K
data = [0.55,-0.03,-0.29,-0.43];
cdata = horzcat(((FS_store-FS_store(1))./FS_store(1))', ((TTP_store-TTP_store(1))./TTP_store(1))', ((RT50_store-RT50_store(1))./RT50_store(1))', ((RT90_store-RT90_store(1))./RT90_store(1))');
cdata1 = vertcat(cdata,data);
xvalues = {'FS','TTP','RT50','RT90'};
yvalues = {'ATP','ka','XB','SRX','Ca','XB + SRX','XB + Ca','Ca + SRX','XB + Ca + SRX','XB + Ca + SRX no coop','data'};

for i = 1:50
    Color(i,1) = 0.6667+0.004*(i-1);
    Color(i,2) = 0.2667+0.01*(i-1);
    Color(i,3) = 0.6-0.003*(i-1);
end

figure
h = heatmap(xvalues,yvalues,cdata1);
h.Colormap = Color;




%% Figure S12
dATP_sim = [0 2 4 6 8 10 20 30 40 50 60 70 80 90 100];

% No SRX state
% pCa = 4
for i = 1:length(dATP_sim)
dATP = dATP_sim(i);
[~, ~, ~, ~, ~, Ftotal_ktr, t_ktr] = myocyte_model(0, 1, 0, 0, dATP, 538, 540, 5.7, 2, 80, 4, 25, 0.4, 85, 900, 4); 
Max_force_pCa_4(i) = max(Ftotal_ktr);
end

% pCa = 4.5
for i = 1:length(dATP_sim)
dATP = dATP_sim(i);
[~, ~, ~, ~, ~, Ftotal_ktr, t_ktr] = myocyte_model(0, 1, 1, 0, dATP, 538, 540, 5.7, 2, 80, 4, 25, 0.4, 85, 900, 4); 
Max_force_pCa_4_5(i) = max(Ftotal_ktr);
end

% pCa = 5
for i = 1:length(dATP_sim)
dATP = dATP_sim(i);
[~, ~, ~, ~, ~, Ftotal_ktr, t_ktr] = myocyte_model(0, 1, 2, 0, dATP, 538, 540, 5.7, 2, 80, 4, 25, 0.4, 85, 900, 4); 
Max_force_pCa_5(i) = max(Ftotal_ktr);
end

% pCa = 5.5
for i = 1:length(dATP_sim)
dATP = dATP_sim(i);
[~, ~, ~, ~, ~, Ftotal_ktr, t_ktr] = myocyte_model(0, 1, 3, 0, dATP, 538, 540, 5.7, 2, 80, 4, 25, 0.4, 85, 900, 4); 
Max_force_pCa_5_5(i) = max(Ftotal_ktr);
end

% pCa = 6
for i = 1:length(dATP_sim)
dATP = dATP_sim(i);
[~, ~, ~, ~, ~, Ftotal_ktr, t_ktr] = myocyte_model(0, 1, 4, 0, dATP, 538, 540, 5.7, 2, 80, 4, 25, 0.4, 85, 900, 4); 
Max_force_pCa_6(i) = max(Ftotal_ktr);
end

% pCa = 6.5
for i = 1:length(dATP_sim)
dATP = dATP_sim(i);
[~, ~, ~, ~, ~, Ftotal_ktr, t_ktr] = myocyte_model(0, 1, 5, 0, dATP, 538, 540, 5.7, 2, 80, 4, 25, 0.4, 85, 900, 4); 
Max_force_pCa_6_5(i) = max(Ftotal_ktr);
end

% pCa = 7
for i = 1:length(dATP_sim)
dATP = dATP_sim(i);
[~, ~, ~, ~, ~, Ftotal_ktr, t_ktr] = myocyte_model(0, 1, 6, 0, dATP, 538, 540, 5.7, 2, 80, 4, 25, 0.4, 85, 900, 4); 
Max_force_pCa_7(i) = max(Ftotal_ktr);
end

figure
hold on
plot(dATP_sim, Max_force_pCa_4./min(Max_force_pCa_4),'linewidth', 3,'color', [0.20, 0.13, 0.53])
plot(dATP_sim, Max_force_pCa_4_5./min(Max_force_pCa_4_5),'linewidth', 3, 'color', [0.52, 0.05, 0.29])
plot(dATP_sim, Max_force_pCa_5./min(Max_force_pCa_5),'linewidth', 3, 'color', [0.87, 0.80, 0.47])
plot(dATP_sim, Max_force_pCa_5_5./min(Max_force_pCa_5_5),'linewidth', 3, 'color', [0.15, 0.58, 0.77])
plot(dATP_sim, Max_force_pCa_6./min(Max_force_pCa_6),'linewidth', 3, 'color', [0.05, 0.42, 0.22])
plot(dATP_sim, Max_force_pCa_6_5./min(Max_force_pCa_6_5),'linewidth', 3, 'color', [0.68, 0.87, 0.94])
plot(dATP_sim, Max_force_pCa_7./min(Max_force_pCa_7),'linewidth', 3, 'color', [0.67, 0.27, 0.60])
xlabel('dATP (%)')
ylabel('Maximum steady state force (normalized)')
legend('pCa = 4','pCa = 4.5','pCa = 5','pCa = 5.5','pCa = 6','pCa = 6.5','pCa = 7')
set(gca, 'fontsize', 14)


% With SRX state
% pCa = 4
for i = 1:length(dATP_sim)
dATP = dATP_sim(i);
[~, ~, ~, ~, ~, Ftotal_ktr, t_ktr] = myocyte_model(0, 1, 0, 0, dATP, 538, 540, 5.7, 2, 80, 4, 25, 350, 85, 900, 4); 
Max_force_pCa_4_SRX(i) = max(Ftotal_ktr);
end

% pCa = 4.5
for i = 1:length(dATP_sim)
dATP = dATP_sim(i);
[~, ~, ~, ~, ~, Ftotal_ktr, t_ktr] = myocyte_model(0, 1, 1, 0, dATP, 538, 540, 5.7, 2, 80, 4, 25, 350, 85, 900, 4); 
Max_force_pCa_4_5_SRX(i) = max(Ftotal_ktr);
end

% pCa = 5
for i = 1:length(dATP_sim)
dATP = dATP_sim(i);
[~, ~, ~, ~, ~, Ftotal_ktr, t_ktr] = myocyte_model(0, 1, 2, 0, dATP, 538, 540, 5.7, 2, 80, 4, 25, 350, 85, 900, 4); 
Max_force_pCa_5_SRX(i) = max(Ftotal_ktr);
end

% pCa = 5.5
for i = 1:length(dATP_sim)
dATP = dATP_sim(i);
[~, ~, ~, ~, ~, Ftotal_ktr, t_ktr] = myocyte_model(0, 1, 3, 0, dATP, 538, 540, 5.7, 2, 80, 4, 25, 350, 85, 900, 4); 
Max_force_pCa_5_5_SRX(i) = max(Ftotal_ktr);
end

% pCa = 6
for i = 1:length(dATP_sim)
dATP = dATP_sim(i);
[~, ~, ~, ~, ~, Ftotal_ktr, t_ktr] = myocyte_model(0, 1, 4, 0, dATP, 538, 540, 5.7, 2, 80, 4, 25, 350, 85, 900, 4); 
Max_force_pCa_6_SRX(i) = max(Ftotal_ktr);
end

% pCa = 6.5
for i = 1:length(dATP_sim)
dATP = dATP_sim(i);
[~, ~, ~, ~, ~, Ftotal_ktr, t_ktr] = myocyte_model(0, 1, 5, 0, dATP, 538, 540, 5.7, 2, 80, 4, 25, 350, 85, 900, 4); 
Max_force_pCa_6_5_SRX(i) = max(Ftotal_ktr);
end

% pCa = 7
for i = 1:length(dATP_sim)
dATP = dATP_sim(i);
[~, ~, ~, ~, ~, Ftotal_ktr, t_ktr] = myocyte_model(0, 1, 6, 0, dATP, 538, 540, 5.7, 2, 80, 4, 25, 350, 85, 900, 4); 
Max_force_pCa_7_SRX(i) = max(Ftotal_ktr);
end

figure
hold on
plot(dATP_sim, Max_force_pCa_4_SRX./min(Max_force_pCa_4_SRX),'linewidth', 3,'color', [0.20, 0.13, 0.53])
plot(dATP_sim, Max_force_pCa_4_5_SRX./min(Max_force_pCa_4_5_SRX),'linewidth', 3, 'color', [0.52, 0.05, 0.29])
plot(dATP_sim, Max_force_pCa_5_SRX./min(Max_force_pCa_5_SRX),'linewidth', 3, 'color', [0.87, 0.80, 0.47])
plot(dATP_sim, Max_force_pCa_5_5_SRX./min(Max_force_pCa_5_5_SRX),'linewidth', 3, 'color', [0.15, 0.58, 0.77])
plot(dATP_sim, Max_force_pCa_6_SRX./min(Max_force_pCa_6_SRX),'linewidth', 3, 'color', [0.05, 0.42, 0.22])
plot(dATP_sim, Max_force_pCa_6_5_SRX./min(Max_force_pCa_6_5_SRX),'linewidth', 3, 'color', [0.68, 0.87, 0.94])
plot(dATP_sim, Max_force_pCa_7_SRX./min(Max_force_pCa_7_SRX),'linewidth', 3, 'color', [0.67, 0.27, 0.60])
xlabel('dATP (%)')
ylabel('Maximum steady state force (normalized)')
legend('pCa = 4','pCa = 4.5','pCa = 5','pCa = 5.5','pCa = 6','pCa = 6.5','pCa = 7')
set(gca, 'fontsize', 14)