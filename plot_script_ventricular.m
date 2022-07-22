%% Figure 4
% Healthy
dATP_sim = [0 2 4 6 8 10 20];
for i = 1:length(dATP_sim)
    dATP = dATP_sim(i);
    [V_LV_store_F, P_LV_store_F, max_force, FS, CO, EF, LVDP, plus_dPdt_F, minus_dPdt_F, work_rate_F, ATP_F, ADP_F, Pi_F, MVO2_F, PCrATP_F, XB_turnover_F, ATPase_F, efficiency_F, SRX_store_F] = cardiovascular_model(1, 0, 0, dATP, 538, 540, 5.7, 2, 80, 4, 25, 350, 280, 800, 4);
    V_LV_store_F_healthy{i} = V_LV_store_F;
    P_LV_store_F_healthy{i} = P_LV_store_F;
    max_force_healthy{i} = max_force;
    FS_healthy{i} = FS;
    CO_healthy{i} = CO;
    EF_healthy{i} = EF;
    LVDP_healthy{i} = LVDP;
    plus_dPdt_F_healthy{i} = plus_dPdt_F;
    minus_dPdt_F_healthy{i} = minus_dPdt_F;
    work_rate_F_healthy{i} = work_rate_F;
    ATP_F_healthy{i} = ATP_F; 
    ADP_F_healthy{i} = ADP_F;
    Pi_F_healthy{i} = Pi_F;
    MVO2_F_healthy{i} = MVO2_F;
    PCrATP_F_healthy{i} = PCrATP_F;
    XB_turnover_F_healthy{i} = XB_turnover_F;
    ATPase_F_healthy{i} = ATPase_F;
    efficiency_F_healthy{i} = efficiency_F;
end

% Failue
for i = 1:length(dATP_sim)
    dATP = dATP_sim(i);
    [V_LV_store_F, P_LV_store_F, max_force, FS, CO, EF, LVDP, plus_dPdt_F, minus_dPdt_F, work_rate_F, ATP_F, ADP_F, Pi_F, MVO2_F, PCrATP_F, XB_turnover_F, ATPase_F, efficiency_F, SRX_store_F] = cardiovascular_model(1, 1, 0, dATP, 538, 540, 5.7, 2, 80, 4, 25, 350, 280, 800, 4);
    V_LV_store_F_failure{i} = V_LV_store_F;
    P_LV_store_F_failure{i} = P_LV_store_F;
    max_force_failure{i} = max_force;
    FS_failure{i} = FS;
    CO_failure{i} = CO;
    EF_failure{i} = EF;
    LVDP_failure{i} = LVDP;
    plus_dPdt_F_failure{i} = plus_dPdt_F;
    minus_dPdt_F_failure{i} = minus_dPdt_F;
    work_rate_F_failure{i} = work_rate_F;
    ATP_F_failure{i} = ATP_F; 
    ADP_F_failure{i} = ADP_F;
    Pi_F_failure{i} = Pi_F;
    MVO2_F_failure{i} = MVO2_F;
    PCrATP_F_failure{i} = PCrATP_F;
    XB_turnover_F_failure{i} = XB_turnover_F;
    ATPase_F_failure{i} = ATPase_F;
    efficiency_F_failure{i} = efficiency_F;
end

% A
figure
hold on
plot(V_LV_store_F_healthy{1}{1},P_LV_store_F_healthy{1}{1}-47,':','linewidth',3,'color',[0.2+0.0111*0 0.1333+0.0889*0 0.5333+0.0111*0])
plot(V_LV_store_F_healthy{2}{1},P_LV_store_F_healthy{2}{1}-47,'linewidth',3,'color',[0.2+0.0111*1 0.1333+0.0889*1 0.5333+0.0111*1])
plot(V_LV_store_F_healthy{3}{1},P_LV_store_F_healthy{3}{1}-47,'linewidth',3,'color',[0.2+0.0111*2 0.1333+0.0889*2 0.5333+0.0111*2])
plot(V_LV_store_F_healthy{4}{1},P_LV_store_F_healthy{4}{1}-47,'linewidth',3,'color',[0.2+0.0111*3 0.1333+0.0889*3 0.5333+0.0111*3])
plot(V_LV_store_F_healthy{5}{1},P_LV_store_F_healthy{5}{1}-47,'linewidth',3,'color',[0.2+0.0111*4 0.1333+0.0889*4 0.5333+0.0111*4])
plot(V_LV_store_F_healthy{6}{1},P_LV_store_F_healthy{6}{1}-47,'linewidth',3,'color',[0.2+0.0111*5 0.1333+0.0889*5 0.5333+0.0111*5])
plot(V_LV_store_F_healthy{7}{1},P_LV_store_F_healthy{7}{1}-47,'linewidth',3,'color',[0.2+0.0111*6 0.1333+0.0889*6 0.5333+0.0111*6])

plot(V_LV_store_F_failure{1}{1},P_LV_store_F_failure{1}{1}-47,':','linewidth',3,'color',[0.6667+0.0299*0 0.2667+0.0888*0 0.6000-0.0222*0])
plot(V_LV_store_F_failure{2}{1},P_LV_store_F_failure{2}{1}-47,'linewidth',3,'color',[0.6667+0.0299*1 0.2667+0.0888*1 0.6000-0.0222*1])
plot(V_LV_store_F_failure{3}{1},P_LV_store_F_failure{3}{1}-47,'linewidth',3,'color',[0.6667+0.0299*2 0.2667+0.0888*2 0.6000-0.0222*2])
plot(V_LV_store_F_failure{4}{1},P_LV_store_F_failure{4}{1}-47,'linewidth',3,'color',[0.6667+0.0299*3 0.2667+0.0888*3 0.6000-0.0222*3])
plot(V_LV_store_F_failure{5}{1},P_LV_store_F_failure{5}{1}-47,'linewidth',3,'color',[0.6667+0.0299*4 0.2667+0.0888*4 0.6000-0.0222*4])
plot(V_LV_store_F_failure{6}{1},P_LV_store_F_failure{6}{1}-47,'linewidth',3,'color',[0.6667+0.0299*5 0.2667+0.0888*5 0.6000-0.0222*5])
plot(V_LV_store_F_failure{7}{1},P_LV_store_F_failure{7}{1}-47,'linewidth',3,'color',[0.6667+0.0299*6 0.2667+0.0888*6 0.6000-0.0222*6])
xlabel('Volume (mL)')
ylabel('Pressure (mmHg)')
legend('ATP','2% dATP','4% dATP','6% dATP','8% dATP','10% dATP','20% dATP','HF','HF + 2% dATP','HF + 4% dATP','HF + 6% dATP','HF + 8% dATP','HF + 10% dATP','HF + 20% dATP')
set(gca, 'fontsize', 14)

% B
percent_dATP = [2 4 6 8 10 20];
figure
subplot(4,4,1)
hold on
plot(0,max_force_healthy{1}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*0 0.1333+0.0889*0 0.5333+0.0111*0])
plot(percent_dATP(1),max_force_healthy{2}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*1 0.1333+0.0889*1 0.5333+0.0111*1])
plot(percent_dATP(2),max_force_healthy{3}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*2 0.1333+0.0889*2 0.5333+0.0111*2])
plot(percent_dATP(3),max_force_healthy{4}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*3 0.1333+0.0889*3 0.5333+0.0111*3])
plot(percent_dATP(4),max_force_healthy{5}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*4 0.1333+0.0889*4 0.5333+0.0111*4])
plot(percent_dATP(5),max_force_healthy{6}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*5 0.1333+0.0889*5 0.5333+0.0111*5])
plot(percent_dATP(6),max_force_healthy{7}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*6 0.1333+0.0889*6 0.5333+0.0111*6])

plot(0,max_force_failure{1}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*0 0.2667+0.0888*0 0.6000-0.0222*0])
plot(percent_dATP(1),max_force_failure{2}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*1 0.2667+0.0888*1 0.6000-0.0222*1])
plot(percent_dATP(2),max_force_failure{3}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*2 0.2667+0.0888*2 0.6000-0.0222*2])
plot(percent_dATP(3),max_force_failure{4}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*3 0.2667+0.0888*3 0.6000-0.0222*3])
plot(percent_dATP(4),max_force_failure{5}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*4 0.2667+0.0888*4 0.6000-0.0222*4])
plot(percent_dATP(5),max_force_failure{6}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*5 0.2667+0.0888*5 0.6000-0.0222*5])
plot(percent_dATP(6),max_force_failure{7}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*6 0.2667+0.0888*6 0.6000-0.0222*6])
xlabel('Percent dATP')
ylabel('Max Force (kPa)')
set(gca, 'fontsize', 14)

subplot(4,4,2)
hold on
plot(0,CO_healthy{1}{1}/454,'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*0 0.1333+0.0889*0 0.5333+0.0111*0])
plot(percent_dATP(1),CO_healthy{2}{1}/454,'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*1 0.1333+0.0889*1 0.5333+0.0111*1])
plot(percent_dATP(2),CO_healthy{3}{1}/454,'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*2 0.1333+0.0889*2 0.5333+0.0111*2])
plot(percent_dATP(3),CO_healthy{4}{1}/454,'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*3 0.1333+0.0889*3 0.5333+0.0111*3])
plot(percent_dATP(4),CO_healthy{5}{1}/454,'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*4 0.1333+0.0889*4 0.5333+0.0111*4])
plot(percent_dATP(5),CO_healthy{6}{1}/454,'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*5 0.1333+0.0889*5 0.5333+0.0111*5])
plot(percent_dATP(6),CO_healthy{7}{1}/454,'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*6 0.1333+0.0889*6 0.5333+0.0111*6])

plot(0,CO_failure{1}{1}/454,'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*0 0.2667+0.0888*0 0.6000-0.0222*0])
plot(percent_dATP(1),CO_failure{2}{1}/454,'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*1 0.2667+0.0888*1 0.6000-0.0222*1])
plot(percent_dATP(2),CO_failure{3}{1}/454,'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*2 0.2667+0.0888*2 0.6000-0.0222*2])
plot(percent_dATP(3),CO_failure{4}{1}/454,'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*3 0.2667+0.0888*3 0.6000-0.0222*3])
plot(percent_dATP(4),CO_failure{5}{1}/454,'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*4 0.2667+0.0888*4 0.6000-0.0222*4])
plot(percent_dATP(5),CO_failure{6}{1}/454,'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*5 0.2667+0.0888*5 0.6000-0.0222*5])
plot(percent_dATP(6),CO_failure{7}{1}/454,'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*6 0.2667+0.0888*6 0.6000-0.0222*6])
xlabel('Percent dATP')
ylabel('SV (mL/beat)')
set(gca, 'fontsize', 14)

subplot(4,4,3)
hold on
plot(0,EF_healthy{1}{1}*100,'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*0 0.1333+0.0889*0 0.5333+0.0111*0])
plot(percent_dATP(1),EF_healthy{2}{1}*100,'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*1 0.1333+0.0889*1 0.5333+0.0111*1])
plot(percent_dATP(2),EF_healthy{3}{1}*100,'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*2 0.1333+0.0889*2 0.5333+0.0111*2])
plot(percent_dATP(3),EF_healthy{4}{1}*100,'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*3 0.1333+0.0889*3 0.5333+0.0111*3])
plot(percent_dATP(4),EF_healthy{5}{1}*100,'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*4 0.1333+0.0889*4 0.5333+0.0111*4])
plot(percent_dATP(5),EF_healthy{6}{1}*100,'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*5 0.1333+0.0889*5 0.5333+0.0111*5])
plot(percent_dATP(6),EF_healthy{7}{1}*100,'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*6 0.1333+0.0889*6 0.5333+0.0111*6])

plot(0,EF_failure{1}{1}*100,'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*0 0.2667+0.0888*0 0.6000-0.0222*0])
plot(percent_dATP(1),EF_failure{2}{1}*100,'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*1 0.2667+0.0888*1 0.6000-0.0222*1])
plot(percent_dATP(2),EF_failure{3}{1}*100,'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*2 0.2667+0.0888*2 0.6000-0.0222*2])
plot(percent_dATP(3),EF_failure{4}{1}*100,'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*3 0.2667+0.0888*3 0.6000-0.0222*3])
plot(percent_dATP(4),EF_failure{5}{1}*100,'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*4 0.2667+0.0888*4 0.6000-0.0222*4])
plot(percent_dATP(5),EF_failure{6}{1}*100,'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*5 0.2667+0.0888*5 0.6000-0.0222*5])
plot(percent_dATP(6),EF_failure{7}{1}*100,'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*6 0.2667+0.0888*6 0.6000-0.0222*6])
xlabel('Percent dATP')
ylabel('EF (%)')
set(gca, 'fontsize', 14)

subplot(4,4,4)
hold on
plot(0,CO_healthy{1}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*0 0.1333+0.0889*0 0.5333+0.0111*0])
plot(percent_dATP(1),CO_healthy{2}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*1 0.1333+0.0889*1 0.5333+0.0111*1])
plot(percent_dATP(2),CO_healthy{3}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*2 0.1333+0.0889*2 0.5333+0.0111*2])
plot(percent_dATP(3),CO_healthy{4}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*3 0.1333+0.0889*3 0.5333+0.0111*3])
plot(percent_dATP(4),CO_healthy{5}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*4 0.1333+0.0889*4 0.5333+0.0111*4])
plot(percent_dATP(5),CO_healthy{6}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*5 0.1333+0.0889*5 0.5333+0.0111*5])
plot(percent_dATP(6),CO_healthy{7}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*6 0.1333+0.0889*6 0.5333+0.0111*6])

plot(0,CO_failure{1}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*0 0.2667+0.0888*0 0.6000-0.0222*0])
plot(percent_dATP(1),CO_failure{2}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*1 0.2667+0.0888*1 0.6000-0.0222*1])
plot(percent_dATP(2),CO_failure{3}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*2 0.2667+0.0888*2 0.6000-0.0222*2])
plot(percent_dATP(3),CO_failure{4}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*3 0.2667+0.0888*3 0.6000-0.0222*3])
plot(percent_dATP(4),CO_failure{5}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*4 0.2667+0.0888*4 0.6000-0.0222*4])
plot(percent_dATP(5),CO_failure{6}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*5 0.2667+0.0888*5 0.6000-0.0222*5])
plot(percent_dATP(6),CO_failure{7}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*6 0.2667+0.0888*6 0.6000-0.0222*6])
xlabel('Percent dATP')
ylabel('CO (mL/min)')
set(gca, 'fontsize', 14)

subplot(4,4,5)
hold on
plot(0,LVDP_healthy{1}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*0 0.1333+0.0889*0 0.5333+0.0111*0])
plot(percent_dATP(1),LVDP_healthy{2}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*1 0.1333+0.0889*1 0.5333+0.0111*1])
plot(percent_dATP(2),LVDP_healthy{3}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*2 0.1333+0.0889*2 0.5333+0.0111*2])
plot(percent_dATP(3),LVDP_healthy{4}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*3 0.1333+0.0889*3 0.5333+0.0111*3])
plot(percent_dATP(4),LVDP_healthy{5}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*4 0.1333+0.0889*4 0.5333+0.0111*4])
plot(percent_dATP(5),LVDP_healthy{6}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*5 0.1333+0.0889*5 0.5333+0.0111*5])
plot(percent_dATP(6),LVDP_healthy{7}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*6 0.1333+0.0889*6 0.5333+0.0111*6])

plot(0,LVDP_failure{1}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*0 0.2667+0.0888*0 0.6000-0.0222*0])
plot(percent_dATP(1),LVDP_failure{2}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*1 0.2667+0.0888*1 0.6000-0.0222*1])
plot(percent_dATP(2),LVDP_failure{3}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*2 0.2667+0.0888*2 0.6000-0.0222*2])
plot(percent_dATP(3),LVDP_failure{4}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*3 0.2667+0.0888*3 0.6000-0.0222*3])
plot(percent_dATP(4),LVDP_failure{5}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*4 0.2667+0.0888*4 0.6000-0.0222*4])
plot(percent_dATP(5),LVDP_failure{6}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*5 0.2667+0.0888*5 0.6000-0.0222*5])
plot(percent_dATP(6),LVDP_failure{7}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*6 0.2667+0.0888*6 0.6000-0.0222*6])
xlabel('Percent dATP')
ylabel('LVDP (mmHg)')
set(gca, 'fontsize', 14)

subplot(4,4,6)
hold on
plot(0,plus_dPdt_F_healthy{1}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*0 0.1333+0.0889*0 0.5333+0.0111*0])
plot(percent_dATP(1),plus_dPdt_F_healthy{2}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*1 0.1333+0.0889*1 0.5333+0.0111*1])
plot(percent_dATP(2),plus_dPdt_F_healthy{3}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*2 0.1333+0.0889*2 0.5333+0.0111*2])
plot(percent_dATP(3),plus_dPdt_F_healthy{4}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*3 0.1333+0.0889*3 0.5333+0.0111*3])
plot(percent_dATP(4),plus_dPdt_F_healthy{5}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*4 0.1333+0.0889*4 0.5333+0.0111*4])
plot(percent_dATP(5),plus_dPdt_F_healthy{6}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*5 0.1333+0.0889*5 0.5333+0.0111*5])
plot(percent_dATP(6),plus_dPdt_F_healthy{7}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*6 0.1333+0.0889*6 0.5333+0.0111*6])

plot(0,plus_dPdt_F_failure{1}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*0 0.2667+0.0888*0 0.6000-0.0222*0])
plot(percent_dATP(1),plus_dPdt_F_failure{2}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*1 0.2667+0.0888*1 0.6000-0.0222*1])
plot(percent_dATP(2),plus_dPdt_F_failure{3}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*2 0.2667+0.0888*2 0.6000-0.0222*2])
plot(percent_dATP(3),plus_dPdt_F_failure{4}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*3 0.2667+0.0888*3 0.6000-0.0222*3])
plot(percent_dATP(4),plus_dPdt_F_failure{5}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*4 0.2667+0.0888*4 0.6000-0.0222*4])
plot(percent_dATP(5),plus_dPdt_F_failure{6}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*5 0.2667+0.0888*5 0.6000-0.0222*5])
plot(percent_dATP(6),plus_dPdt_F_failure{7}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*6 0.2667+0.0888*6 0.6000-0.0222*6])
xlabel('Percent dATP')
ylabel('+dP/dt (mmHg/s)')
set(gca, 'fontsize', 14)

subplot(4,4,7)
hold on
plot(0,minus_dPdt_F_healthy{1}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*0 0.1333+0.0889*0 0.5333+0.0111*0])
plot(percent_dATP(1),minus_dPdt_F_healthy{2}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*1 0.1333+0.0889*1 0.5333+0.0111*1])
plot(percent_dATP(2),minus_dPdt_F_healthy{3}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*2 0.1333+0.0889*2 0.5333+0.0111*2])
plot(percent_dATP(3),minus_dPdt_F_healthy{4}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*3 0.1333+0.0889*3 0.5333+0.0111*3])
plot(percent_dATP(4),minus_dPdt_F_healthy{5}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*4 0.1333+0.0889*4 0.5333+0.0111*4])
plot(percent_dATP(5),minus_dPdt_F_healthy{6}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*5 0.1333+0.0889*5 0.5333+0.0111*5])
plot(percent_dATP(6),minus_dPdt_F_healthy{7}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*6 0.1333+0.0889*6 0.5333+0.0111*6])

plot(0,minus_dPdt_F_failure{1}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*0 0.2667+0.0888*0 0.6000-0.0222*0])
plot(percent_dATP(1),minus_dPdt_F_failure{2}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*1 0.2667+0.0888*1 0.6000-0.0222*1])
plot(percent_dATP(2),minus_dPdt_F_failure{3}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*2 0.2667+0.0888*2 0.6000-0.0222*2])
plot(percent_dATP(3),minus_dPdt_F_failure{4}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*3 0.2667+0.0888*3 0.6000-0.0222*3])
plot(percent_dATP(4),minus_dPdt_F_failure{5}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*4 0.2667+0.0888*4 0.6000-0.0222*4])
plot(percent_dATP(5),minus_dPdt_F_failure{6}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*5 0.2667+0.0888*5 0.6000-0.0222*5])
plot(percent_dATP(6),minus_dPdt_F_failure{7}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*6 0.2667+0.0888*6 0.6000-0.0222*6])
xlabel('Percent dATP')
ylabel('-dP/dt (mmHg/s)')
set(gca, 'fontsize', 14)

subplot(4,4,8)
hold on
plot(0,work_rate_F_healthy{1}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*0 0.1333+0.0889*0 0.5333+0.0111*0])
plot(percent_dATP(1),work_rate_F_healthy{2}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*1 0.1333+0.0889*1 0.5333+0.0111*1])
plot(percent_dATP(2),work_rate_F_healthy{3}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*2 0.1333+0.0889*2 0.5333+0.0111*2])
plot(percent_dATP(3),work_rate_F_healthy{4}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*3 0.1333+0.0889*3 0.5333+0.0111*3])
plot(percent_dATP(4),work_rate_F_healthy{5}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*4 0.1333+0.0889*4 0.5333+0.0111*4])
plot(percent_dATP(5),work_rate_F_healthy{6}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*5 0.1333+0.0889*5 0.5333+0.0111*5])
plot(percent_dATP(6),work_rate_F_healthy{7}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*6 0.1333+0.0889*6 0.5333+0.0111*6])

plot(0,work_rate_F_failure{1}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*0 0.2667+0.0888*0 0.6000-0.0222*0])
plot(percent_dATP(1),work_rate_F_failure{2}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*1 0.2667+0.0888*1 0.6000-0.0222*1])
plot(percent_dATP(2),work_rate_F_failure{3}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*2 0.2667+0.0888*2 0.6000-0.0222*2])
plot(percent_dATP(3),work_rate_F_failure{4}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*3 0.2667+0.0888*3 0.6000-0.0222*3])
plot(percent_dATP(4),work_rate_F_failure{5}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*4 0.2667+0.0888*4 0.6000-0.0222*4])
plot(percent_dATP(5),work_rate_F_failure{6}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*5 0.2667+0.0888*5 0.6000-0.0222*5])
plot(percent_dATP(6),work_rate_F_failure{7}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*6 0.2667+0.0888*6 0.6000-0.0222*6])
xlabel('Percent dATP')
ylabel('Work Rate (mmHg*mL/s)')
set(gca, 'fontsize', 14)

subplot(4,4,9)
hold on
plot(0,ATP_F_healthy{1}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*0 0.1333+0.0889*0 0.5333+0.0111*0])
plot(percent_dATP(1),ATP_F_healthy{2}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*1 0.1333+0.0889*1 0.5333+0.0111*1])
plot(percent_dATP(2),ATP_F_healthy{3}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*2 0.1333+0.0889*2 0.5333+0.0111*2])
plot(percent_dATP(3),ATP_F_healthy{4}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*3 0.1333+0.0889*3 0.5333+0.0111*3])
plot(percent_dATP(4),ATP_F_healthy{5}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*4 0.1333+0.0889*4 0.5333+0.0111*4])
plot(percent_dATP(5),ATP_F_healthy{6}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*5 0.1333+0.0889*5 0.5333+0.0111*5])
plot(percent_dATP(6),ATP_F_healthy{7}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*6 0.1333+0.0889*6 0.5333+0.0111*6])

plot(0,ATP_F_failure{1}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*0 0.2667+0.0888*0 0.6000-0.0222*0])
plot(percent_dATP(1),ATP_F_failure{2}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*1 0.2667+0.0888*1 0.6000-0.0222*1])
plot(percent_dATP(2),ATP_F_failure{3}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*2 0.2667+0.0888*2 0.6000-0.0222*2])
plot(percent_dATP(3),ATP_F_failure{4}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*3 0.2667+0.0888*3 0.6000-0.0222*3])
plot(percent_dATP(4),ATP_F_failure{5}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*4 0.2667+0.0888*4 0.6000-0.0222*4])
plot(percent_dATP(5),ATP_F_failure{6}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*5 0.2667+0.0888*5 0.6000-0.0222*5])
plot(percent_dATP(6),ATP_F_failure{7}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*6 0.2667+0.0888*6 0.6000-0.0222*6])
xlabel('Percent dATP')
ylabel('ATP (mM)')
set(gca, 'fontsize', 14)

subplot(4,4,10)
hold on
plot(0,ADP_F_healthy{1}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*0 0.1333+0.0889*0 0.5333+0.0111*0])
plot(percent_dATP(1),ADP_F_healthy{2}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*1 0.1333+0.0889*1 0.5333+0.0111*1])
plot(percent_dATP(2),ADP_F_healthy{3}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*2 0.1333+0.0889*2 0.5333+0.0111*2])
plot(percent_dATP(3),ADP_F_healthy{4}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*3 0.1333+0.0889*3 0.5333+0.0111*3])
plot(percent_dATP(4),ADP_F_healthy{5}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*4 0.1333+0.0889*4 0.5333+0.0111*4])
plot(percent_dATP(5),ADP_F_healthy{6}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*5 0.1333+0.0889*5 0.5333+0.0111*5])
plot(percent_dATP(6),ADP_F_healthy{7}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*6 0.1333+0.0889*6 0.5333+0.0111*6])

plot(0,ADP_F_failure{1}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*0 0.2667+0.0888*0 0.6000-0.0222*0])
plot(percent_dATP(1),ADP_F_failure{2}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*1 0.2667+0.0888*1 0.6000-0.0222*1])
plot(percent_dATP(2),ADP_F_failure{3}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*2 0.2667+0.0888*2 0.6000-0.0222*2])
plot(percent_dATP(3),ADP_F_failure{4}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*3 0.2667+0.0888*3 0.6000-0.0222*3])
plot(percent_dATP(4),ADP_F_failure{5}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*4 0.2667+0.0888*4 0.6000-0.0222*4])
plot(percent_dATP(5),ADP_F_failure{6}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*5 0.2667+0.0888*5 0.6000-0.0222*5])
plot(percent_dATP(6),ADP_F_failure{7}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*6 0.2667+0.0888*6 0.6000-0.0222*6])
xlabel('Percent dATP')
ylabel('ADP (\muM)')
set(gca, 'fontsize', 14)

subplot(4,4,11)
hold on
plot(0,Pi_F_healthy{1}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*0 0.1333+0.0889*0 0.5333+0.0111*0])
plot(percent_dATP(1),Pi_F_healthy{2}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*1 0.1333+0.0889*1 0.5333+0.0111*1])
plot(percent_dATP(2),Pi_F_healthy{3}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*2 0.1333+0.0889*2 0.5333+0.0111*2])
plot(percent_dATP(3),Pi_F_healthy{4}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*3 0.1333+0.0889*3 0.5333+0.0111*3])
plot(percent_dATP(4),Pi_F_healthy{5}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*4 0.1333+0.0889*4 0.5333+0.0111*4])
plot(percent_dATP(5),Pi_F_healthy{6}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*5 0.1333+0.0889*5 0.5333+0.0111*5])
plot(percent_dATP(6),Pi_F_healthy{7}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*6 0.1333+0.0889*6 0.5333+0.0111*6])

plot(0,Pi_F_failure{1}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*0 0.2667+0.0888*0 0.6000-0.0222*0])
plot(percent_dATP(1),Pi_F_failure{2}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*1 0.2667+0.0888*1 0.6000-0.0222*1])
plot(percent_dATP(2),Pi_F_failure{3}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*2 0.2667+0.0888*2 0.6000-0.0222*2])
plot(percent_dATP(3),Pi_F_failure{4}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*3 0.2667+0.0888*3 0.6000-0.0222*3])
plot(percent_dATP(4),Pi_F_failure{5}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*4 0.2667+0.0888*4 0.6000-0.0222*4])
plot(percent_dATP(5),Pi_F_failure{6}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*5 0.2667+0.0888*5 0.6000-0.0222*5])
plot(percent_dATP(6),Pi_F_failure{7}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*6 0.2667+0.0888*6 0.6000-0.0222*6])
xlabel('Percent dATP')
ylabel('Pi (mM)')
set(gca, 'fontsize', 14)

subplot(4,4,12)
hold on
plot(0,MVO2_F_healthy{1}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*0 0.1333+0.0889*0 0.5333+0.0111*0])
plot(percent_dATP(1),MVO2_F_healthy{2}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*1 0.1333+0.0889*1 0.5333+0.0111*1])
plot(percent_dATP(2),MVO2_F_healthy{3}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*2 0.1333+0.0889*2 0.5333+0.0111*2])
plot(percent_dATP(3),MVO2_F_healthy{4}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*3 0.1333+0.0889*3 0.5333+0.0111*3])
plot(percent_dATP(4),MVO2_F_healthy{5}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*4 0.1333+0.0889*4 0.5333+0.0111*4])
plot(percent_dATP(5),MVO2_F_healthy{6}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*5 0.1333+0.0889*5 0.5333+0.0111*5])
plot(percent_dATP(6),MVO2_F_healthy{7}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*6 0.1333+0.0889*6 0.5333+0.0111*6])

plot(0,MVO2_F_failure{1}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*0 0.2667+0.0888*0 0.6000-0.0222*0])
plot(percent_dATP(1),MVO2_F_failure{2}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*1 0.2667+0.0888*1 0.6000-0.0222*1])
plot(percent_dATP(2),MVO2_F_failure{3}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*2 0.2667+0.0888*2 0.6000-0.0222*2])
plot(percent_dATP(3),MVO2_F_failure{4}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*3 0.2667+0.0888*3 0.6000-0.0222*3])
plot(percent_dATP(4),MVO2_F_failure{5}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*4 0.2667+0.0888*4 0.6000-0.0222*4])
plot(percent_dATP(5),MVO2_F_failure{6}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*5 0.2667+0.0888*5 0.6000-0.0222*5])
plot(percent_dATP(6),MVO2_F_failure{7}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*6 0.2667+0.0888*6 0.6000-0.0222*6])
xlabel('Percent dATP')
ylabel('MVO_{2} (uM O_{2}/min/g tissue)')
set(gca, 'fontsize', 14)

subplot(4,4,13)
hold on
plot(0,PCrATP_F_healthy{1}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*0 0.1333+0.0889*0 0.5333+0.0111*0])
plot(percent_dATP(1),PCrATP_F_healthy{2}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*1 0.1333+0.0889*1 0.5333+0.0111*1])
plot(percent_dATP(2),PCrATP_F_healthy{3}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*2 0.1333+0.0889*2 0.5333+0.0111*2])
plot(percent_dATP(3),PCrATP_F_healthy{4}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*3 0.1333+0.0889*3 0.5333+0.0111*3])
plot(percent_dATP(4),PCrATP_F_healthy{5}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*4 0.1333+0.0889*4 0.5333+0.0111*4])
plot(percent_dATP(5),PCrATP_F_healthy{6}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*5 0.1333+0.0889*5 0.5333+0.0111*5])
plot(percent_dATP(6),PCrATP_F_healthy{7}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*6 0.1333+0.0889*6 0.5333+0.0111*6])

plot(0,PCrATP_F_failure{1}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*0 0.2667+0.0888*0 0.6000-0.0222*0])
plot(percent_dATP(1),PCrATP_F_failure{2}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*1 0.2667+0.0888*1 0.6000-0.0222*1])
plot(percent_dATP(2),PCrATP_F_failure{3}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*2 0.2667+0.0888*2 0.6000-0.0222*2])
plot(percent_dATP(3),PCrATP_F_failure{4}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*3 0.2667+0.0888*3 0.6000-0.0222*3])
plot(percent_dATP(4),PCrATP_F_failure{5}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*4 0.2667+0.0888*4 0.6000-0.0222*4])
plot(percent_dATP(5),PCrATP_F_failure{6}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*5 0.2667+0.0888*5 0.6000-0.0222*5])
plot(percent_dATP(6),PCrATP_F_failure{7}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*6 0.2667+0.0888*6 0.6000-0.0222*6])
xlabel('Percent dATP')
ylabel('CrP/ATP ratio')
set(gca, 'fontsize', 14)

subplot(4,4,14)
hold on
plot(0,XB_turnover_F_healthy{1}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*0 0.1333+0.0889*0 0.5333+0.0111*0])
plot(percent_dATP(1),XB_turnover_F_healthy{2}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*1 0.1333+0.0889*1 0.5333+0.0111*1])
plot(percent_dATP(2),XB_turnover_F_healthy{3}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*2 0.1333+0.0889*2 0.5333+0.0111*2])
plot(percent_dATP(3),XB_turnover_F_healthy{4}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*3 0.1333+0.0889*3 0.5333+0.0111*3])
plot(percent_dATP(4),XB_turnover_F_healthy{5}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*4 0.1333+0.0889*4 0.5333+0.0111*4])
plot(percent_dATP(5),XB_turnover_F_healthy{6}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*5 0.1333+0.0889*5 0.5333+0.0111*5])
plot(percent_dATP(6),XB_turnover_F_healthy{7}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*6 0.1333+0.0889*6 0.5333+0.0111*6])

plot(0,XB_turnover_F_failure{1}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*0 0.2667+0.0888*0 0.6000-0.0222*0])
plot(percent_dATP(1),XB_turnover_F_failure{2}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*1 0.2667+0.0888*1 0.6000-0.0222*1])
plot(percent_dATP(2),XB_turnover_F_failure{3}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*2 0.2667+0.0888*2 0.6000-0.0222*2])
plot(percent_dATP(3),XB_turnover_F_failure{4}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*3 0.2667+0.0888*3 0.6000-0.0222*3])
plot(percent_dATP(4),XB_turnover_F_failure{5}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*4 0.2667+0.0888*4 0.6000-0.0222*4])
plot(percent_dATP(5),XB_turnover_F_failure{6}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*5 0.2667+0.0888*5 0.6000-0.0222*5])
plot(percent_dATP(6),XB_turnover_F_failure{7}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*6 0.2667+0.0888*6 0.6000-0.0222*6])
xlabel('Percent dATP')
ylabel('XB Cycling (1/s)')
set(gca, 'fontsize', 14)

subplot(4,4,15)
hold on
plot(0,ATPase_F_healthy{1}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*0 0.1333+0.0889*0 0.5333+0.0111*0])
plot(percent_dATP(1),ATPase_F_healthy{2}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*1 0.1333+0.0889*1 0.5333+0.0111*1])
plot(percent_dATP(2),ATPase_F_healthy{3}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*2 0.1333+0.0889*2 0.5333+0.0111*2])
plot(percent_dATP(3),ATPase_F_healthy{4}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*3 0.1333+0.0889*3 0.5333+0.0111*3])
plot(percent_dATP(4),ATPase_F_healthy{5}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*4 0.1333+0.0889*4 0.5333+0.0111*4])
plot(percent_dATP(5),ATPase_F_healthy{6}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*5 0.1333+0.0889*5 0.5333+0.0111*5])
plot(percent_dATP(6),ATPase_F_healthy{7}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*6 0.1333+0.0889*6 0.5333+0.0111*6])

plot(0,ATPase_F_failure{1}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*0 0.2667+0.0888*0 0.6000-0.0222*0])
plot(percent_dATP(1),ATPase_F_failure{2}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*1 0.2667+0.0888*1 0.6000-0.0222*1])
plot(percent_dATP(2),ATPase_F_failure{3}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*2 0.2667+0.0888*2 0.6000-0.0222*2])
plot(percent_dATP(3),ATPase_F_failure{4}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*3 0.2667+0.0888*3 0.6000-0.0222*3])
plot(percent_dATP(4),ATPase_F_failure{5}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*4 0.2667+0.0888*4 0.6000-0.0222*4])
plot(percent_dATP(5),ATPase_F_failure{6}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*5 0.2667+0.0888*5 0.6000-0.0222*5])
plot(percent_dATP(6),ATPase_F_failure{7}{1},'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*6 0.2667+0.0888*6 0.6000-0.0222*6])
xlabel('Percent dATP')
ylabel('ATPase (M/s/L)')
set(gca, 'fontsize', 14)

subplot(4,4,16)
hold on
plot(0,efficiency_F_healthy{1}{1}./1000,'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*0 0.1333+0.0889*0 0.5333+0.0111*0])
plot(percent_dATP(1),efficiency_F_healthy{2}{1}./1000,'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*1 0.1333+0.0889*1 0.5333+0.0111*1])
plot(percent_dATP(2),efficiency_F_healthy{3}{1}./1000,'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*2 0.1333+0.0889*2 0.5333+0.0111*2])
plot(percent_dATP(3),efficiency_F_healthy{4}{1}./1000,'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*3 0.1333+0.0889*3 0.5333+0.0111*3])
plot(percent_dATP(4),efficiency_F_healthy{5}{1}./1000,'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*4 0.1333+0.0889*4 0.5333+0.0111*4])
plot(percent_dATP(5),efficiency_F_healthy{6}{1}./1000,'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*5 0.1333+0.0889*5 0.5333+0.0111*5])
plot(percent_dATP(6),efficiency_F_healthy{7}{1}./1000,'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*6 0.1333+0.0889*6 0.5333+0.0111*6])

plot(0,efficiency_F_failure{1}{1}./1000,'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*0 0.2667+0.0888*0 0.6000-0.0222*0])
plot(percent_dATP(1),efficiency_F_failure{2}{1}./1000,'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*1 0.2667+0.0888*1 0.6000-0.0222*1])
plot(percent_dATP(2),efficiency_F_failure{3}{1}./1000,'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*2 0.2667+0.0888*2 0.6000-0.0222*2])
plot(percent_dATP(3),efficiency_F_failure{4}{1}./1000,'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*3 0.2667+0.0888*3 0.6000-0.0222*3])
plot(percent_dATP(4),efficiency_F_failure{5}{1}./1000,'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*4 0.2667+0.0888*4 0.6000-0.0222*4])
plot(percent_dATP(5),efficiency_F_failure{6}{1}./1000,'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*5 0.2667+0.0888*5 0.6000-0.0222*5])
plot(percent_dATP(6),efficiency_F_failure{7}{1}./1000,'o','markersize',10,'linewidth',3,'color',[0.6667+0.0299*6 0.2667+0.0888*6 0.6000-0.0222*6])
xlabel('Percent dATP')
ylabel('Efficiency (mL^{2}*mmHg*ms/M)')
set(gca, 'fontsize', 14)


%% Figure S8
dATP_sim = [0 1 2 3 4 5 6 7 8 9 10 15 20];
for i = 1:length(dATP_sim)
    dATP = dATP_sim(i);
    [V_LV_store_F, P_LV_store_F, max_force, FS, CO, EF, LVDP, plus_dPdt_F, minus_dPdt_F, work_rate_F, ATP_F, ADP_F, Pi_F, MVO2_F, PCrATP_F, XB_turnover_F, ATPase_F, efficiency_F, SRX_store_F] = cardiovascular_model(1, 0, 0, dATP, 538, 540, 5.7, 2, 80, 4, 25, 350, 280, 800, 4);
    SRX_store_F_model(i) = SRX_store_F;
end

for i = 1:length(dATP_sim)
    SRX_store_F_model_mat = cell2mat(SRX_store_F_model(i));
    SRX_max(i) = max(SRX_store_F_model_mat);
end

figure
plot(dATP_sim, 1-SRX_max,'-o','linewidth',3,'color',[0 0 0],'markersize',6)
xlabel('Percent dATP')
ylabel('Fraction SRX')
set(gca, 'fontsize', 14)

%% Figure S9
% ATP
[V_LV_store_F_ATP, P_LV_store_F_ATP, max_force_ATP, FS_ATP, CO_ATP, EF_ATP, LVDP_ATP, plus_dPdt_F_ATP, minus_dPdt_F_ATP, work_rate_F_ATP, ATP_F_ATP, ADP_F_ATP, Pi_F_ATP, MVO2_F_ATP, PCrATP_F_ATP, XB_turnover_F_ATP, ATPase_F_ATP, efficiency_F_ATP, SRX_store_F_ATP] = cardiovascular_model(0, 0, 0, 0, 250, 304.7, 4, 2, 80, 4, 25, 0.4, 280, 800, 4);

% ka
[V_LV_store_F_ka, P_LV_store_F_ka, max_force_ka, FS_ka, CO_ka, EF_ka, LVDP_ka, plus_dPdt_F_ka, minus_dPdt_F_ka, work_rate_F_ka, ATP_F_ka, ADP_F_ka, Pi_F_ka, MVO2_F_ka, PCrATP_F_ka, XB_turnover_F_ka, ATPase_F_ka, efficiency_F_ka, SRX_store_F_ka] = cardiovascular_model(0, 0, 0, 2, 538, 304.7, 4, 2, 80, 4, 25, 0.4, 280, 800, 4);

% XB
[V_LV_store_F_XB, P_LV_store_F_XB, max_force_XB, FS_XB, CO_XB, EF_XB, LVDP_XB, plus_dPdt_F_XB, minus_dPdt_F_XB, work_rate_F_XB, ATP_F_XB, ADP_F_XB, Pi_F_XB, MVO2_F_XB, PCrATP_F_XB, XB_turnover_F_XB, ATPase_F_XB, efficiency_F_XB, SRX_store_F_XB] = cardiovascular_model(0, 0, 0, 2, 538, 540, 5.7, 2, 80, 4, 25, 0.4, 280, 800, 4);

% SRX
[V_LV_store_F_SRX, P_LV_store_F_SRX, max_force_SRX, FS_SRX, CO_SRX, EF_SRX, LVDP_SRX, plus_dPdt_F_SRX, minus_dPdt_F_SRX, work_rate_F_SRX, ATP_F_SRX, ADP_F_SRX, Pi_F_SRX, MVO2_F_SRX, PCrATP_F_SRX, XB_turnover_F_SRX, ATPase_F_SRX, efficiency_F_SRX, SRX_store_F_SRX] = cardiovascular_model(0, 0, 0, 2, 250, 304.7, 4, 2, 80, 4, 25, 350, 280, 800, 4);

% Ca
[V_LV_store_F_Ca, P_LV_store_F_Ca, max_force_Ca, FS_Ca, CO_Ca, EF_Ca, LVDP_Ca, plus_dPdt_F_Ca, minus_dPdt_F_Ca, work_rate_F_Ca, ATP_F_Ca, ADP_F_Ca, Pi_F_Ca, MVO2_F_Ca, PCrATP_F_Ca, XB_turnover_F_Ca, ATPase_F_Ca, efficiency_F_Ca, SRX_store_F_Ca] = cardiovascular_model(1, 0, 0, 2, 250, 304.7, 4, 2, 80, 4, 25, 0.4, 280, 800, 4);

% XB + SRX
[V_LV_store_F_XB_SRX, P_LV_store_F_XB_SRX, max_force_XB_SRX, FS_XB_SRX, CO_XB_SRX, EF_XB_SRX, LVDP_XB_SRX, plus_dPdt_F_XB_SRX, minus_dPdt_F_XB_SRX, work_rate_F_XB_SRX, ATP_F_XB_SRX, ADP_F_XB_SRX, Pi_F_XB_SRX, MVO2_F_XB_SRX, PCrATP_F_XB_SRX, XB_turnover_F_XB_SRX, ATPase_F_XB_SRX, efficiency_F_XB_SRX, SRX_store_F_XB_SRX] = cardiovascular_model(0, 0, 0, 2, 538, 540, 5.7, 2, 80, 4, 25, 350, 280, 800, 4);

% XB + Ca
[V_LV_store_F_XB_Ca, P_LV_store_F_XB_Ca, max_force_XB_Ca, FS_XB_Ca, CO_XB_Ca, EF_XB_Ca, LVDP_XB_Ca, plus_dPdt_F_XB_Ca, minus_dPdt_F_XB_Ca, work_rate_F_XB_Ca, ATP_F_XB_Ca, ADP_F_XB_Ca, Pi_F_XB_Ca, MVO2_F_XB_Ca, PCrATP_F_XB_Ca, XB_turnover_F_XB_Ca, ATPase_F_XB_Ca, efficiency_F_XB_Ca, SRX_store_F_XB_Ca] = cardiovascular_model(1, 0, 0, 2, 538, 540, 5.7, 2, 80, 4, 25, 0.4, 280, 800, 4);

% Ca + SRX
[V_LV_store_F_Ca_SRX, P_LV_store_F_Ca_SRX, max_force_Ca_SRX, FS_Ca_SRX, CO_Ca_SRX, EF_Ca_SRX, LVDP_Ca_SRX, plus_dPdt_F_Ca_SRX, minus_dPdt_F_Ca_SRX, work_rate_F_Ca_SRX, ATP_F_Ca_SRX, ADP_F_Ca_SRX, Pi_F_Ca_SRX, MVO2_F_Ca_SRX, PCrATP_F_Ca_SRX, XB_turnover_F_Ca_SRX, ATPase_F_Ca_SRX, efficiency_F_Ca_SRX, SRX_store_F_Ca_SRX] = cardiovascular_model(1, 0, 0, 2, 250, 304.7, 4, 2, 80, 4, 25, 350, 280, 800, 4);

% XB + Ca + SRX
[V_LV_store_F_XB_Ca_SRX, P_LV_store_F_XB_Ca_SRX, max_force_XB_Ca_SRX, FS_XB_Ca_SRX, CO_XB_Ca_SRX, EF_XB_Ca_SRX, LVDP_XB_Ca_SRX, plus_dPdt_F_XB_Ca_SRX, minus_dPdt_F_XB_Ca_SRX, work_rate_F_XB_Ca_SRX, ATP_F_XB_Ca_SRX, ADP_F_XB_Ca_SRX, Pi_F_XB_Ca_SRX, MVO2_F_XB_Ca_SRX, PCrATP_F_XB_Ca_SRX, XB_turnover_F_XB_Ca_SRX, ATPase_F_XB_Ca_SRX, efficiency_F_XB_Ca_SRX, SRX_store_F_XB_Ca_SRX] = cardiovascular_model(1, 0, 0, 2, 538, 540, 5.7, 2, 80, 4, 25, 350, 280, 800, 4);

% From Nowakowski et al. 2013
data_ATP = [0	70.65	30	104.20	3113	-2730	0	1.76	0	0.00	0];
data_dATP = [0	80.68	41	138.20	4197	-3410	0	1.60	0	0.00	0];
data = (data_dATP - data_ATP)./data_ATP;

figure
subplot(5,2,1)
hold on
plot(0, (max_force_ka{1}-max_force_ATP{1})./max_force_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(1, (max_force_XB{1}-max_force_ATP{1})./max_force_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(2, (max_force_SRX{1}-max_force_ATP{1})./max_force_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(3, (max_force_Ca{1}-max_force_ATP{1})./max_force_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(4, (max_force_XB_SRX{1}-max_force_ATP{1})./max_force_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(5, (max_force_XB_Ca{1}-max_force_ATP{1})./max_force_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(6, (max_force_Ca_SRX{1}-max_force_ATP{1})./max_force_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(7, (max_force_XB_Ca_SRX{1}-max_force_ATP{1})./max_force_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
yline((max_force_ATP{1}-max_force_ATP{1})./max_force_ATP{1},':','color',[0.2, 0.1333, 0.5333],'linewidth',3);
ylabel('Maximum Force')
set(gca,'FontSize',14,'XTick',[],'XTickLabel',{'','','','','','','',''})
xlim([0 7])

subplot(5,2,2)
hold on
plot(0, (EF_ka{1}-EF_ATP{1})./EF_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(1, (EF_XB{1}-EF_ATP{1})./EF_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(2, (EF_SRX{1}-EF_ATP{1})./EF_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(3, (EF_Ca{1}-EF_ATP{1})./EF_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(4, (EF_XB_SRX{1}-EF_ATP{1})./EF_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(5, (EF_XB_Ca{1}-EF_ATP{1})./EF_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(6, (EF_Ca_SRX{1}-EF_ATP{1})./EF_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(7, (EF_XB_Ca_SRX{1}-EF_ATP{1})./EF_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
yline((EF_ATP{1}-EF_ATP{1})./EF_ATP{1},':','color',[0.2, 0.1333, 0.5333],'linewidth',3);
yline(data(2),':','color',[0.2667, 0.6667, 0.6],'linewidth',3);
ylabel('EF')
set(gca,'FontSize',14,'XTick',[],'XTickLabel',{'','','','','','','',''})
xlim([0 7])

subplot(5,2,3)
hold on
plot(0, (CO_ka{1}-CO_ATP{1})./CO_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(1, (CO_XB{1}-CO_ATP{1})./CO_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(2, (CO_SRX{1}-CO_ATP{1})./CO_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(3, (CO_Ca{1}-CO_ATP{1})./CO_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(4, (CO_XB_SRX{1}-CO_ATP{1})./CO_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(5, (CO_XB_Ca{1}-CO_ATP{1})./CO_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(6, (CO_Ca_SRX{1}-CO_ATP{1})./CO_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(7, (CO_XB_Ca_SRX{1}-CO_ATP{1})./CO_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
yline((CO_ATP{1}-CO_ATP{1})./CO_ATP{1},':','color',[0.2, 0.1333, 0.5333],'linewidth',3);
yline(data(3),':','color',[0.2667, 0.6667, 0.6],'linewidth',3);
ylabel('CO')
set(gca,'FontSize',14,'XTick',[],'XTickLabel',{'','','','','','','',''})
xlim([0 7])

subplot(5,2,4)
hold on
plot(0, (LVDP_ka{1}-LVDP_ATP{1})./LVDP_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(1, (LVDP_XB{1}-LVDP_ATP{1})./LVDP_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(2, (LVDP_SRX{1}-LVDP_ATP{1})./LVDP_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(3, (LVDP_Ca{1}-LVDP_ATP{1})./LVDP_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(4, (LVDP_XB_SRX{1}-LVDP_ATP{1})./LVDP_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(5, (LVDP_XB_Ca{1}-LVDP_ATP{1})./LVDP_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(6, (LVDP_Ca_SRX{1}-LVDP_ATP{1})./LVDP_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(7, (LVDP_XB_Ca_SRX{1}-LVDP_ATP{1})./LVDP_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
yline((LVDP_ATP{1}-LVDP_ATP{1})./LVDP_ATP{1},':','color',[0.2, 0.1333, 0.5333],'linewidth',3);
yline(data(4),':','color',[0.2667, 0.6667, 0.6],'linewidth',3);
ylabel('LVDP')
set(gca,'FontSize',14,'XTick',[],'XTickLabel',{'','','','','','','',''})
xlim([0 7])

subplot(5,2,5)
hold on
plot(0, (plus_dPdt_F_ka{1}-plus_dPdt_F_ATP{1})./plus_dPdt_F_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(1, (plus_dPdt_F_XB{1}-plus_dPdt_F_ATP{1})./plus_dPdt_F_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(2, (plus_dPdt_F_SRX{1}-plus_dPdt_F_ATP{1})./plus_dPdt_F_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(3, (plus_dPdt_F_Ca{1}-plus_dPdt_F_ATP{1})./plus_dPdt_F_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(4, (plus_dPdt_F_XB_SRX{1}-plus_dPdt_F_ATP{1})./plus_dPdt_F_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(5, (plus_dPdt_F_XB_Ca{1}-plus_dPdt_F_ATP{1})./plus_dPdt_F_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(6, (plus_dPdt_F_Ca_SRX{1}-plus_dPdt_F_ATP{1})./plus_dPdt_F_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(7, (plus_dPdt_F_XB_Ca_SRX{1}-plus_dPdt_F_ATP{1})./plus_dPdt_F_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
yline((plus_dPdt_F_ATP{1}-plus_dPdt_F_ATP{1})./plus_dPdt_F_ATP{1},':','color',[0.2, 0.1333, 0.5333],'linewidth',3);
yline(data(5),':','color',[0.2667, 0.6667, 0.6],'linewidth',3);
ylabel('+dP/dt')
set(gca,'FontSize',14,'XTick',[],'XTickLabel',{'','','','','','','',''})
xlim([0 7])

subplot(5,2,6)
hold on
plot(0, (minus_dPdt_F_ka{1}-minus_dPdt_F_ATP{1})./minus_dPdt_F_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(1, (minus_dPdt_F_XB{1}-minus_dPdt_F_ATP{1})./minus_dPdt_F_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(2, (minus_dPdt_F_SRX{1}-minus_dPdt_F_ATP{1})./minus_dPdt_F_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(3, (minus_dPdt_F_Ca{1}-minus_dPdt_F_ATP{1})./minus_dPdt_F_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(4, (minus_dPdt_F_XB_SRX{1}-minus_dPdt_F_ATP{1})./minus_dPdt_F_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(5, (minus_dPdt_F_XB_Ca{1}-minus_dPdt_F_ATP{1})./minus_dPdt_F_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(6, (minus_dPdt_F_Ca_SRX{1}-minus_dPdt_F_ATP{1})./minus_dPdt_F_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(7, (minus_dPdt_F_XB_Ca_SRX{1}-minus_dPdt_F_ATP{1})./minus_dPdt_F_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
yline((minus_dPdt_F_ATP{1}-minus_dPdt_F_ATP{1})./minus_dPdt_F_ATP{1},':','color',[0.2, 0.1333, 0.5333],'linewidth',3);
yline(data(6),':','color',[0.2667, 0.6667, 0.6],'linewidth',3);
ylabel('-dP/dt')
set(gca,'FontSize',14,'XTick',[],'XTickLabel',{'','','','','','','',''})
xlim([0 7])

subplot(5,2,7)
hold on
plot(0, (work_rate_F_ka{1}-work_rate_F_ATP{1})./work_rate_F_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(1, (work_rate_F_XB{1}-work_rate_F_ATP{1})./work_rate_F_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(2, (work_rate_F_SRX{1}-work_rate_F_ATP{1})./work_rate_F_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(3, (work_rate_F_Ca{1}-work_rate_F_ATP{1})./work_rate_F_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(4, (work_rate_F_XB_SRX{1}-work_rate_F_ATP{1})./work_rate_F_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(5, (work_rate_F_XB_Ca{1}-work_rate_F_ATP{1})./work_rate_F_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(6, (work_rate_F_Ca_SRX{1}-work_rate_F_ATP{1})./work_rate_F_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(7, (work_rate_F_XB_Ca_SRX{1}-work_rate_F_ATP{1})./work_rate_F_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
yline((work_rate_F_ATP{1}-work_rate_F_ATP{1})./work_rate_F_ATP{1},':','color',[0.2, 0.1333, 0.5333],'linewidth',3);
ylabel('Work Rate')
set(gca,'FontSize',14,'XTick',[],'XTickLabel',{'','','','','','','',''})
xlim([0 7])

subplot(5,2,8)
hold on
plot(0, (PCrATP_F_ka{1}-PCrATP_F_ATP{1})./PCrATP_F_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(1, (PCrATP_F_XB{1}-PCrATP_F_ATP{1})./PCrATP_F_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(2, (PCrATP_F_SRX{1}-PCrATP_F_ATP{1})./PCrATP_F_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(3, (PCrATP_F_Ca{1}-PCrATP_F_ATP{1})./PCrATP_F_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(4, (PCrATP_F_XB_SRX{1}-PCrATP_F_ATP{1})./PCrATP_F_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(5, (PCrATP_F_XB_Ca{1}-PCrATP_F_ATP{1})./PCrATP_F_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(6, (PCrATP_F_Ca_SRX{1}-PCrATP_F_ATP{1})./PCrATP_F_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(7, (PCrATP_F_XB_Ca_SRX{1}-PCrATP_F_ATP{1})./PCrATP_F_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
yline((PCrATP_F_ATP{1}-PCrATP_F_ATP{1})./PCrATP_F_ATP{1},':','color',[0.2, 0.1333, 0.5333],'linewidth',3);
yline(data(8),':','color',[0.2667, 0.6667, 0.6],'linewidth',3);
ylabel('CrP/ATP ratio')
set(gca,'FontSize',14,'XTick',[],'XTickLabel',{'','','','','','','',''})
xlim([0 7])

subplot(5,2,9)
hold on
plot(0, (MVO2_F_ka{1}-MVO2_F_ATP{1})./MVO2_F_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(1, (MVO2_F_XB{1}-MVO2_F_ATP{1})./MVO2_F_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(2, (MVO2_F_SRX{1}-MVO2_F_ATP{1})./MVO2_F_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(3, (MVO2_F_Ca{1}-MVO2_F_ATP{1})./MVO2_F_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(4, (MVO2_F_XB_SRX{1}-MVO2_F_ATP{1})./MVO2_F_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(5, (MVO2_F_XB_Ca{1}-MVO2_F_ATP{1})./MVO2_F_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(6, (MVO2_F_Ca_SRX{1}-MVO2_F_ATP{1})./MVO2_F_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(7, (MVO2_F_XB_Ca_SRX{1}-MVO2_F_ATP{1})./MVO2_F_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
yline((MVO2_F_ATP{1}-MVO2_F_ATP{1})./MVO2_F_ATP{1},':','color',[0.2, 0.1333, 0.5333],'linewidth',3);
ylabel('MVO_{2}')
set(gca,'FontSize',14,'XTickLabel',{'k_a','XB','SRX','Ca','XB + SRX','XB + Ca','Ca + SRX','XB + Ca + SRX'})
xlim([0 7])

subplot(5,2,10)
hold on
plot(0, (efficiency_F_ka{1}-efficiency_F_ATP{1})./efficiency_F_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(1, (efficiency_F_XB{1}-efficiency_F_ATP{1})./efficiency_F_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(2, (efficiency_F_SRX{1}-efficiency_F_ATP{1})./efficiency_F_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(3, (efficiency_F_Ca{1}-efficiency_F_ATP{1})./efficiency_F_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(4, (efficiency_F_XB_SRX{1}-efficiency_F_ATP{1})./efficiency_F_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(5, (efficiency_F_XB_Ca{1}-efficiency_F_ATP{1})./efficiency_F_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(6, (efficiency_F_Ca_SRX{1}-efficiency_F_ATP{1})./efficiency_F_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
plot(7, (efficiency_F_XB_Ca_SRX{1}-efficiency_F_ATP{1})./efficiency_F_ATP{1},'x','color',[0.5 0.5 0.5],'linewidth',3,'markersize',10)
yline((efficiency_F_ATP{1}-efficiency_F_ATP{1})./efficiency_F_ATP{1},':','color',[0.2, 0.1333, 0.5333],'linewidth',3);
ylabel('Efficiency')
set(gca,'FontSize',14,'XTickLabel',{'k_a','XB','SRX','Ca','XB + SRX','XB + Ca','Ca + SRX','XB + Ca + SRX'})
xlim([0 7])

%% Figure S10

% Simulations for this figure are run for Figure S9
figure
hold on
plot(V_LV_store_F_ATP{1},P_LV_store_F_ATP{1}-47,'linewidth',3,'color',[0.2+0.0111*0 0.1333+0.0889*0 0.5333+0.0111*0])
plot(V_LV_store_F_XB_SRX{1},P_LV_store_F_XB_SRX{1}-47,'--','linewidth',3,'color',[0.5,0.5,0.5])
plot(V_LV_store_F_Ca{1},P_LV_store_F_Ca{1}-47,'-.','linewidth',3,'color',[0.5,0.5,0.5])
plot(V_LV_store_F_XB_Ca_SRX{1},P_LV_store_F_XB_Ca_SRX{1}-47,'linewidth',3,'color',[0.2+0.0111*6 0.1333+0.0889*6 0.5333+0.0111*6])
xlabel('Volume (mL)')
ylabel('Pressure (mmHg)')
legend('ATP','dATP XB','dATP Ca','dATP XB + Ca')
set(gca, 'fontsize', 14)

%% Figure S11
dATP_sim = [0 2 4 6 8 10];
for i = 1:length(dATP_sim)
    dATP = dATP_sim(i);
    [V_LV_store_F, P_LV_store_F, max_force, FS, CO, EF, LVDP, plus_dPdt_F, minus_dPdt_F, work_rate_F, ATP_F, ADP_F, Pi_F, MVO2_F, PCrATP_F, XB_turnover_F, ATPase_F, efficiency_F, SRX_store_F] = cardiovascular_model(0, 0, 0, dATP, 538, 540, 5.7, 2, 80, 4, 25, 350, 280, 800, 4);
    V_LV_store_F_XB_SRX{i} = V_LV_store_F;
    P_LV_store_F_XB_SRX{i} = P_LV_store_F;
    max_force_XB_SRX{i} = max_force;
    FS_XB_SRX{i} = FS;
    CO_XB_SRX{i} = CO;
    EF_XB_SRX{i} = EF;
    LVDP_XB_SRX{i} = LVDP;
    plus_dPdt_F_XB_SRX{i} = plus_dPdt_F;
    minus_dPdt_F_XB_SRX{i} = minus_dPdt_F;
    work_rate_F_XB_SRX{i} = work_rate_F;
    ATP_F_XB_SRX{i} = ATP_F; 
    ADP_F_XB_SRX{i} = ADP_F;
    Pi_F_XB_SRX{i} = Pi_F;
    MVO2_F_XB_SRX{i} = MVO2_F;
    PCrATP_F_XB_SRX{i} = PCrATP_F;
    XB_turnover_F_XB_SRX{i} = XB_turnover_F;
    ATPase_F_XB_SRX{i} = ATPase_F;
    efficiency_F_XB_SRX{i} = efficiency_F;
end

% A
figure
hold on
plot(V_LV_store_F_XB_SRX{1}{1},P_LV_store_F_XB_SRX{1}{1}-47,':','linewidth',3,'color',[0.2+0.0111*0 0.1333+0.0889*0 0.5333+0.0111*0])
plot(V_LV_store_F_XB_SRX{2}{1},P_LV_store_F_XB_SRX{2}{1}-47,'linewidth',3,'color',[0.2+0.0111*1 0.1333+0.0889*1 0.5333+0.0111*1])
plot(V_LV_store_F_XB_SRX{3}{1},P_LV_store_F_XB_SRX{3}{1}-47,'linewidth',3,'color',[0.2+0.0111*2 0.1333+0.0889*2 0.5333+0.0111*2])
plot(V_LV_store_F_XB_SRX{4}{1},P_LV_store_F_XB_SRX{4}{1}-47,'linewidth',3,'color',[0.2+0.0111*3 0.1333+0.0889*3 0.5333+0.0111*3])
plot(V_LV_store_F_XB_SRX{5}{1},P_LV_store_F_XB_SRX{5}{1}-47,'linewidth',3,'color',[0.2+0.0111*4 0.1333+0.0889*4 0.5333+0.0111*4])
plot(V_LV_store_F_XB_SRX{6}{1},P_LV_store_F_XB_SRX{6}{1}-47,'linewidth',3,'color',[0.2+0.0111*5 0.1333+0.0889*5 0.5333+0.0111*5])
xlabel('Volume (mL)')
ylabel('Pressure (mmHg)')
legend('ATP','2% dATP','4% dATP','6% dATP','8% dATP','10% dATP')
set(gca, 'fontsize', 14)

% B
percent_dATP = [2 4 6 8 10 20];
figure
subplot(4,4,1)
hold on
plot(0,max_force_XB_SRX{1}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*0 0.1333+0.0889*0 0.5333+0.0111*0])
plot(percent_dATP(1),max_force_XB_SRX{2}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*1 0.1333+0.0889*1 0.5333+0.0111*1])
plot(percent_dATP(2),max_force_XB_SRX{3}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*2 0.1333+0.0889*2 0.5333+0.0111*2])
plot(percent_dATP(3),max_force_XB_SRX{4}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*3 0.1333+0.0889*3 0.5333+0.0111*3])
plot(percent_dATP(4),max_force_XB_SRX{5}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*4 0.1333+0.0889*4 0.5333+0.0111*4])
plot(percent_dATP(5),max_force_XB_SRX{6}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*5 0.1333+0.0889*5 0.5333+0.0111*5])
xlabel('Percent dATP')
ylabel('Max Force (kPa)')

subplot(4,4,2)
hold on
plot(0,CO_XB_SRX{1}{1}/454,'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*0 0.1333+0.0889*0 0.5333+0.0111*0])
plot(percent_dATP(1),CO_XB_SRX{2}{1}/454,'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*1 0.1333+0.0889*1 0.5333+0.0111*1])
plot(percent_dATP(2),CO_XB_SRX{3}{1}/454,'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*2 0.1333+0.0889*2 0.5333+0.0111*2])
plot(percent_dATP(3),CO_XB_SRX{4}{1}/454,'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*3 0.1333+0.0889*3 0.5333+0.0111*3])
plot(percent_dATP(4),CO_XB_SRX{5}{1}/454,'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*4 0.1333+0.0889*4 0.5333+0.0111*4])
plot(percent_dATP(5),CO_XB_SRX{6}{1}/454,'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*5 0.1333+0.0889*5 0.5333+0.0111*5])
xlabel('Percent dATP')
ylabel('SV (mL/beat)')

subplot(4,4,3)
hold on
plot(0,EF_XB_SRX{1}{1}*100,'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*0 0.1333+0.0889*0 0.5333+0.0111*0])
plot(percent_dATP(1),EF_XB_SRX{2}{1}*100,'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*1 0.1333+0.0889*1 0.5333+0.0111*1])
plot(percent_dATP(2),EF_XB_SRX{3}{1}*100,'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*2 0.1333+0.0889*2 0.5333+0.0111*2])
plot(percent_dATP(3),EF_XB_SRX{4}{1}*100,'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*3 0.1333+0.0889*3 0.5333+0.0111*3])
plot(percent_dATP(4),EF_XB_SRX{5}{1}*100,'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*4 0.1333+0.0889*4 0.5333+0.0111*4])
plot(percent_dATP(5),EF_XB_SRX{6}{1}*100,'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*5 0.1333+0.0889*5 0.5333+0.0111*5])
xlabel('Percent dATP')
ylabel('EF (%)')

subplot(4,4,4)
hold on
plot(0,CO_XB_SRX{1}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*0 0.1333+0.0889*0 0.5333+0.0111*0])
plot(percent_dATP(1),CO_XB_SRX{2}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*1 0.1333+0.0889*1 0.5333+0.0111*1])
plot(percent_dATP(2),CO_XB_SRX{3}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*2 0.1333+0.0889*2 0.5333+0.0111*2])
plot(percent_dATP(3),CO_XB_SRX{4}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*3 0.1333+0.0889*3 0.5333+0.0111*3])
plot(percent_dATP(4),CO_XB_SRX{5}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*4 0.1333+0.0889*4 0.5333+0.0111*4])
plot(percent_dATP(5),CO_XB_SRX{6}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*5 0.1333+0.0889*5 0.5333+0.0111*5])
xlabel('Percent dATP')
ylabel('CO (mL/min)')

subplot(4,4,5)
hold on
plot(0,LVDP_XB_SRX{1}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*0 0.1333+0.0889*0 0.5333+0.0111*0])
plot(percent_dATP(1),LVDP_XB_SRX{2}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*1 0.1333+0.0889*1 0.5333+0.0111*1])
plot(percent_dATP(2),LVDP_XB_SRX{3}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*2 0.1333+0.0889*2 0.5333+0.0111*2])
plot(percent_dATP(3),LVDP_XB_SRX{4}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*3 0.1333+0.0889*3 0.5333+0.0111*3])
plot(percent_dATP(4),LVDP_XB_SRX{5}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*4 0.1333+0.0889*4 0.5333+0.0111*4])
plot(percent_dATP(5),LVDP_XB_SRX{6}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*5 0.1333+0.0889*5 0.5333+0.0111*5])
xlabel('Percent dATP')
ylabel('LVDP (mmHg)')

subplot(4,4,6)
hold on
plot(0,plus_dPdt_F_XB_SRX{1}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*0 0.1333+0.0889*0 0.5333+0.0111*0])
plot(percent_dATP(1),plus_dPdt_F_XB_SRX{2}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*1 0.1333+0.0889*1 0.5333+0.0111*1])
plot(percent_dATP(2),plus_dPdt_F_XB_SRX{3}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*2 0.1333+0.0889*2 0.5333+0.0111*2])
plot(percent_dATP(3),plus_dPdt_F_XB_SRX{4}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*3 0.1333+0.0889*3 0.5333+0.0111*3])
plot(percent_dATP(4),plus_dPdt_F_XB_SRX{5}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*4 0.1333+0.0889*4 0.5333+0.0111*4])
plot(percent_dATP(5),plus_dPdt_F_XB_SRX{6}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*5 0.1333+0.0889*5 0.5333+0.0111*5])
xlabel('Percent dATP')
ylabel('+dP/dt (mmHg/s)')

subplot(4,4,7)
hold on
plot(0,minus_dPdt_F_XB_SRX{1}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*0 0.1333+0.0889*0 0.5333+0.0111*0])
plot(percent_dATP(1),minus_dPdt_F_XB_SRX{2}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*1 0.1333+0.0889*1 0.5333+0.0111*1])
plot(percent_dATP(2),minus_dPdt_F_XB_SRX{3}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*2 0.1333+0.0889*2 0.5333+0.0111*2])
plot(percent_dATP(3),minus_dPdt_F_XB_SRX{4}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*3 0.1333+0.0889*3 0.5333+0.0111*3])
plot(percent_dATP(4),minus_dPdt_F_XB_SRX{5}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*4 0.1333+0.0889*4 0.5333+0.0111*4])
plot(percent_dATP(5),minus_dPdt_F_XB_SRX{6}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*5 0.1333+0.0889*5 0.5333+0.0111*5])
xlabel('Percent dATP')
ylabel('-dP/dt (mmHg/s)')

subplot(4,4,8)
hold on
plot(0,work_rate_F_XB_SRX{1}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*0 0.1333+0.0889*0 0.5333+0.0111*0])
plot(percent_dATP(1),work_rate_F_XB_SRX{2}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*1 0.1333+0.0889*1 0.5333+0.0111*1])
plot(percent_dATP(2),work_rate_F_XB_SRX{3}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*2 0.1333+0.0889*2 0.5333+0.0111*2])
plot(percent_dATP(3),work_rate_F_XB_SRX{4}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*3 0.1333+0.0889*3 0.5333+0.0111*3])
plot(percent_dATP(4),work_rate_F_XB_SRX{5}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*4 0.1333+0.0889*4 0.5333+0.0111*4])
plot(percent_dATP(5),work_rate_F_XB_SRX{6}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*5 0.1333+0.0889*5 0.5333+0.0111*5])
xlabel('Percent dATP')
ylabel('Work Rate (mmHg*mL/s)')

subplot(4,4,9)
hold on
plot(0,ATP_F_XB_SRX{1}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*0 0.1333+0.0889*0 0.5333+0.0111*0])
plot(percent_dATP(1),ATP_F_XB_SRX{2}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*1 0.1333+0.0889*1 0.5333+0.0111*1])
plot(percent_dATP(2),ATP_F_XB_SRX{3}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*2 0.1333+0.0889*2 0.5333+0.0111*2])
plot(percent_dATP(3),ATP_F_XB_SRX{4}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*3 0.1333+0.0889*3 0.5333+0.0111*3])
plot(percent_dATP(4),ATP_F_XB_SRX{5}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*4 0.1333+0.0889*4 0.5333+0.0111*4])
plot(percent_dATP(5),ATP_F_XB_SRX{6}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*5 0.1333+0.0889*5 0.5333+0.0111*5])
xlabel('Percent dATP')
ylabel('ATP (mM)')

subplot(4,4,10)
hold on
plot(0,ADP_F_XB_SRX{1}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*0 0.1333+0.0889*0 0.5333+0.0111*0])
plot(percent_dATP(1),ADP_F_XB_SRX{2}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*1 0.1333+0.0889*1 0.5333+0.0111*1])
plot(percent_dATP(2),ADP_F_XB_SRX{3}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*2 0.1333+0.0889*2 0.5333+0.0111*2])
plot(percent_dATP(3),ADP_F_XB_SRX{4}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*3 0.1333+0.0889*3 0.5333+0.0111*3])
plot(percent_dATP(4),ADP_F_XB_SRX{5}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*4 0.1333+0.0889*4 0.5333+0.0111*4])
plot(percent_dATP(5),ADP_F_XB_SRX{6}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*5 0.1333+0.0889*5 0.5333+0.0111*5])
xlabel('Percent dATP')
ylabel('ADP (\muM)')

subplot(4,4,11)
hold on
plot(0,Pi_F_XB_SRX{1}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*0 0.1333+0.0889*0 0.5333+0.0111*0])
plot(percent_dATP(1),Pi_F_XB_SRX{2}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*1 0.1333+0.0889*1 0.5333+0.0111*1])
plot(percent_dATP(2),Pi_F_XB_SRX{3}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*2 0.1333+0.0889*2 0.5333+0.0111*2])
plot(percent_dATP(3),Pi_F_XB_SRX{4}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*3 0.1333+0.0889*3 0.5333+0.0111*3])
plot(percent_dATP(4),Pi_F_XB_SRX{5}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*4 0.1333+0.0889*4 0.5333+0.0111*4])
plot(percent_dATP(5),Pi_F_XB_SRX{6}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*5 0.1333+0.0889*5 0.5333+0.0111*5])
xlabel('Percent dATP')
ylabel('Pi (mM)')

subplot(4,4,12)
hold on
plot(0,MVO2_F_XB_SRX{1}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*0 0.1333+0.0889*0 0.5333+0.0111*0])
plot(percent_dATP(1),MVO2_F_XB_SRX{2}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*1 0.1333+0.0889*1 0.5333+0.0111*1])
plot(percent_dATP(2),MVO2_F_XB_SRX{3}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*2 0.1333+0.0889*2 0.5333+0.0111*2])
plot(percent_dATP(3),MVO2_F_XB_SRX{4}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*3 0.1333+0.0889*3 0.5333+0.0111*3])
plot(percent_dATP(4),MVO2_F_XB_SRX{5}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*4 0.1333+0.0889*4 0.5333+0.0111*4])
plot(percent_dATP(5),MVO2_F_XB_SRX{6}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*5 0.1333+0.0889*5 0.5333+0.0111*5])
xlabel('Percent dATP')
ylabel('MVO_{2} (uM O_{2}/min/g tissue)')

subplot(4,4,13)
hold on
plot(0,PCrATP_F_XB_SRX{1}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*0 0.1333+0.0889*0 0.5333+0.0111*0])
plot(percent_dATP(1),PCrATP_F_XB_SRX{2}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*1 0.1333+0.0889*1 0.5333+0.0111*1])
plot(percent_dATP(2),PCrATP_F_XB_SRX{3}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*2 0.1333+0.0889*2 0.5333+0.0111*2])
plot(percent_dATP(3),PCrATP_F_XB_SRX{4}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*3 0.1333+0.0889*3 0.5333+0.0111*3])
plot(percent_dATP(4),PCrATP_F_XB_SRX{5}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*4 0.1333+0.0889*4 0.5333+0.0111*4])
plot(percent_dATP(5),PCrATP_F_XB_SRX{6}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*5 0.1333+0.0889*5 0.5333+0.0111*5])
xlabel('Percent dATP')
ylabel('CrP/ATP ratio')

subplot(4,4,14)
hold on
plot(0,XB_turnover_F_XB_SRX{1}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*0 0.1333+0.0889*0 0.5333+0.0111*0])
plot(percent_dATP(1),XB_turnover_F_XB_SRX{2}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*1 0.1333+0.0889*1 0.5333+0.0111*1])
plot(percent_dATP(2),XB_turnover_F_XB_SRX{3}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*2 0.1333+0.0889*2 0.5333+0.0111*2])
plot(percent_dATP(3),XB_turnover_F_XB_SRX{4}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*3 0.1333+0.0889*3 0.5333+0.0111*3])
plot(percent_dATP(4),XB_turnover_F_XB_SRX{5}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*4 0.1333+0.0889*4 0.5333+0.0111*4])
plot(percent_dATP(5),XB_turnover_F_XB_SRX{6}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*5 0.1333+0.0889*5 0.5333+0.0111*5])
xlabel('Percent dATP')
ylabel('XB Cycling (1/s)')

subplot(4,4,15)
hold on
plot(0,ATPase_F_XB_SRX{1}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*0 0.1333+0.0889*0 0.5333+0.0111*0])
plot(percent_dATP(1),ATPase_F_XB_SRX{2}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*1 0.1333+0.0889*1 0.5333+0.0111*1])
plot(percent_dATP(2),ATPase_F_XB_SRX{3}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*2 0.1333+0.0889*2 0.5333+0.0111*2])
plot(percent_dATP(3),ATPase_F_XB_SRX{4}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*3 0.1333+0.0889*3 0.5333+0.0111*3])
plot(percent_dATP(4),ATPase_F_XB_SRX{5}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*4 0.1333+0.0889*4 0.5333+0.0111*4])
plot(percent_dATP(5),ATPase_F_XB_SRX{6}{1},'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*5 0.1333+0.0889*5 0.5333+0.0111*5])
xlabel('Percent dATP')
ylabel('ATPase (M/s/L)')

subplot(4,4,16)
hold on
plot(0,efficiency_F_XB_SRX{1}{1}./1000,'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*0 0.1333+0.0889*0 0.5333+0.0111*0])
plot(percent_dATP(1),efficiency_F_XB_SRX{2}{1}./1000,'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*1 0.1333+0.0889*1 0.5333+0.0111*1])
plot(percent_dATP(2),efficiency_F_XB_SRX{3}{1}./1000,'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*2 0.1333+0.0889*2 0.5333+0.0111*2])
plot(percent_dATP(3),efficiency_F_XB_SRX{4}{1}./1000,'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*3 0.1333+0.0889*3 0.5333+0.0111*3])
plot(percent_dATP(4),efficiency_F_XB_SRX{5}{1}./1000,'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*4 0.1333+0.0889*4 0.5333+0.0111*4])
plot(percent_dATP(5),efficiency_F_XB_SRX{6}{1}./1000,'o','markersize',10,'linewidth',3,'color',[0.2+0.0111*5 0.1333+0.0889*5 0.5333+0.0111*5])
xlabel('Percent dATP')
ylabel('Efficiency (mL^{2}*mmHg*ms/M)')
