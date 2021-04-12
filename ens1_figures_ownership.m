% Paper 1 - Alexandrino da Silva & Guillen (2021) - Sovereign Holdings 
% clear all

% FIGURES

clear all
load ens1_FC_off_200_61_1_betta770_a0880_gamma1300_kappa0_var26


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ens1_FC_off_41_61_1_betta772_a0894_gamma130_kappa0
figure (1)
subplot(2,2,1)
plot(dn,f(1,:),'b:',     dn, f(21,:),'k',    dn,f(41,:),'k--',    'LineWidth',2)
title('A: Default','FontWeight','bold')
xlabel('Current Debt','FontWeight','bold')
ylabel('Default=1','FontWeight','bold')

subplot(2,2,2)
plot(dn,dnpc(1,:),'b:',     dn, dnpc(21,:),'k',    dn,dnpc(41,:),'k--',    'LineWidth',2)
title('B: Debt issuance','FontWeight','bold')
xlabel('Current Debt','FontWeight','bold')
ylabel('Next period Debt','FontWeight','bold')
legend('Low Output','Median Output', 'High Output')

reref = (1-f).*rerc + f.*rera;
subplot(2,2,3)
plot(dn,reref(1,:),'b:',     dn, reref(21,:),'k',    dn,reref(41,:),'k--',    'LineWidth',2)
title('C: Real exchange rate','FontWeight','bold')
xlabel('Current Debt','FontWeight','bold')
ylabel('RER','FontWeight','bold')

subplot(2,2,4)
plot(dn,ppc(1,:)-1,'b:',     dn, ppc(21,:)-1,'k',    dn,ppc(41,:)-1,'k--',    'LineWidth',2)
title('D: Inflation','FontWeight','bold')
xlabel('Current Debt','FontWeight','bold')
ylabel('Inflation','FontWeight','bold')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure (2)
plot(dn,qn(21,:),'b',     dn, qnfree(21,:),'k--',    'LineWidth',2)
title('Bond price schedule for median output','FontWeight','bold')
xlabel('Next period Debt','FontWeight','bold')
ylabel('Price','FontWeight','bold')
legend('Risky bond','Risk-free bond')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




qnFC = qn;
rerefFC = reref;
dnpcFC = dnpc;
% clear all
clearvars -except qnFC rerefFC dnpcFC
load ens1_LC_on_200_61_121_betta770_a0880_gamma1300_kappa0_var26




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ens1_FC_off_41_61_1_betta772_a0894_gamma130_kappa0
figure (3)
subplot(2,2,1)
plot(dn,f(1,:),'b:',     dn, f(21,:),'k',    dn,f(41,:),'k--',    'LineWidth',2)
title('A: Default','FontWeight','bold')
xlabel('Current Debt','FontWeight','bold')
ylabel('Default=1','FontWeight','bold')

subplot(2,2,2)
plot(dn,dnpc(1,:),'b:',     dn, dnpc(21,:),'k',    dn,dnpc(41,:),'k--',    'LineWidth',2)
title('B: Debt issuance','FontWeight','bold')
xlabel('Current Debt','FontWeight','bold')
ylabel('Next period Debt','FontWeight','bold')
legend('Low Output','Median Output', 'High Output')

reref = (1-f).*rerc + f.*rera;
subplot(2,2,3)
plot(dn,reref(1,:),'b:',     dn, reref(21,:),'k',    dn,reref(41,:),'k--',    'LineWidth',2)
title('C: Real exchange rate','FontWeight','bold')
xlabel('Current Debt','FontWeight','bold')
ylabel('RER','FontWeight','bold')

subplot(2,2,4)
plot(dn,ppc(1,:)-1,'b:',     dn, ppc(21,:)-1,'k',    dn,ppc(41,:)-1,'k--',    'LineWidth',2)
title('D: Inflation','FontWeight','bold')
xlabel('Current Debt','FontWeight','bold')
ylabel('Inflation','FontWeight','bold')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure (4)
% plot(dn,qn(11,:),'b--',     dn, qnFC(11,:),'k--',  dn,qn(21,:),'b',     dn, qnFC(21,:),'k',   'LineWidth',2)
plot( dn,qn(21,:),'b',     dn, qnFC(21,:),'k--',   'LineWidth',2)
title('Bond price schedule for median output','FontWeight','bold')
xlabel('Next period Debt','FontWeight','bold')
ylabel('Price','FontWeight','bold')
legend('LC bond','FC bond')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure (5)
% 
% subplot(3,1,1)
% plot(dn,qn(1,:),'b',     dn, qnfree(1,:),'k--',    'LineWidth',2)
% title('Bond price schedule for low output','FontWeight','bold')
% xlabel('Next period Debt','FontWeight','bold')
% ylabel('Price','FontWeight','bold')
% legend('Risky bond','Risk-Free bond')
% 
% subplot(3,1,2)
% plot(dn,qn(21,:),'b',     dn, qnfree(21,:),'k--',    'LineWidth',2)
% title('Bond price schedule for median output','FontWeight','bold')
% xlabel('Next period Debt','FontWeight','bold')
% ylabel('Price','FontWeight','bold')
% 
% subplot(3,1,3)
% plot(dn,qn(41,:),'b',     dn, qnfree(41,:),'k--',    'LineWidth',2)
% title('Bond price schedule for high output','FontWeight','bold')
% xlabel('Next period Debt','FontWeight','bold')
% ylabel('Price','FontWeight','bold')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure (6)
% subplot(2,1,1)
% plot(dn,qnfree(1,:),'b:',     dn, qnfree(21,:),'k',    dn, qnfree(41,:),'k--', 'LineWidth',2)
% title('Bond price ','FontWeight','bold')
% xlabel('Next period Debt','FontWeight','bold')
% ylabel('Price','FontWeight','bold')
% 
% subplot(2,1,2)
% plot(dn,qntilfree(1,:),'b:',     dn, qntilfree(21,:),'k--',    dn, qntilfree(41,:),'k--', 'LineWidth',2)
% title('Bond price adjusted by RER \q~','FontWeight','bold')
% xlabel('Next period Debt','FontWeight','bold')
% ylabel('Price','FontWeight','bold')
% legend('Low Output','Median Output', 'High Output')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ens1_FC_off_41_61_1_betta772_a0894_gamma130_kappa0
figure (7)

subplot(2,1,1)
plot(dn,dnpc(21,:),'b:',     dn, dnpcFC(21,:),'k',        'LineWidth',2)
title('A: Debt issuance','FontWeight','bold')
xlabel('Current Debt','FontWeight','bold')
ylabel('Next period Debt','FontWeight','bold')
legend('LC','FC')

reref = (1-f).*rerc + f.*rera;
subplot(2,1,2)
plot(dn,reref(21,:),'b:',     dn, rerefFC(21,:),'k',       'LineWidth',2)
title('B: Real exchange rate','FontWeight','bold')
xlabel('Current Debt','FontWeight','bold')
ylabel('RER','FontWeight','bold')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






