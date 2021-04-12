% Paper 1 - Alexandrino da Silva & Guillen (2021) - Sovereign Holdings 
% clear all

% STATISTICS

% clear all, format compact

 %load simu_ens1_nominal_teste02 
load ens1_FC_off_200_61_1_betta770_a0880_gamma1300_kappa0_var26
%load ens1_LC_on_200_61_121_betta770_a0880_gamma1300_kappa0_var26.mat

denomination
nominal
sigg
gamma 
betta
a0
meta
ny
ndn
ndr
np

x = find(B==0&F==0);
TBYT = ((YT - C))./ YT; 
% TBY = ((YT - C))./ Y; 
% ou 
TBY = (RER.*(YT - C))./ Y; %quase não afeta resultados

gdp = Ytilde(x);
gdptrade = YTtilde(x);
tbyt = TBYT(x);
tby = TBY(x);
constrade = C(x);
rer = RER(x);
infl = (PI(x)-1)*100;
einfl = (EPI(x)-1)*100;

if strcmp(denomination,'LC')
dy = D(x)./( PI(x) .* Ytilde(x) )*100;
end

if strcmp(denomination,'FC')
    dy =D(x).*RER(x)./(Ytilde(x) )*100;
end


default_frequency = mean(F)*100
default_frequency_good_standing = mean(F(B==0))*100

mean_dy     = mean(dy)
mean_tby   = mean(tby)*100    %%%%%%%%%%%%%%%%%%%
mean_rer = mean(rer)
mean_infl = mean(infl)
mean_loggdp      = mean(log(gdp))*100;
mean_loggdptrade = mean(log(gdptrade))*100;
mean_logconstrade     = mean(log(constrade))*100;
mean_PM  = mean(PM(x))
mean_PMN = mean(PMN(x))
mean_PMNe  = mean(PMNe(x))
mean_PMNNe = mean(PMNNe(x))

std_dy      = std(dy)
std_tby    = std(tby)*100    %%%%%%%%%%%%%%%%%%%
std_rer = std(rer)*100
std_gdp = std(gdp)*100;
std_gdptrade = std(gdptrade)*100;
std_infl  = std(infl)
std_loggdp             = std(log(gdp))*100
std_loggdptrade        = std(log(gdptrade))*100;
std_logconstrade            = std(log(constrade))*100;
ratio_std_constrade_gdptrade = std_logconstrade/std_loggdptrade
std_PM   = std(PM(x))
std_PMN  = std(PMN(x))
std_PMNe   = std(PMNe(x))
std_PMNNe  = std(PMNNe(x))

corr_dy_loggdp      = corr( [ dy log(gdp)])
corr_tby_loggdp    = corr([ tby log(gdp)])
corr_logconstrade_loggdp = corr( [ log(constrade) log(gdp)])
corr_infl_loggdp     = corr([infl log(gdp)])
corr_rer_loggdp  = corr( [ rer log(gdp)])
corr_pm_loggdp      = corr([PM(x) log(gdp)])
corr_pmn_loggdp     = corr([PMN(x) log(gdp)])
corr_pmne_loggdp      = corr([PMNe(x) log(gdp)])
corr_pmnne_loggdp     = corr([PMNNe(x) log(gdp)])


corr_dy_rer      = corr( [ dy rer])
corr_tby_rer    = corr([ tby rer])
corr_logconstrade_rer = corr( [ log(constrade) rer])
corr_loggdp_rer     = corr([log(gdp) rer])
corr_infl_rer  = corr( [ infl rer])
corr_pm_rer      = corr([PM(x) rer])
corr_pmn_rer     = corr([PMN(x) rer])
corr_pmne_rer      = corr([PMNe(x) rer])
corr_pmnne_rer     = corr([PMNNe(x) rer])


corr_dy_infl      = corr( [ dy infl])
corr_tby_infl    = corr([ tby infl])
corr_logconstrade_infl = corr( [ log(constrade) infl])
corr_loggdp_infl     = corr([log(gdp) infl])
corr_rer_infl  = corr( [ rer infl])
corr_pm_infl      = corr([PM(x) infl])
corr_pmn_infl     = corr([PMN(x) infl])
corr_pmne_infl      = corr([PMNe(x) infl])
corr_pmnne_infl     = corr([PMNNe(x) infl])


% corr_dy_loggdptrade      = corr( [ dy log(gdptrade)])
% % corr_edy_loggdptrade     = corr( [ edy log(gdptrade)])
% % corr_cambio_loggdptrade  = corr( [ cambio log(gdptrade)])
% corr_tby_loggdptrade    = corr([ tby log(gdptrade)])
% corr_pm_loggdptrade      = corr([PM(x) log(gdptrade)])
% corr_pmn_loggdptrade     = corr([PMN(x) log(gdptrade)])
% corr_logcons_loggdptrade = corr( [ log(cons) log(gdptrade)])
% corr_infl_loggdptrade     = corr([infl log(gdptrade)])
% corr_rer_loggdptrade  = corr( [ rer log(gdptrade)])


% filename = [nome, '.xlsx']
% % filename = aaa;
% 
% AAA = [
%     convergiu;
%     niter;
%     
%     default_frequency_good_standing; 
%     
%     mean_dy; 
%     mean_tby;
%     mean_infl;
%     mean_rer;
%     mean_PM; 
%     mean_PMN; 
%     mean_PMNe; 
%     mean_PMNNe; 
% 
%     std_dy; 
%     std_tby;
%     std_infl;
%     std_rer;
%     std_loggdp; 
%     std_loggdptrade;
%     std_logconstrade; 
%     ratio_std_constrade_gdptrade ;
%     std_PM; 
%     std_PMN;
%     std_PMNe; 
%     std_PMNNe;
% 
%     corr_dy_loggdp(1,2); 
%     corr_tby_loggdp(1,2); 
%     corr_infl_loggdp(1,2);
%     corr_rer_loggdp(1,2);   
%     corr_logconstrade_loggdp(1,2);
%     corr_pm_loggdp(1,2); 
%     corr_pmn_loggdp(1,2);
%     corr_pmne_loggdp(1,2); 
%     corr_pmnne_loggdp(1,2);
%     
%     corr_dy_rer(1,2); 
%     corr_tby_rer(1,2); 
%     corr_loggdp_rer(1,2);
%     corr_infl_rer(1,2);   
%     corr_logconstrade_rer(1,2);
%     corr_pm_rer(1,2); 
%     corr_pmn_rer(1,2);
%     corr_pmne_rer(1,2); 
%     corr_pmnne_rer(1,2);
%  
%     corr_dy_infl(1,2); 
%     corr_tby_infl(1,2); 
%     corr_loggdp_infl(1,2);
%     corr_rer_infl(1,2);   
%     corr_logconstrade_infl(1,2);
%     corr_pm_infl(1,2); 
%     corr_pmn_infl(1,2);
%     corr_pmne_infl(1,2); 
%     corr_pmnne_infl(1,2);
% 
% 
%     ];

% xlswrite(filename, AAA)

% clear TBY TBYT tbyt tby gdp infl cons rer gdptrade dy D PI Ytilde YTtilde Y YT YA YTA B D C F RER NER Qfree

% eval(['save ens1_' nome])