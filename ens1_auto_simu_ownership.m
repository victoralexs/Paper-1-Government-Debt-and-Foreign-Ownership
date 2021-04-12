% Paper 1 - Alexandrino da Silva & Guillen (2021) - Sovereign Holdings 
% clear all

% SIMULATIONS

% clear all, format compact

load ens1_FC_off_200_61_1_betta770_a0880_gamma1300_kappa0_var26
%load ens1_LC_on_200_61_121_betta770_a0880_gamma1300_kappa0_var26.mat
 
cpai = cumsum(pai,2);

T = 1e5;
Tburn = 1e4;

[xx,i] = min(abs(y-1));
j = floor(nd/2);
b = 0; 

%Y* é PIB/Produto dividido por níveld e preços P

Y = zeros(Tburn+T,1);
Ytilde = zeros(Tburn+T,1);
YTtilde = zeros(Tburn+T,1);
YA = zeros(Tburn+T,1);
YT = zeros(Tburn+T,1);
YTA = zeros(Tburn+T,1);

B = zeros(Tburn+T,1);
D = zeros(Tburn+T,1);
C = zeros(Tburn+T,1);
F = zeros(Tburn+T,1);

NER = zeros(Tburn+T,1);
NERP = zeros(Tburn+T,1);
RER = zeros(Tburn+T,1);
RERP = zeros(Tburn+T,1);

Q = zeros(Tburn+T,1);  
Qfree = zeros(Tburn+T,1);  
QNe = zeros(Tburn+T,1);  
QNefree = zeros(Tburn+T,1);  

PM = zeros(Tburn+T,1); 
PMN = zeros(Tburn+T,1); 
PMNe = zeros(Tburn+T,1); 
PMNNe = zeros(Tburn+T,1); 

STATE = zeros(Tburn+T,1);
rstar_pct = rstar*100; 

PI = zeros(Tburn+T,1); 
EPI = zeros(Tburn+T,1); 
EYT = zeros(Tburn+T,1); 

RR = zeros(Tburn+T,1); 

for t=1:T+Tburn  

STATE(t,1) = sub2ind([ny nd ],i,j);

rr = rand;
RR(t,1) = rr;

F(t,1)      = f(i,j); 
D(t,1)      = dn(j); 
%Y* é PIB/Produto dividido por níveld e preços P
Y(t,1)      = rerc(i,j)*y(i)  + pndivpc(i,j)*yn;
YA(t,1)     = rera(i,j)*ya(i) + pndivpa(i,j)*yn;
YT(t,1)     = y(i);
YTA(t,1)    = ya(i);


if (b==0) & (F(t)==0) ; 
B(t,1) = 0; 
Ytilde(t,1) = Y(t);
YTtilde(t,1) = YT(t);
C(t,1) = cc(i,j);
RER(t,1) = rerc(i,j);
RERP(t,1) = rerc(i,j);
PI(t,1)      = ppc(i,j);
EPI(t,1)      = eppc(i,j);
EYT(t,1)      = eyt(i);
Q(t,1)= qn(i,dnpix(i,j));
Qfree(t,1)= qnfree(i,dnpix(i,j));
QNe(t,1)= qnn(i,dnpix(i,j));    %%%%%%%%%%%%%%%%%%%%%%
QNefree(t,1)= qnnfree(i,dnpix(i,j));    %%%%%%%%%%%%%%%%%%%%%%
jp  = dnpix(i,j); 
end 

if (b==0) & (F(t) ==1);
B(t,1) = 0; 
Ytilde(t,1) = YA(t);
YTtilde(t,1) = YTA(t);
C(t,1) = ya(i);
RER(t,1) = rera(i,j);
RERP(t,1) = rerc(i,j);
PI(t,1)      = 1;
EPI(t,1)      = 1;
EYT(t,1)      = eyt(i);
Q(t,1) = 0;
Qfree(t,1) = 0;
QNe(t,1)= 0;    %%%%%%%%%%%%%%%%%%%%%%
QNefree(t,1)= 0;    %%%%%%%%%%%%%%%%%%%%%%
jp = nd0; 
end 

if (b==1) & (rr>theta);  
B(t,1) = 1; 
Ytilde(t,1) = YA(t);
YTtilde(t,1) = YTA(t);
C(t,1) = ya(i);
RER(t,1) = rera(i,j);
RERP(t,1) = rerc(i,j);
PI(t,1)      = 1;
EPI(t,1)      = 1;
EYT(t,1)      = eyt(i);
Q(t,1) = 0;
Qfree(t,1) = 0;
QNe(t,1)= 0;    %%%%%%%%%%%%%%%%%%%%%%
QNefree(t,1)= 0;    %%%%%%%%%%%%%%%%%%%%%%
jp = nd0; 
end 

if (b==1) & (rr<=theta) ;
B(t,1) = 0; 
Ytilde(t,1) = Y(t);
YTtilde(t,1) = YT(t);
C(t,1) = cc(i,nd0);
RER(t,1) = rerc(i,nd0);
RERP(t,1) = rerc(i,nd0);
YT(t,1) = y(i);
PI(t,1)      = ppc(i,nd0);
EPI(t,1)      = eppc(i,nd0);
EYT(t,1)      = eyt(i);
Q(t,1)= qn(i,dnpix(i,nd0));
Qfree(t,1)= qnfree(i,dnpix(i,nd0));
QNe(t,1)= qnn(i,dnpix(i,nd0));    %%%%%%%%%%%%%%%%%%%%%%
QNefree(t,1)= qnnfree(i,dnpix(i,nd0));    %%%%%%%%%%%%%%%%%%%%%%
jp = dnpix(i,nd0);
end 

find(cpai(i,:)>rand); 
i = ans(1);
b = B(t) + F(t); 
j=jp; 
end

for t=1:T+Tburn-1
if strcmp(denomination,'FC') 
    PM(t,1)     = ( (  1/Q(t,1)    ) -1 )*100 - rstar_pct;
    PMNe(t,1)   = ( (  1/QNe(t,1)  ) -1 )*100 - rstar_pct;
end

if strcmp(denomination,'LC') 
%     PM(t,1)     =   (  Qfree(t,1)/Q(t,1) - 1)*(1+rstar)*100   ;
        PM(t,1)     =   (  1/Q(t,1) - 1/Qfree(t,1))*100   ;

    PMN(t,1)    = ( (  1/Q(t,1)  ) -1 )*100 - rstar_pct; 
    if isnan(PM(t,1))==1;
        PM(t,1)=inf;
    end
    
    PMNe(t,1)   =   (  QNefree(t,1)/QNe(t,1) - 1)*(1+rstar)*100   ; 
    PMNNe(t,1)  = ( (  1/QNe(t,1)  ) -1 )*100 - rstar_pct; 
    if isnan(PMNe(t,1))==1;
        PMNe(t,1)=inf;
    end
    
end
end 

STATE = STATE(Tburn+1:end);
Y = Y(Tburn+1:end);
YA = YA(Tburn+1:end);
YT = YT(Tburn+1:end);
YTA = YTA(Tburn+1:end);
Ytilde = Ytilde(Tburn+1:end);
YTtilde = YTtilde(Tburn+1:end);
B = B(Tburn+1:end);
D = D(Tburn+1:end);
C = C(Tburn+1:end);
F = F(Tburn+1:end);
RER = RER(Tburn+1:end); 
RERP = RERP(Tburn+1:end); 
NER = NER(Tburn+1:end); 
NERP = NERP(Tburn+1:end); 
Q = Q(Tburn+1:end);
Qfree = Qfree(Tburn+1:end);
QNe = QNe(Tburn+1:end);
QNefree = QNefree(Tburn+1:end);
PM = PM(Tburn+1:end);
PMN = PMN(Tburn+1:end);
PMNe = PMNe(Tburn+1:end);
PMNNe = PMNNe(Tburn+1:end);

PI = PI(Tburn+1:end);   
EPI = EPI(Tburn+1:end);   
EYT = EYT(Tburn+1:end);   

RR = RR(Tburn+1:end);   

for i=1:nd
lad(i,1) = mean(D==dn(i));
end


for i=1:np
lap(i,1) = mean(PI==p(i));
end
% 
% subplot(2,1,1)
% plot(dn,lad,'x-' ,'linewidth', 2)
% xlabel('d')
% ylabel('lad')
% shg
% 
% subplot(2,1,2)
% plot(p,lap,'x-' ,'linewidth', 2)
% xlabel('p')
% ylabel('lap')
% shg


% eval(['save simu_' denomination '_' outputloss ] )
% eval(['save simu_ens1_nominal_teste02'] )

eval(['save ens1_' nome])



