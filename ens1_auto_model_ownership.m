% Paper 1 - Alexandrino da Silva & Guillen (2021) - Sovereign Holdings 
% clear all

% MAIN MODEL

format compact
outputloss = 'flat';  
denomination = 'LC'
nominal = 'on'
%ownership = 'domestic' % [domestic,external]

% load tpm_ny41_rho70_sig104_w4.mat ygrid pai  
load tpm.mat ygrid pai  
rho=0.70;
epsilonvar=0.026^2;


ny = numel(ygrid);
y = ( exp(ygrid(:)) );

nitermax = 200; % number of iterations

% parâmetros fixos
rstar = 0.04;   
theta = 0.5;  
sigg  = 2;  
meta  = 1.0; 
omega = 0.5; %0.23
yn=3.3;
cn=yn;
mu = 0; % Foreign share of sovereign bonds



gamma = 1.3 %gamma_auto; % Change this
kappa  = 0 %kappa_auto; % Change this

if strcmp(outputloss,'flat') % Returns 1 if the terms are equal, 0 otherwise
betta = 0.77 %beta_auto;    % 0.544   %0.80    # You may change this
a0    = -0.88 %-a0_auto;  %-0.849  %-0.89     # You mau change this
a1 = 1;
a2 = 0;
end

ya =  y - max(0,a0+ a1*y + a2*y.^2);        % y in autarky

% Creating grid "dn" for bonds "d"

dnupper = 0.6;         
dnlower = 0;
ndn = 61; 
dn = dnlower:(dnupper-dnlower)/(ndn-1):dnupper;
dn = dn(:);

% force the element closest to zero to be exactly zero
[xx,nd0]=min(abs(dn)); 
dn(nd0) = 0; 

% If is FC, bond is not subject to inflation

if strcmp(denomination,'FC')
    nominal = 'off'
end

% Creating grid for prices, when bond is nominal

if strcmp(nominal,'on')
pupper = 1.6; 
plower = 1.0;
np = 121; 
p = plower:(pupper-plower)/(np-1):pupper;
p = p(:); 
end

if strcmp(nominal,'off')
p=1;
np=1;
end

nd = ndn; % calling nd the lenght of bond grid
n = ny*nd; % total number of states (bonds and shocks)
nc = nd*np; % total number of nominal bonds

epsilon = -rho*repmat(log(y),1,ny)+repmat(log(y)',ny,1);
M = [    exp(-rstar - kappa*epsilon -0.5*(kappa^2)*epsilonvar)    ];

    
% pra facilitar adaptar o código
ndr=1;

ytry = ( repmat(y,nd,nc) );
dntry = repmat(dn',[ny ndr]);
dntry = reshape(dntry,n,1);
dntry = repmat(dntry,1,nc);

dnptry = ( repmat(dn',n,ndr*np) ); 

ptry = repmat(p',nd,1);  
ptry = reshape(ptry,1,nc);
ptry = repmat(ptry,n,1);

% qntiltry = . q~ = q/RER     = zeros(n,nc);
qntiltry     = ones(n,nc)/(1+rstar)/omega;


ctry        = ones(n,nc);
utry        = zeros(n,nc);

% Atry        = ones(n,nc);
% Btry        = ones(n,nc);

% Autharky

ca=ya;
% ua = ( ca.^(1-sigg)-1)  / (1-sigg);  
ua = (   (  (ca.^omega) * (cn^(1-omega))   ).^(1-sigg)  - 1)  / (1-sigg);


cc = zeros(ny,nd);

vc = zeros(ny,nd);  
vcnew = vc;         
vg = zeros(ny,nd);  
vb = zeros(ny,1); 
vr =  zeros(ny,1); 
f = zeros(ny,nd);  
qntil = ones(ny,nd)/(1+rstar)/omega; 

cpix = ones(ny,nd);     
cpixnew = ones(ny,nd);  
dnpix = ones(ny,nd);     
dnpixnew = ones(ny,nd);  
ppix = ones(ny,nd);     
ppixnew = ones(ny,nd);  

ctodn = ( repmat(1:ndn,1,ndr*np) );
repmat(1:np,nd,1);
ctop = ( reshape(ans,1,nc) );
 
 
dist = 1;
niter=1;
tol = 1e-5;
while [dist>tol niter<nitermax];
    tic

qntiltry = repmat(qntil,nd,np);
    
Evg = pai*vg;
Evr = pai*vr;
Evb = pai*vb;


a=0;
distX=1;
tolX=1e-5;
damp=0.8;

if strcmp(denomination,'LC')
    
    if omega==0.5
        Atry = ytry + qntiltry.*dnptry;
        Btry = (1/omega) * ((1/yn)^(1-omega)) * (dntry./ptry);
        ctry = ( ( (Btry.^2 + 4*Atry).^0.5 - Btry ).^2 )/4;
        ctry(imag(ctry)~=0)=NaN;
        clear Atry Btry
   
    else
         
    while [distX>tolX a<200]
clear Atry Btry
ctrynew = ytry + qntiltry.*dnptry - (1/omega).*((ctry/yn).^(1-omega)).*(dntry./ptry);
ctrynew = (1-damp)*ctrynew+ damp*ctry;
distX = max( abs( ctrynew(:) - ctry(:) ));
ctry=ctrynew;
clear ctrynew
a=a+1;
controlX(a)=distX;
    end
ctry(imag(ctry)~=0)=NaN;
a
    if distX>tolX
warning('No solution to non-linear equation')
    break
    end   
end
end

if  strcmp(denomination,'FC')
        ctry =  dnptry .* qntiltry - dntry + ytry;
end

% utry = (ctry.^(1-sigg) -1)  / (1-sigg)  - (gamma/2)*( (ptry - meta).^2 ) ;
utry = (  (  (ctry.^omega) * (cn^(1-omega))  ).^(1-sigg)  -   1  )  / (1-sigg)  - (gamma/2)*( (ptry - meta).^2 ) ;
utry(ctry<=0) = -inf;
utry(ctry==NaN) = -inf;
[vcnew(:), cpixnew(:)]  = max(  utry + betta * repmat(Evg,nd,np),[],2 ) ;

vbnew = ua+betta*(theta*Evr + (1-theta) * Evb);
f = vcnew<repmat(vbnew, 1,nd);

dnpixnew = ctodn(cpixnew);
ppixnew = ctop(cpixnew);

auxcc=cpixnew(:);
for i=1:n
    auxcc2(i)=ctry( i, auxcc(i));
end
cc=reshape(auxcc2,ny,ndn);


if  strcmp(denomination,'LC')
% qntilnew = pai*( (1/omega)*(1 - f).*( (cc/yn).^(1-omega) )./p(ppixnew) )/(1+rstar); 
qntilnew = (pai.*M)*( (1/omega)*(1 - f).*( (cc/yn).^(1-omega) )./p(ppixnew) ); 
end


if  strcmp(denomination,'FC')
% qntilnew = pai*( (1 - f) )/(1+rstar); 
qntilnew = (pai.*M)*( (1 - f) ); 
end


dist1 = max(abs(qntilnew(:)-qntil(:)))
dist2 = max(abs(vcnew(:)-vc(:)));
dist3 = max(abs(vbnew(:)-vb(:)));
dist4 = max(abs(cpixnew(:)-cpix(:))) ;
dist5 = max(abs(dnpixnew(:)-dnpix(:))) ;
dist6 = max(abs(ppixnew(:)-ppix(:))) ;

dist  = max([dist1,dist2,dist3,dist4])

count1 = size( find ( abs(  qntilnew(:)-qntil(:)  ) >tol ) )
count2 = size( find ( abs(  vcnew(:)-vc(:)        ) >tol ) )
count3 = size( find ( abs(  vbnew(:)-vb(:)        ) >tol ) )
count4 = size( find ( abs(  cpixnew(:)-cpix(:)    ) >tol ) )
count5 = size( find ( abs(  dnpixnew(:)-dnpix(:)  ) >tol ) )
count6 = size( find ( abs(  ppixnew(:)-ppix(:)    ) >tol ) )



qntilold = qntil;
qntil = qntilnew;


vc = vcnew;
vb = vbnew;
vg = max(vc,repmat(vb, 1, nd));
vr = vg(:,nd0); 

cpix = cpixnew;
dnpixold = dnpix;
dnpix = dnpixnew;
ppixold=ppix;
ppix = ppixnew;


control1(niter)=dist1;
control2(niter)=dist2;
control3(niter)=dist3;
control4(niter)=dist4;
control5(niter)=dist5;
control6(niter)=dist6;
controla(niter)=a;

niter=niter+1

toc
end


% qntilfree = pai*( (1/omega)*((cc/yn).^(1-omega))./p(ppixnew) )/(1+rstar); 
qntilfree = (pai.*M)*( (1/omega)*((cc/yn).^(1-omega))./p(ppixnew) ); 
qnntilfree =  exp(-rstar)*pai*( (1/omega)*((cc/yn).^(1-omega))./p(ppixnew) ); 


if  strcmp(denomination,'LC')
qnntil = exp(-rstar)*pai*( (1/omega)*(1 - f).*( (cc/yn).^(1-omega) )./p(ppixnew) ); 
end

if  strcmp(denomination,'FC')
qnntil = exp(-rstar)*pai*( (1 - f) ); 
end

dnpix = ctodn(cpix);
ppix  = ctop(cpix);

ppix=reshape(ppix,ny,ndn,ndr);
dnpix=reshape(dnpix,ny,ndn,ndr);

qntil = reshape(qntil,ny,ndn,ndr);
qntilfree = reshape(qntilfree,ny,ndn,ndr);
qnntil = reshape(qnntil,ny,ndn,ndr);

f  = reshape(f,ny,ndn,ndr);
ca = repmat(ca,[1,ndn,ndr]);

dnpc = dn(dnpix);
ppc = p(ppix);

eppc = pai*ppc;
eyt = pai*y;

rerc         = omega*((cc/yn).^(omega-1));
relc         = ( omega /(1-omega) )./(cc/yn);
pndivpc       = (1-omega) * ( (cc/yn).^omega ) ;

rera         = omega*((ca/yn).^(omega-1));
rela         = ( omega /(1-omega) )./(ca/yn);
pndivpa       = (1-omega)*((ca/yn).^omega) ;

if  strcmp(denomination,'LC')
qn = rerc.*qntil;
qnn = rerc.*qnntil;
end

if  strcmp(denomination,'FC')
qn=qntil;
qnn=qnntil;
end

qnfree = rerc.*qntilfree;
qnnfree = rerc.*qnntilfree;

maxppc = max(ppc(:))
maxdnpc = max(dnpc(:))

if dist<tol
    convergiu=1;
else
    convergiu=0;
end

convergiu

clear ans *try *new *tryix  auxy auxdn auxdr I
% eval(['save ens1_nominal_teste02'])

nome = [denomination, '_', nominal, '_', num2str(ny), '_', num2str(ndn),'_', num2str(np),...
'_betta', num2str(betta*1000), '_a0', num2str(-a0*1000), '_gamma', num2str(gamma*1000),...
'_kappa', num2str(kappa*100),'_var',num2str((epsilonvar^0.5)*1000)]

eval(['save ens1_' nome])




