

clear all


pars
load('NKsstate.mat')


sstate=NKsstate; % default to NK model steadystate

% targets to be the same

    % 1) Labour share
    % 2) Capital ouput ratio
    % 3) Prices and interest rates
    % 4) frictions
    
    
 % matching labour share and ratios using wage foc's and ls equation
    
 lyle_ratio=p.shareE/(1-p.shareE);

 mc=exp(sstate.mc);
 
 scaley=1/mc*sstate.sy/(1+lyle_ratio); % thetay*(1-alphay)
 
 scalen=1/(1-mc)*sstate.sy/(1+1/lyle_ratio); % thetan*(1-alphan)
    
 p.thetan=scalen;
 p.thetay=1;
 p.alphay=1-scaley;
 
 sstate.le=log(0.2); % expansionary labour
 sstate.ly=log(0.8); % production labour
 sstate.Mg=0; % measure of goods
 
 p.zn=exp(sstate.Mg)/exp(sstate.le)^p.thetan; % expansionary technology
 
 % match capital output ratio
 
 sstate.ki=log((sstate.pk/(mc*p.thetay*p.alphay) * exp(sstate.ly)^(-p.thetay*(1-p.alphay)) )^(1/(p.thetay*p.alphay-1))); % firm capital
 
 sstate.K=log(exp(sstate.Mg)*exp(sstate.ki)); % total capital
 
 sstate.yi=log((exp(sstate.ki)^p.alphay * exp(sstate.ly)^(1-p.alphay))^p.thetay); % firm  output
 
 sstate.Y=sstate.yi+sstate.Mg; % total output
 
 sstate.N=log(exp(sstate.ly)*exp(sstate.Mg)+exp(sstate.le));
 
%  wages
 
	sstate.w=log(exp(sstate.mc)*p.thetay*(1-p.alphay)*exp(sstate.yi)/exp(sstate.ly)); 
         
%  dividend
  
    piy=exp(sstate.Mg)*(exp(sstate.yi)*exp(sstate.mc)-exp(sstate.w)*exp(sstate.ly)-sstate.pk*exp(sstate.ki));
    pin=exp(sstate.Mg)*(exp(sstate.yi)*(1-exp(sstate.mc)))-exp(sstate.w)*exp(sstate.le);
    
    sstate.PId=log(piy+pin);
 
% investment    
    
    sstate.Inv=log(p.delta*exp(sstate.K));

%  no arbritrage

    sstate.qk=log(1+p.tau*(exp(sstate.Inv)/exp(sstate.K)-p.delta));

	sstate.ra=(sstate.pk+exp(sstate.qk)*(1-p.delta)) / exp(sstate.qk) -1;

% stock price

	sstate.q=log(exp(sstate.PId)/sstate.ra);

% fund value

	sstate.a=log(exp(sstate.q)+exp(sstate.K)*exp(sstate.qk));

% deposit

    sstate.d=-sstate.ra*exp(sstate.a);
        
% consumption

    adjcost=p.chi0*abs(sstate.d) + p.chi1*abs(sstate.d)^p.chi2*exp(sstate.a)^(1-p.chi2);

% government

    sstate.G=log(exp(sstate.Y)*p.gshr);
    
	%sstate.C=log(exp(sstate.w)*exp(sstate.N)-sstate.d-adjcost);

    sstate.C=log(exp(sstate.Y)-exp(sstate.Inv)-exp(sstate.G));
        
 % other variables
 
    sstate.sy=exp(sstate.N)*exp(sstate.w)/(exp(sstate.Y));
 
    sstate.pitw=0;

YNsstate=sstate;

save('YNsstate.mat','YNsstate')



%% dynamics (SGU Form)

xss=[sstate.K; sstate.qk; sstate.q; sstate.a; sstate.int; sstate.w; 0];

yss=[sstate.pit; sstate.pitw; sstate.mc; sstate.N;
     sstate.PId; sstate.G; 
    sstate.Inv; sstate.ra; sstate.C; sstate.pk;
    sstate.sy; sstate.qk; sstate.Y; sstate.d;
    sstate.Mg; sstate.le; sstate.ly; sstate.ki; sstate.yi];
    
p.numcontrols=size(yss,1);
p.numstates=size(xss,1);

F = @(a,b,c,d)Fsys_YN(a,b,c,d,xss,yss,p);

[Fss,LHS,RHS]=F(xss*0,0*xss,0*yss,0*yss);

if abs(norm(Fss,'inf'))>10e-4
    [Fss LHS RHS]
    error('steady state not solved')
end
    
p.overrideEigen=true;

 %% Solve RE via Schmitt-Grohe-Uribe Form
[hx,gx,F1,F2,F3,F4,par] = SGU_solver(F,p,p);


%% Produce IRFs
disp('Calculating IRFs.');

mpar.maxlag=25; % Quarters

x0=zeros(p.numstates,1);
x0(end)=0.01;

MX=[eye(length(x0));gx];
IRF_state_sparse=[];
x=x0;

% y=ybar+gx*(x-xbar)
% x'=xbar+hx(x-xbar)


for t=1:mpar.maxlag
    IRF_state_sparse(:,t)=(MX*x)';
    x=hx*x;
end


irflist=char('RB','PI','Y','I','LS','W','N','PId','ra');
irfind=[5,8,20,14,18,6,11,12,15]';

figure(102)
clf

for i=1:9
    
subplot(3,3,i)

plot(IRF_state_sparse(irfind(i),:))
eval(sprintf('title("%s")',irflist(i,:)))
hline=refline(0,0);
hline.Color='black';
xlim([0,24])
end
saveas(gcf,'IRF6.jpg')


 
    