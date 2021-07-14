
close all
clear all


% NK model with two sectors. capital, capital adjustment cost and full modelling of
% illiquid fund with frictions. [ if liquidity fricionts are zero same as
% BaskNKcap

% Parameters:


p.bbeta=1-0.025/4;
p.ssigma=4.0;	% crra
p.eeps=20;	% epsilon
p.ppsi=1.0;	% labour elast (inverse frish)

p.mu=p.eeps/(p.eeps-1);

p.aalpha=1-2/3/p.mu;	% capital elasticity
p.ddelta=0.054/4;
p.tthetay=1.00;	% scale Y
p.tthetan=0.90; % scale N

p.ttau=11.4;

p.phipi=1.5; % taylor rule

p.Np=4; % average price duration

p.ttheta=p.Np/(p.Np+1); % calvo param

p.pphi=p.ttheta*(p.eeps-1)/((1-p.ttheta)*(1-p.bbeta*p.ttheta)); % price adjustment cost

p.pitstar=0.00/4; % inflation target

p.intstar=(1/p.bbeta)*(1+p.pitstar)-1; % r*

% adjustment costs 
p.cchi0=0.04;
p.cchi1=0.90;
p.cchi2=1.3;

% stocashtic parameters

p.rhoint=0.8;
p.rhozy=0.9;

p.seint=0.01;
p.sea=0.01;

% bisection method

pkl=1/p.bbeta+p.ddelta-1-0.02;
pkh=0.15;

diff=10;
tol1=10e-8;
iter=0;

while abs(diff)>tol1

iter=iter+1;

%  pk 

	pk=0.5*(pkl+pkh);

%  inflation

   sstate.pit=p.pitstar;
   
%  interest rate

   sstate.int=p.intstar; %(16)

%  Expansion labour

	sstate.le=log(1/15);
    
%  marginal cost

	sstate.mc=log( 1/p.eeps*((p.eeps-1)  +p.pphi*(sstate.pit)*(sstate.pit+1)-p.pphi*p.bbeta*(sstate.pit+1)*sstate.pit)); % marginal cost
    
%  Measure of goods

    p.zn=1; 

    sstate.N=log(p.zn*exp(sstate.le)^p.tthetan);

%  Production labour

	sstate.ly=log(	exp(sstate.mc)*(1-p.aalpha)*p.tthetay*exp(sstate.le)/((1-exp(sstate.mc))*p.tthetan*exp(sstate.N))	);

% labour

    sstate.l=log(exp(sstate.N)*exp(sstate.ly)+exp(sstate.le));
    
%  technology

	sstate.zy=0;
    
%  capital

	sstate.ki=log((pk/(exp(sstate.mc)*p.aalpha*p.tthetay)*exp(sstate.ly)^-((1-p.aalpha)*p.tthetay))^(1/(p.aalpha*p.tthetay-1)));

	sstate.k=log(exp(sstate.N)*exp(sstate.ki));

   %  firm output

	sstate.y=log((exp(sstate.ki)^p.aalpha*exp(sstate.ly)^(1-p.aalpha))^p.tthetay); % firm level production

%  wages
 
	sstate.w=log(exp(sstate.mc)*p.tthetay*(1-p.aalpha)*exp(sstate.y)/exp(sstate.ly)); 
         
%  dividend
  
    piy=exp(sstate.N)*(exp(sstate.y)*exp(sstate.mc)-exp(sstate.w)*exp(sstate.ly)-pk*exp(sstate.ki));
    pin=exp(sstate.N)*(exp(sstate.y)*(1-exp(sstate.mc)))-exp(sstate.w)*exp(sstate.le)-p.pphi/2*(sstate.pit)^2;
    
    sstate.PId=log(piy+pin);
 
% investment    
    
    sstate.Inv=log(p.ddelta*exp(sstate.k));

%  no arbritrage

    sstate.qk=log(1+p.ttau*(exp(sstate.Inv)/exp(sstate.k)-p.ddelta));

	sstate.ra=(pk+exp(sstate.qk)*(1-p.ddelta)) / exp(sstate.qk) -1;

% stock price

	sstate.q=log(exp(sstate.PId)/sstate.ra);

% fund value

	sstate.a=log(exp(sstate.q)+exp(sstate.k)*exp(sstate.qk));

% deposit

    sstate.d=-sstate.ra*exp(sstate.a);
    
% consumption

    adjcost=p.cchi0*abs(sstate.d) + p.cchi1*abs(sstate.d)^p.cchi2*exp(sstate.a)^(1-p.cchi2);

	sstate.c=log(exp(sstate.w)*exp(sstate.l)-sstate.d-adjcost);

 % other variables
 
    sstate.Y=log(exp(sstate.N)*exp(sstate.y));

    sstate.sy=exp(sstate.l)*exp(sstate.w)/(exp(sstate.Y));
 
% difference

	diff= (1+sstate.ra) - 1/p.bbeta - euler2(sstate.d,p,sstate);

if diff>0

	pkh=pk;

else

	pkl=pk;

end


if iter>200
break
end

end



sstate.pk=pk;
 
 %p;
 %sstate;
 
 p.lss=sstate.l;
p.wss=sstate.w;
 
 save('sstate','sstate')
  save('p','p')

%% SGU Form

xss=[sstate.k; sstate.qk; sstate.q; sstate.a; sstate.int; 0];

yss=[sstate.pit; sstate.le; sstate.mc; sstate.N; sstate.ly;
    sstate.l; sstate.ki; sstate.y; sstate.w; sstate.PId;
    sstate.Inv; sstate.ra; sstate.c; sstate.pk;
    sstate.sy; sstate.qk; sstate.Y; sstate.d];
    
p.numcontrols=size(yss,1);
p.numstates=size(xss,1);

F = @(a,b,c,d)Fsys(a,b,c,d,xss,yss,p);

[Fss,LHS,RHS]=F(xss*0,0*xss,0*yss,0*yss);

if abs(norm(Fss,'inf'))>10e-4
    [Fss LHS RHS]
    error('steady state not solved')
end
    
p.overrideEigen=true;

%% Solve RE via Schmitt-Grohe-Uribe Form
[hx,gx,F1,F2,F3,F4,par] = SGU_solver(F,p,p);

%%


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


irflist=char('RB','PI','Y','I','LS','N');
irfind=[5,7,23,17,21,12]';

figure(100)
clf

for i=1:6
    
subplot(3,2,i)

plot(IRF_state_sparse(irfind(i),:))
eval(sprintf('title("%s")',irflist(i,:)))
hline=refline(0,0);
hline.Color='black'
end
saveas(gcf,'IRF6.jpg')

     
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% other functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
  	function x = euler2(d,p,sstate)
	
	xd=p.cchi0+p.cchi2*p.cchi1*abs(d)^(p.cchi2-1)*exp(sstate.a)^(1-p.cchi2);
    xa=(1-p.cchi2)*p.cchi1*abs(d)^p.cchi2*exp(sstate.a)^(-p.cchi2);
	x=xa/(1+xd);
end
  