clear all


pars



%% steady state

% bisection method

pkl=1/p.beta+p.delta-1-0.02;
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
    
%  marginal cost

	sstate.mc=log(1/(p.mu_p-(sstate.pit-p.beta*sstate.pit)/p.kappa_p));
    
%  Production labour

	sstate.N=1;

%  technology

	sstate.zy=0;
    
%  capital

	sstate.K=log((pk/(exp(sstate.mc)*p.alpha*p.thetay)*exp(sstate.N)^-((1-p.alpha)*p.thetay))^(1/(p.alpha*p.thetay-1)));

   %  firm output

	sstate.Y=log((exp(sstate.K)^p.alpha*exp(sstate.N)^(1-p.alpha))^p.thetay); % firm level production

%  wages
 
	sstate.w=log(exp(sstate.mc)*p.thetay*(1-p.alpha)*exp(sstate.Y)/exp(sstate.N)); 
         
%  dividend
  
    piy=exp(sstate.Y)-exp(sstate.w)*exp(sstate.N)-pk*exp(sstate.K);
    
    sstate.PId=log(piy);
 
% investment    
    
    sstate.Inv=log(p.delta*exp(sstate.K));

%  no arbritrage

    sstate.qk=log(1+p.tau*(exp(sstate.Inv)/exp(sstate.K)-p.delta));

	sstate.ra=(pk+exp(sstate.qk)*(1-p.delta)) / exp(sstate.qk) -1;

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
 
% difference

	diff= (1+sstate.ra) - 1/p.beta - euler2(sstate.d,p,sstate);

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
sstate.pitw=0;

NKsstate=sstate;
save('NKsstate.mat','NKsstate')



%% dynamics (SGU Form)

xss=[sstate.K; sstate.qk; sstate.q; sstate.a; sstate.int; sstate.w; 0];

yss=[sstate.pit; sstate.pitw; sstate.mc; sstate.N;
     sstate.PId; sstate.G;
    sstate.Inv; sstate.ra; sstate.C; sstate.pk;
    sstate.sy; sstate.qk; sstate.Y; sstate.d];
    
p.numcontrols=size(yss,1);
p.numstates=size(xss,1);

F = @(a,b,c,d)Fsys_NK(a,b,c,d,xss,yss,p);

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

figure(101)
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% other functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
  	function x = euler2(d,p,sstate)
	
	xd=p.chi0+p.chi2*p.chi1*abs(d)^(p.chi2-1)*exp(sstate.a)^(1-p.chi2);
    xa=(1-p.chi2)*p.chi1*abs(d)^p.chi2*exp(sstate.a)^(-p.chi2);
	x=xa/(1+xd);
end
  