function [Difference,LHS,RHS]  = Fsys(x,xminus,y,yminus,xss,yss,p)

%xss=[sstate.K; sstate.qk; sstate.q; sstate.a; sstate.int; sstate.w; 0];

%yss=[sstate.pit; sstate.pitw; sstate.mc; sstate.N;
 %    sstate.PId; sstate.G;
 %   sstate.Inv; sstate.ra; sstate.C; sstate.pk;
 %   sstate.sy; sstate.qk; sstate.Y; sstate.d];

  yminus=exp(yss+yminus);
  y=exp(yss+y);
  xminus=exp(xss+xminus);
  x=exp(xss+x);
 
    controls=char('pit','pitw','mc','N','PId','G','Inv','ra','C','pk','sy','Eqk','Y','d');
    states=char('K','qk','q','a','int','W','eint');
    delog=char('int','pk','ra','pit','d','sy','eint','pitw');

% controls t


for i=1:p.numcontrols
    eval(sprintf('%sminus = yminus(%f);',strtrim(controls(i,:)),i))
end

% controls t+1

for i=1:p.numcontrols
    eval(sprintf('%s = y(%f);',strtrim(controls(i,:)),i))
end

% state t-1

for i=1:p.numstates
    eval(sprintf('%sminus = xminus(%f);',strtrim(states(i,:)),i))
end

% state t

for i=1:p.numstates
    eval(sprintf('%s = x(%f);',strtrim(states(i,:)),i))
end
    
% variabels needing logging

for i=1:size(delog,1)
        eval(sprintf('%s = log(%s);',strtrim(delog(i,:)),delog(i,:)))
        eval(sprintf('%sminus = log(%sminus);',strtrim(delog(i,:)),strtrim(delog(i,:))))
end

 % other
    wss=exp(xss(6));
    Nss=exp(yss(4));
    Gss=exp(yss(6));
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LHS=zeros(p.numcontrols+p.numstates,1);
RHS=LHS;

% (1) Euler Bonds

LHS(1)=Cminus^-p.sigma;
RHS(1)=p.beta*(1+int)/(1+pit)*C^-p.sigma;

 % (2) Euler Capital
 
     %Adjustment frictions
     chid=p.chi0+p.chi2*p.chi1*abs(dminus)^(p.chi2-1)*aminus^(1-p.chi2);
    chidprime=p.chi0+p.chi2*p.chi1*abs(d)^(p.chi2-1)*a^(1-p.chi2);
     chiaprime=(1-p.chi2)*p.chi1*abs(d)^p.chi2*exp(a)^(-p.chi2);
     chi=p.chi0*abs(dminus) + p.chi1*abs(dminus)^p.chi2*aminus^(1-p.chi2);
 

 LHS(2)=(1+chid)*Cminus^-p.sigma;
 RHS(2)=p.beta*((1+ra)*(1+chidprime)-chiaprime)*C^-p.sigma;
 
 % (3) Wages / Labour supply
 
 muwhat=p.psi*log(Nminus/Nss) - log(W/wss);
 
 LHS(3)=pitwminus;
 RHS(3)=p.beta*pitw+p.kappa_w*(muwhat);
 

 % (4) production
 
 LHS(4)=Yminus;
 RHS(4)=(Kminus^p.alpha*(Nminus)^(1-p.alpha))^p.thetay;
 
 % (5) philips curve
 
 LHS(5)=pitminus;
 RHS(5)=p.beta*pit+p.kappa_p*(log(p.mu_p)+log(mcminus));
 
 % (6) labour demand
 
 LHS(6)=W;
 RHS(6)=mcminus*p.thetay*(1-p.alpha)*Yminus/Nminus;

 % (7) capital demand
 
 LHS(7)=pkminus;
 RHS(7)=mcminus*p.thetay*p.alpha*Yminus/Kminus;
 
%  (8) Dividend
 
 divy=Yminus-W*Nminus-Kminus*pkminus;
 
 LHS(8)=PIdminus;
 RHS(8)=divy;
 
% (9) capital price

 LHS(9)=qk;
 RHS(9)=1+p.tau*(Invminus/Kminus-p.delta);
 
% (10) expected capital price
 
 LHS(10)=Eqkminus;
 RHS(10)=qk;
 
% (11) expected capital price
 
 LHS(11)=Eqkminus;
 RHS(11)=(1/(1+ra))*(pk-p.tau/2*(Inv/K-p.delta)^2+p.tau*(Inv/K-p.delta)*Inv/K+(1-p.delta)*Eqk);
 
 % (12) capital lom
 
 LHS(12)=K;
 RHS(12)=(1-p.delta)*Kminus+Invminus-p.tau/2*(Invminus/Kminus-p.delta)^2*Kminus;
 
% (13) fund value

 LHS(13)=a;
 RHS(13)=K*qk+q;

% (14) ra
 
 LHS(14)=raminus*aminus;
 RHS(14)=Kminus*qk*(1-p.delta)-Kminus*qkminus+pkminus*Kminus+PIdminus +q-qminus;

% (15) fund lom
 
 LHS(15)=a;
 RHS(15)=(1+raminus)*aminus+dminus;
 
% (16) Government spending

LHS(16)=Gminus;
RHS(16)=Gss;
 
% (17) resource constraint
 
 LHS(17)=Yminus;
 RHS(17)=Cminus + Gminus + Invminus + p.tau/2*(Invminus/Kminus-p.delta)^2*Kminus+0*chi;
 
% (18) Interest rate

 LHS(18)=int;
 RHS(18)=p.rhoint*intminus+ (1-p.rhoint)*(p.phipi*(pitminus-p.pitstar)+p.intstar) + eintminus;
 
% (19)

    LHS(19)=pitw;
    RHS(19)=log(W/Wminus);
 
 % (20) labor share
 
 LHS(20)=syminus;
 RHS(20)=W*Nminus/Yminus;

 % (21 shock)
 
 RHS(21)=eint;
 LHS(21)=0*eintminus;

%% Difference
Difference=((LHS-RHS));

end
