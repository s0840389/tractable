

%% households

p.beta=0.99; % discount
p.sigma=4.0;	% crra
p.psi=2.0;	% labour elast (inverse frish)

p.epsilon_p=7;	% epsilon prices
p.mu_p=p.epsilon_p/(p.epsilon_p-1); % price markup

p.epsilon_w=11;	% epsilon wagws
p.mu_w=p.epsilon_w/(p.epsilon_w-1); % wage markup

p.shareE=0.2; % share of expansionary workers

%% firms

p.alpha=1-2/3*p.mu_p;	% capital elasticity
p.delta=0.07/4;
p.thetay=1.00;	% scale Y
p.thetan=0.90; % scale N

%% government

p.phipi=1.5; % taylor rule

p.rhoint=0.8; % smoothing

p.pitstar=0.00/4; % inflation target

p.intstar=(1/p.beta)*(1+p.pitstar)-1; % r*

p.gshr=0.2; % gov purchases


%% frictions

p.tau=25; % investment adjustment costs

p.Np=3; % average price duration

p.theta=p.Np/(p.Np+1); % calvo param


p.phi=p.theta*(p.epsilon_p-1)/((1-p.theta)*(1-p.beta*p.theta)); % price adjustment cost

p.kappa_p=(1-p.theta)*(1-p.beta*p.theta)/p.theta;

p.Nw=2; % average wage duration

p.thetaw=p.Nw/(p.Nw+1); % calvo wage param

p.kappa_w=(1-p.thetaw)*(1-p.beta*p.thetaw)/(p.thetaw*(1+p.epsilon_w*p.psi));

p.chi0=0.05;
p.chi1=0.8;
p.chi2=1.4;

%% stocashtic parameters

p.rhozy=0.9;

p.se_int=0.01;
p.se_zy=0.01;

%% other parameters

p.Nss=1; % steady state labour
p.Eshr=0.2; % share of expansionary labour
