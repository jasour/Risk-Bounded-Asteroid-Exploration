%% Risk Bounded Planning For Asteroid Exploration in Presence of Uncertain Model of Asteroid
% Set of all Impulsive Thrusts that satisfies probability of success,
% i.e.,{v: Probability(Radius of Periapsis >= Impact Radius) >= Delta}
% Ashkan Jasour, rarnop.mit.edu 
%%
clc;clear all;close all
%%
nx=3; % number of uncertain parameters + design parameter
d=5; % relaxation order: SDP uses 2*d number of the moments of uncertainties to calculate the Risk.

%% Safety Constraint 
%uncertain parameters: Density and Impact Radius
%design parameter: Impulsive Thrust
mpol('x',1,nx); 

r0=3; % initial orbit of satelite
v0=1; % initial velocity of satelite
rho=1; % nominal density of asteroide
drho= rho + x(1); % uncertain density of asteroide where x(1) is a random variable 
R=1; % nominal impact radius of asteroide
dR= R + x(2); % uncertain impact radius of asteroide where x(2) is a random variable 
mu=4/3*3.6*pi*(rho)*R^3;% gravitational parameter of asteroide
dmu=4/3*3.6*pi*drho*dR^3;% Uncertain gravitational parameter of asteroide
vd=x(3);% design parameter

% safety constraint: Radius of Periapsis >= Impact Radius
% rp= r0*2*dmu/(2*dmu-r0*(v0+vd)^2)-1; Radius of Periapsis
% rp >= R -----> (r0^2+R*r0)*(v0+vd)^2-2*R*mu>=0;
g=(r0^2+dR*r0)*(v0+vd)^2-2*dR*dmu; % g>=0: safety constraint in terms of uncertain parameters

%% Moments Information

%moments of Lebesgue Measure over [-1,1]^2 to calculate the integral
u=1;l=-1; yL=[2];for i=1:2*d ;yL(i+1,1)= ( u^(i+1) - l^(i+1) )/(i+1);end 
vpow=[];for k = 0:2*d; vpow = [vpow;genpow(nx,k)]; end; 
yL=prod(yL(vpow+1),2);

%moments of Uniform probability distribution over [0,1] for uncertain variable x1
u=0.5;l=-0.5;yx1=[1];for i=1:2*d ;yx1(i+1,1)=(1/(u-l))*((u^(i+1) - l^(i+1))/(i+1));end 
u=0.5;l=-0.5;yx2=[1];for i=1:2*d ;yx2(i+1,1)=(1/(u-l))*((u^(i+1) - l^(i+1))/(i+1));end 


%% GloptiPoly
% Gloptipoly solves the moment SDP.
% We use the solution of the dual SDP (i.e., c_mom) to calculate the risk.
% c_mom is the coefficients of the polynomial indicator function.
[c_mom]=func_Glopti(nx,g,d,yL);

syms v
vpow=[];for k = 0:2*d; vpow = [vpow;genpow(nx,k)]; end; 
yx1x2x3=yx1(vpow(:,1)+1).*yx2(vpow(:,2)+1).*v.^vpow(:,3);

% Obtained Success Curve
P=sum(yx1x2x3.*c_mom);

v=[-1:0.1:1]; plot(v,eval(P),'LineWidth',3); hold on
Delta=0.3; plot(v,1-Delta*ones(size(v,2)),'--');
ylabel('Probability of Success ');
xlabel('Impulsive Thrust');
str1 = '$ 1-\Delta_{risk} $';text(0.94,0.75,str1,'HorizontalAlignment','right','Interpreter','latex','FontSize',30) 


