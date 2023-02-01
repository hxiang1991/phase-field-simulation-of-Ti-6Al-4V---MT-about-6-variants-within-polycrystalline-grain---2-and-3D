clear; clc

format long

%% == parent beta grains and alpha variants
grainBs= 8; variants = 4;

%% == simulation system
nx= 128; ny= nx; nz= nx; nxyz= nx* ny* nz;
lengthx= nx* 2.5e-8;                                                                                 % length in x direction of simulation, unit: m
dx= lengthx/nx;                                                                                      % mesh size, unit: m 

%% == simulation parameters
Kb= 1.3806e-23;                                                                                      % Boltzmann constant
T0= 1075;                                                                                            % temperature
gamma= 50* 1e-3;                                                                                     % interfacial energy, unit: J/m2
Vmol= 1e-5;                                                                                          % molar volume, unit: m3/mol
thick= 5* dx;                                                                                        % thickness of interface, unit: m
ttime= 0;                                                                                            % initial time
dt0= 2e-03;                                                                                          % delta time 
M0= 1.6e-07;                                                                                         % interfacial mobility, unit: J/(m3* s)
nstep= 1e04; nprint= 1e02;                                                                           % for loop step and output step

nucleationstep= 600;                                                                                 % nucleation step
ampnoise= 2.5e-3;                                                                                      % amplititude of noise for nucleation

rho= 0.25;                                                                                           % misfit strain relaxtion parameter at grain boundary

kappa0= 3* gamma* thick/8^0.5;                                                                       % gradient coefficient, unit: J/m

%% == Gibbs free energy constant
diffGibbs= -3.1375e02;                                                                               % driving force, unit: J/mol
g_barrier= 3* Vmol* gamma/(4*2^0.5* thick);                                                          % gibbs energy barrier, unit: J/mol

a0= 32* g_barrier;                                                                                   % gibbs free energy coefficient
b0= 3* a0- 12* diffGibbs;
c0= 2* a0- 12* diffGibbs;

%% == dimensionless form
delta_Gm= -diffGibbs; delta_M= 1e-18;

a= a0/delta_Gm; b= b0/delta_Gm; c= c0/delta_Gm;

kappa= kappa0* Vmol/(delta_Gm* dx^2);
dt= dt0* delta_Gm* delta_M/ dx^2;
M= M0* dx^2/(delta_M* Vmol);     

%% == elastic strain parameters
ang= 2* pi* [0.9543 0.6491 0.8516; 0.0222 0.8006 0.8781;...                                           % orientations of parent beta grains
             0.4079 0.4361 0.0402; 0.8290 0.6900 0.3390;...
             0.7425 0.0658 0.2568; 0.0940 0.0014 0.2869;...
             0.9527 0.6490 0.1307; 0.0886 0.0105 0.2120]; 

rot0= [0 1 0; -1/2^0.5 0 1/2^.5; 1/2^0.5 0 1/2^.5];

e0(:,:,1) = [-0.0490, -0.0031,       0; -0.0031,  0.0670,       0;       0,       0, -0.0003];        % case 4
e0(:,:,2) = [-0.0490,       0, -0.0031;       0, -0.0003,       0; -0.0031,       0,  0.0670];
e0(:,:,3) = [ 0.0334, -0.0222, -0.0253; -0.0222, -0.0100,  0.0412; -0.0253,  0.0412, -0.0056];
e0(:,:,4) = [ 0.0334, -0.0253, -0.0222; -0.0253, -0.0056,  0.0412; -0.0222,  0.0412, -0.0100];

c11= 97.7e9/(delta_Gm/Vmol); c12= 82.7e9/(delta_Gm/Vmol); c44= 37.5e9/(delta_Gm/Vmol);                 % elastic constants in dimensionless form

%% == plastic strain/stress parameters
% == parameters for modified Hall-Petch function;                                                     _pl: plastic
% == refer: https://doi.org/10.1016/j.matdes.2018.09.028
D_alpha_pl= 0.3;                                                                                      % average thickness of alpha lath from SLM,unit: mm,refer:https://doi.org/10.1016/j.actamat.2014.11.028
sigmaA_pl= 550; sigmaB_pl= 1350;                                                                      % lattice friction stress of alpha and beta phase, unit: MPa
kappaHP_pl= 300;                                                                                      % Hall-Petch coefficient, unit: MPa
kappa_pl= 0.23; n_pl= 0.4;                                                                            % fitting parameters from experimental data
epsilon_pl= 1;                                                                                        % strain rate, assume to be 1
mu_pl= (54- 0.03* T0)* 1e9;                                                                           % shear modulus, unit: Pa;
b_pl= 2.9e-10;                                                                                        % unit: m

k_pl= 1e-03;                                                                                          % constant for fourth-order plastic kinetic coefficient tensor

G_pl= (kappa_pl* mu_pl* b_pl^3/(Kb* T0* log(1e7)))^n_pl;

% == modified Johnson-Cook(MJC) model for plastic strain hardening
% == refer: https://doi.org/10.3390/ma11060938
JCa_pl= 0.92e03; JCb_pl= 0.4e03;                                                                      % fitting parameters for MJC model, unit: MPa
JCn_pl= 0.578; JCm1_pl= 0.1578; JCm2_pl= 0.633; 
Tr= 298; Tm= 1878;                                                                                    % room temperature and melt temperature, unit: K

T= (T0- Tr)/(Tm- Tr);

e1_pl= zeros(nx,ny,nz); e2_pl= zeros(nx,ny,nz); e3_pl= zeros(nx,ny,nz);                               % initial plastic strain 
e4_pl= zeros(nx,ny,nz); e5_pl= zeros(nx,ny,nz); e6_pl= zeros(nx,ny,nz); 

%% == external applied loading
s0yy= -50e06* ones(nx, ny, nz);                                                                         % unit: Pa
s0yy= s0yy/(delta_Gm/Vmol);

%% == pre-set matrix
tmpphiAplot= zeros(nx, ny, nz); phiAplot= zeros(nx, ny, nz, variants, grainBs); 
sumphiAplot1= zeros(nx, ny, nz); sumphiAplot2= zeros(nx, ny, nz);

[sym_x, sym_y, sym_z]= ndgrid(1: nx, 1: ny, 1: nz);                                                         % 2d spacial coordinate
tmp_x= reshape(sym_x, nxyz, 1);
tmp_y= reshape(sym_y, nxyz, 1);
tmp_z= reshape(sym_z, nxyz, 1);
sym_cor_mat = [tmp_x tmp_y tmp_z];

cc= hsv(variants* grainBs);
vflag=nstep/nprint; VolF = zeros(vflag+ 1, variants+ 1);
VolF(2:vflag+ 1,1)= (1: vflag)* nprint;





