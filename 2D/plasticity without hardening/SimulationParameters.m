clear; clc

format long

%% == parent beta grains and alpha variants
grainBs= 16; variants = 2;

%% == simulation system
nx= 1024; ny= nx; nxy= nx* ny;
lengthx= nx* 2.5e-8;                                          % length in x direction of simulation, unit: m
dx= lengthx/nx;                                               % mesh size, unit: m 

%% == simulation parameters
gamma= 50* 1e-3;                                              % interfacial energy, unit: J/m2
Vmol= 1e-5;                                                   % molar volume, unit: m3/mol
thick= 5* dx;                                                 % thickness of interface, unit: m
ttime= 0;                                                     % initial time
dt0= 2e-03;                                                   % delta time 
M0= 1.6e-07;                                                  % interfacial mobility, unit: J/(m3* s)
nstep= 1e04; nprint1= 1e02;                                   % for loop step and output step

kappa0= 3* gamma* thick/8^0.5;                                % gradient coefficient, unit: J/m

%% == Gibbs free energy constant
diffGibbs= -3.1375e02;                                        % driving force, unit: J/mol
g_barrier= 3* Vmol* gamma/(4*2^0.5* thick);                   % gibbs energy barrier, unit: J/mol

a0= 32* g_barrier;                                            % gibbs free energy coefficient
b0= 3* a0- 12* diffGibbs;
c0= 2* a0- 12* diffGibbs;

%% == dimensionless form
delta_Gm= -diffGibbs; delta_M= 1e-18;

a= a0/delta_Gm; b= b0/delta_Gm; c= c0/delta_Gm;

kappa= kappa0* Vmol/(delta_Gm* dx^2);
dt= dt0* delta_Gm* delta_M/ dx^2;
M= M0* dx^2/(delta_M* Vmol);     

%% == elastic strain parameters
ang= 2* pi* [0.6978 0.3171 0.9502 0.0344...                   % orientation of parent beta grains
             0.4387 0.3816 0.7655 0.7952...
             0.1869 0.4898 0.4456 0.6463...
             0.7094 0.7547 0.2760 0.6797];                                 

e0(:,:,1)= [-0.0490, 0; 0,  0.0670];                          % transformation strain for 2 variants
e0(:,:,2)= [ 0.0670, 0; 0, -0.0490]; 

c11= 97.7e9/(delta_Gm/Vmol);                                  % elastic constants in dimensionless form
c12= 82.7e9/(delta_Gm/Vmol);
c44= 37.5e9/(delta_Gm/Vmol);

e11= zeros(nx,ny); e22= zeros(nx,ny);                         % pre-set blank array for strain and stress
e12= zeros(nx,ny); e21= zeros(nx,ny); 
s11= zeros(nx,ny); s22= zeros(nx,ny);
s12= zeros(nx,ny); s21= zeros(nx,ny);

%% == plastic strain/stress parameters
% == parameters for modified Hall-Petch function;              _pl: plastic
% == refer: https://doi.org/10.1016/j.matdes.2018.09.028
% D_alpha_pl= 0.3;                                               % average thickness of alpha lath from SLM,unit: mm,refer:https://doi.org/10.1016/j.actamat.2014.11.028
% sigmaA_pl= 550; sigmaB_pl= 1350;                               % lattice friction stress of alpha and beta phase, unit: MPa
% kappaHP_pl= 300;                                               % Hall-Petch coefficient, unit: MPa
% kappa_pl= 0.23; n_pl= 0.4;                                     % fitting parameters from experimental data
% epsilon_pl= 1;                                                 % strain rate, assume to be 1
% mu_pl= (54- 0.03* T0)* 1e9;                                    % shear modulus, unit: Pa;
% b_pl= 2.9e-10;                                                 % unit: m
 
k_pl= 1e-03;                                                   % constant for fourth-order plastic kinetic coefficient tensor
 
% G_pl= (kappa_pl* mu_pl* b_pl^3/(Kb* T0* log(1e7)))^n_pl;
% 
% % == modified Johnson-Cook(MJC) model for plastic strain hardening
% % == refer: https://doi.org/10.3390/ma11060938
% JCa_pl= 0.92e03; JCb_pl= 0.4e03;                               % fitting parameters for MJC model, unit: MPa
% JCn_pl= 0.578; JCm1_pl= 0.1578; JCm2_pl= 0.633; 
% Tr= 298; Tm= 1878;                                             % room temperature and melt temperature, unit: K
% 
% T= (T0- Tr)/(Tm- Tr);

Sigmay= 880e06;                                                % initial yield stress, unit: Pa
Sigmay= Sigmay/(delta_Gm/Vmol);                                % dimensionless form

Sigmay2= Sigmay^2;

e1_pl= zeros(nx,ny);                                           % initial plastic strain 
e2_pl= zeros(nx,ny); 
e3_pl= zeros(nx,ny); 

%% == pre-set matrix
tmpphiAplot= zeros(nx, ny); phiAplot= zeros(nx, ny, variants, grainBs); 
sumphiAplot1= zeros(nx, ny); sumphiAplot2= zeros(nx, ny);

[sym_x, sym_y] = ndgrid(1: nx, 1: ny);                          % 2d spacial coordinate
tmp_x = reshape(sym_x, nx* ny, 1);
tmp_y = reshape(sym_y, nx* ny, 1);
sym_cor_mat = [tmp_x tmp_y zeros(nxy,1)];

vflag=nstep/nprint1; VolF = zeros(vflag+ 1, variants+ 1);
VolF(2:vflag+ 1,1)= (1: vflag)* nprint1;





