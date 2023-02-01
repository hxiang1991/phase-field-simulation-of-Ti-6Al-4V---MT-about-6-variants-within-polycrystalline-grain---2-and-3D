function [Sigmay2, e0_pl]= yieldstress(sumphiA2,sigmaA_pl,sigmaB_pl,kappaHP_pl,D_alpha_pl,G_pl,delta_Gm,Vmol,...
                                       JCa_pl,JCb_pl,JCn_pl,JCm1_pl,JCm2_pl,T,T0,Tr,e1_pl,e2_pl,e3_pl,e4_pl,e5_pl,e6_pl)
                           
%% == yield stress from modified Hall-Petch function with volume fraction of alpha variant
sigmay_HP= (sigmaA_pl* sumphiA2+ sigmaB_pl*(1- sumphiA2)+ kappaHP_pl/D_alpha_pl^0.5)* G_pl;

%% == plastic strain hardening by modified Johnson-Cook (MJC) model
%  == equivalent plastic strain e0_pl;
%  == refer: https://doi.org/10.1016/B978-0-12-385204-5.00008-2; eqn 8.2

% e0_pl= 2^0.5/3* ((e1_pl- e2_pl).^2+ 6* e3_pl.^2).^0.5;

%  == refer: https://www.researchgate.net/post/How-can-i-calculate-the-equivalent-plastic-strain
Dev_e1_pl= e1_pl-(e1_pl+ e2_pl+ e3_pl)/3;                   % components of the deviator of plastic strain
Dev_e2_pl= e2_pl-(e1_pl+ e2_pl+ e3_pl)/3;
Dev_e3_pl= e3_pl-(e1_pl+ e2_pl+ e3_pl)/3;
Dev_e4_pl= e4_pl; Dev_e5_pl= e5_pl; Dev_e6_pl= e6_pl;

y2= (Dev_e1_pl.^2+ Dev_e2_pl.^2+ Dev_e3_pl.^2+...
 2* (Dev_e4_pl.^2+ Dev_e5_pl.^2+ Dev_e6_pl.^2))/2;
e0_pl= (4/3* y2).^0.5;                                      % equivalent plastic strain

%  == yield stress from plastic strain hardening based on MJC model
sigmay_MJC= (JCa_pl+ JCb_pl* (1+ JCm1_pl* log(T0/Tr))* e0_pl.^JCn_pl)* (1- T^JCm2_pl);

%% == total yield stress
Sigmay= sigmay_HP+ sigmay_MJC;   Sigmay= Sigmay* 1e06;       % unit: Pa

Sigmay= Sigmay/(delta_Gm/Vmol)/1.5;
Sigmay2= Sigmay.^2;

end   % end function