clear; clc

nx= 1024; ny= nx; 
No= [2 6 14 50 100 200];                % Number of selected frame
len= length(No);

load('BetaMap_16p2v1024grid.mat');
phiB2= phiB.^2; sumphiB2= sum(phiB2, 3);
maxsumphiB2= max(sumphiB2,[],'all');
minsumphiB2= min(sumphiB2,[],'all');

rho= 0.25; 
psi= (rho* (maxsumphiB2- sumphiB2)+ (sumphiB2- minsumphiB2))/(maxsumphiB2- minsumphiB2);
% inrange= (psi~= 1);
inrange= (psi<0.85);

E_elastic0= zeros(nx, ny); sumphiAplot20= zeros(nx, ny); 
output_Vons0= zeros(nx, ny); e0_pl0= zeros(nx, ny); 

for i= 1: len

    t= No(i);

    load(['01_ElasticEnergy_ttime_', num2str(t*0.1004, '%6.4f'),'s.mat']);  
    load(['01_sumphiAplot2_ttime_', num2str(t*0.1004, '%6.4f'),'s.mat']);
    load(['01_VonMisesStress_ttime_', num2str(t*0.1004, '%6.4f'),'s.mat']);
    load(['EquivalentPlasticStrain_ttime_', num2str(t*0.1004/2, '%6.4f'),'s.mat']);

    filename1=strcat('ElasticEnergy_ttime_',num2str(t*0.1004,'_GBs%6.4f'),'.tiff');
    filename2=strcat('sumphiAplot2_ttime_',num2str(t*0.1004,'_GBs%6.4f'),'.tiff');
    filename3=strcat('VonMisesStress_ttime_',num2str(t*0.1004,'_GBs%6.4f'),'.tiff');
    filename4=strcat('EquivalentPlasticStrain_ttime_',num2str(t*0.1004,'_GBs%6.4f'),'.tiff');

    maxE= max(E_elastic, [], 'all'); E_elastic0(inrange)= maxE;
    E_elastic(inrange)= 0;
    f1= figure('visible', 'off');clf; 
    imagesc(E_elastic0); grid off; axis equal; axis([1 nx 1 ny]);colormap('white');hold on
    imagesc(E_elastic); grid off; axis equal; axis([1 nx 1 ny]);colormap('jet'); alpha 0.75; colorbar;
    title(['time: ', num2str(t*0.1004, '%6.4f'),'s']);
    frame1= getframe(f1); imwrite(frame1.cdata, filename1);

    maxphiAplot2= max(sumphiAplot2, [], 'all'); sumphiAplot20(inrange)= maxphiAplot2;
    sumphiAplot2(inrange)= 0;
    f2= figure('visible', 'off');clf; 
    imagesc(sumphiAplot20); grid off; axis equal; axis([1 nx 1 ny]);colormap('white');hold on
    imagesc(sumphiAplot2); grid off; axis equal; axis([1 nx 1 ny]);colormap('cool');alpha 0.75;c=colorbar;c.Limits= [0 2];
    title(['time: ', num2str(t*0.1004, '%6.4f'),'s']);
    frame2= getframe(f2); imwrite(frame2.cdata, filename2); hold off

    maxVons= max(output_Vons, [], 'all'); output_Vons0(inrange)= maxVons;
    output_Vons(inrange)= 0;
    f3= figure('visible', 'off');clf; 
    imagesc(output_Vons0); grid off; axis equal; axis([1 nx 1 ny]);colormap('white');hold on
    imagesc(output_Vons); grid off; axis equal; axis([1 nx 1 ny]);colormap('jet');alpha 0.75;colorbar;
    title(['time: ', num2str(t*0.1004, '%6.4f'),'s']);
    frame3= getframe(f3); imwrite(frame3.cdata, filename3);

    maxe0_pl= max(e0_pl, [], 'all');e0_pl0(inrange)= maxe0_pl;
    e0_pl(inrange)= 0;
    f4= figure('visible', 'off');clf; 
    imagesc(e0_pl0); grid off; axis equal; axis([1 nx 1 ny]);colormap('white');hold on
    imagesc(e0_pl); grid off; axis equal; axis([1 nx 1 ny]);colormap('jet');alpha 0.75;colorbar;
    title(['time: ', num2str(t*0.1004, '%6.4f'),'s']);
    frame4= getframe(f4); imwrite(frame4.cdata, filename4);
end
