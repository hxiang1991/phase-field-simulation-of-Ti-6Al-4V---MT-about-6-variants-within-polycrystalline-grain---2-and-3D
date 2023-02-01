clear; clc
format short

nx= 1024; ny= nx; 
grainBs= 16; variants= 2;

rho= 0.25; 

load('BetaMap_16p2v1024grid.mat');
phiB2= phiB.^2; sumphiB2= sum(phiB2, 3);
maxsumphiB2= max(sumphiB2,[],'all');
minsumphiB2= min(sumphiB2,[],'all');

psi= (rho* (maxsumphiB2- sumphiB2)+ (sumphiB2- minsumphiB2))/(maxsumphiB2- minsumphiB2);
inrange= (psi~= 1);

video1= VideoWriter('Microstructure evolution_16p2v.avi'); 
video1.FrameRate = 5; video1.Quality= 100; open(video1);
 
video2= VideoWriter('equivalent stress.avi'); 
video2.FrameRate = 5; video2.Quality= 100; open(video2);

video3= VideoWriter('Elastic energy.avi'); 
video3.FrameRate = 5; video3.Quality= 100; open(video3);

for i= 1: 100
    
    load(['01_ElasticEnergy_ttime_', num2str(i*0.1004*2, '%6.4f'),'s.mat']);  
    load(['01_sumphiAplot2_ttime_', num2str(i*0.1004*2, '%6.4f'),'s.mat']);
    load(['01_VonMisesStress_ttime_', num2str(i*0.1004*2, '%6.4f'),'s.mat']);

    sumphiAplot2(inrange)= 0;
    f1= figure('visible', 'off');clf; 
    imagesc(sumphiAplot2); grid off; axis equal; axis([1 nx 1 ny]);colormap('cool');colorbar;
    title(['time: ', num2str(i*0.1004*2, '%6.4f'),'s']);
    frame1 = getframe(gcf); writeVideo(video1,frame1); close(f1);

    f2= figure('visible', 'off');clf
    imagesc(output_Vons); grid off; axis equal; axis([1 nx 1 ny]); colormap('jet'); colorbar;
    title(['time: ', num2str(i*0.1004*2, '%6.4f'),'s']);
    frame2 = getframe(gcf); writeVideo(video2,frame2); close(f2);

    f3= figure('visible', 'off');clf
    imagesc(E_elastic); grid off; axis equal; axis([1 nx 1 ny]); colormap('jet'); colorbar;
    title(['time: ', num2str(i*0.1004*2, '%6.4f'),'s']);
    frame3 = getframe(gcf); writeVideo(video3,frame3); close(f3);

end
close(video1); close(video2); close(video3);

