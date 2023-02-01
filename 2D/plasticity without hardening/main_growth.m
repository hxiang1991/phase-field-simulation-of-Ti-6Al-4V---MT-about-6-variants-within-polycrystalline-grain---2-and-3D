clear; clc

format long;
 
%% == simulation paraters
SimulationParameters;

%% == map of Parent beta grains and a seed of alpha grain
% [phiB]= ParentGrains(nx, ny, grainBs);
load('BetaMap_16p2v1024grid.mat');

new_phiB= zeros(nx,ny,variants,grainBs);
for g= 1: grainBs
    for v= 1: variants

        new_phiB(:,:,v,g)= phiB(:,:,g);

    end
end

phiA= zeros(nx,ny,variants,grainBs);
for t= 1: grainBs
    if phiB(nx/2-1:nx/2-1,ny/2-1:ny/2-1,t)== ones(1,1)

       phiA(nx/2-1:nx/2-1,ny/2-1:ny/2-1,1,t)= 1;
       break

    end
end

if max(phiA,[],'all')== 0

    nstep= 0;

end

%% == pre-fft parameters
tmpkx= 2*pi*[0: nx/2 -nx/2+1: -1]/nx; tmpky= tmpkx;
[kx, ky]= ndgrid(tmpkx, tmpky);k2= kx.^2+ ky.^2; 
kx= kx./k2.^0.5; ky= ky./k2.^0.5;
kx(isnan(kx))= 0; ky(isnan(ky))= 0;

%% -- green's operator and eigen strain matrix of different variants in different parent beta grains
[tot_Cpq,hom_Cpq,het_Cpq,eigen,gomega,tot_Cpqrs,Spq] = greenmatrix_stiffnessmatrix_eigenstrain(grainBs,variants,phiB,ang,e0,nx,ny,kx,ky,c11,c12,c44,k_pl);

%% -- denominator in Allen-Cahn function
denom = 1+ dt* M* kappa* k2;

% %% -- video clips
% video1= VideoWriter('Microstructure evolution_16p2v.avi'); 
% video1.FrameRate = 10; video1.Quality= 100; open(video1);
%  
% video2= VideoWriter('equivalent stress.avi'); 
% video2.FrameRate = 10; video2.Quality= 100; open(video2);
% 
% video3= VideoWriter('equivalent plastic strain.avi'); 
% video3.FrameRate = 10; video3.Quality= 100; open(video3);
% 
% video4= VideoWriter('Elastic energy.avi'); 
% video4.FrameRate = 10; video4.Quality= 100; open(video4);

%% == microstructure evolution		
for istep= 1: nstep

    ttime= ttime + dt;

	phiA2= phiA.^2; sumphiA2= sum(phiA2, [3,4]); phiAk= fft(fft(phiA, [], 1), [], 2);
    
    if max(sumphiA2,[],'all')> 1

        sumphiA2= sumphiA2/max(sumphiA2,[],'all'); 

    end    
	
	dfdphiA= a* phiA- b* phiA2+ c* phiA.* sumphiA2;
	dfdphiAk= fft(fft(dfdphiA, [], 1), [], 2);
    
%     [Sigmay2, e0_pl]= yieldstress(sumphiA2,sigmaA_pl,sigmaB_pl,kappaHP_pl,D_alpha_pl,G_pl,delta_Gm,Vmol,JCa_pl,JCb_pl,JCn_pl,JCm1_pl,JCm2_pl,T,T0,Tr,e1_pl,e2_pl,e3_pl);
	
    [deldphiAk,E_elastic,e1_pl,e2_pl,e3_pl,output_s11,output_s22,output_s12] = solve_elasticity(phiA,eigen,e1_pl,e2_pl,e3_pl,grainBs,variants,...
                                                                                                tot_Cpq,hom_Cpq,het_Cpq,gomega,tot_Cpqrs,Spq,Sigmay2,nx,ny,kx,ky,dt,delta_Gm,Vmol);

    output_Vons= (output_s11.^2+ output_s22.^2- output_s11.* output_s22+ 3* output_s12.^2).^0.5;

    Dev_e1_pl= e1_pl-(e1_pl+ e2_pl)/2;                          % components of the deviator of plastic strain
    Dev_e2_pl= e2_pl-(e1_pl+ e2_pl)/2;
    Dev_e3_pl= e3_pl;

    y2= (Dev_e1_pl.^2+ Dev_e2_pl.^2+ 2* Dev_e3_pl.^2)/2;
    e0_pl= (4/3* y2).^0.5;                                      % equivalent plastic strain

    phiAk= (phiAk- dt* M.* (dfdphiAk+ deldphiAk))./denom;
	
	phiA= new_phiB.* real(ifft(ifft(phiAk, [], 1), [], 2));
	
	inrange= (phiA> 1); phiA(inrange)= 1; inrange= (phiA< 0); phiA(inrange)= 0;

	tmpphiA2= phiA.^2;
    for g= 1: grainBs
        for v= 1: variants
            
            tmpphiA2_vg= tmpphiA2(:,:,v,g);
            inrange= (tmpphiA2_vg> 0.5); tmpphiAplot(inrange)= 1;
            phiAplot(:,:,v,g)= tmpphiAplot;

            sumphiAplot1= sumphiAplot1+ phiAplot(:,:,v,g);
            sumphiAplot2= sumphiAplot2+ v* phiAplot(:,:,v,g);

            tmpphiAplot= zeros(nx, ny);

        end
    end

    if (mod(istep, nprint1)== 0)
        
%         f1= figure('visible', 'off');clf
%         imagesc(sumphiAplot2); grid off; axis equal; axis([1 nx 1 ny]); colormap('jet'); colorbar;
%         frame1 = getframe(gcf); writeVideo(video1,frame1); close(f1);
% 
%         f2= figure('visible', 'off');clf
%         imagesc(output_Vons); grid off; axis equal; axis([1 nx 1 ny]); colormap('jet'); colorbar;
%         frame2 = getframe(gcf); writeVideo(video2,frame2); close(f2);
% 
%         f3= figure('visible', 'off');clf
%         imagesc(e0_pl); grid off; axis equal; axis([1 nx 1 ny]); colormap('jet'); colorbar;
%         frame3 = getframe(gcf); writeVideo(video3,frame3); close(f3);
% 
%         f4= figure('visible', 'off');clf
%         imagesc(E_elastic); grid off; axis equal; axis([1 nx 1 ny]); colormap('jet'); colorbar;
%         frame4 = getframe(gcf); writeVideo(video4,frame4); close(f4);

        VolF(istep/nprint1+ 1,2)= sum(phiAplot(:,:,1,:), 'all')* 100/nxy;
        VolF(istep/nprint1+ 1,3)= sum(phiAplot(:,:,2,:), 'all')* 100/nxy;

%     end
% 
%     if (mod(istep, nprint2)== 0)

        filename1= ['phiA2_ttime_',num2str(ttime,'%6.4f'),'s.mat']; save(filename1, 'tmpphiA2');
        filename2= ['sumphiAplot2_ttime_',num2str(ttime,'%6.4f'),'s.mat']; save(filename2, 'sumphiAplot2');
        filename3= ['VonMisesStress_ttime_',num2str(ttime,'%6.4f'),'s.mat']; save(filename3, 'output_Vons');
        filename4= ['EquivalentPlasticStrain_ttime_',num2str(ttime,'%6.4f'),'s.mat']; save(filename4, 'e0_pl');
        filename5= ['ElasticEnergy_ttime_',num2str(ttime,'%6.4f'),'s.mat']; save(filename5, 'E_elastic');

    end
	
	maxphiA2= max(tmpphiA2,[],'all'); maxsumphiAplot= max(sumphiAplot1,[],'all'); 
	
	phiAplot= zeros(nx, ny, variants, grainBs); 
    sumphiAplot1= zeros(nx, ny); sumphiAplot2= zeros(nx, ny);	
	
	if (maxsumphiAplot> 1| maxphiA2< 2e-2| isinf(deldphiAk))			
		break				
	end

end

filename= 'VolumnFraction.mat'; save(filename, 'VolF'); 
% close(video1); close(video2); close(video3); close(video4);
