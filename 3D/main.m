clear; clc

format long;
 
%% == simulation paraters
SimulationParameters;

%% == map of Parent beta grains and a seed of alpha grain
load('8betagrains_128_3dgrid.mat');

%% == interpolation function of relaxation of misfit strain
phiB2= phiB.^2; sumphiB2= sum(phiB2, 4);
maxsumphiB2= max(sumphiB2,[],'all');
minsumphiB2= min(sumphiB2,[],'all');

psi= (rho* (maxsumphiB2- sumphiB2)+ (sumphiB2- minsumphiB2))/(maxsumphiB2- minsumphiB2);

StrucphiB0= ones(nx,ny,nz,grainBs);
StrucphiB0(phiB==0)= 0;
StrucphiB= zeros(nx,ny,nz,variants,grainBs);
for g= 1: grainBs
    for v= 1: variants

        StrucphiB(:,:,:,v,g)= StrucphiB0(:,:,:,g);

    end
end

phiA= zeros(nx,ny,nz,variants,grainBs);

%% == pre-fft parameters
tmpkx= 2*pi*[0: nx/2 -nx/2+1: -1]/nx; tmpky= tmpkx; tmpkz= tmpkx;
[kx, ky, kz]= ndgrid(tmpkx, tmpky, tmpkz);k2= kx.^2+ ky.^2+ kz.^2; 
kx= kx./k2.^0.5; ky= ky./k2.^0.5; kz= kz./k2.^0.5;
kx(isnan(kx))= 0; ky(isnan(ky))= 0; kz(isnan(kz))= 0;

%% -- green's operator and eigen strain matrix of different variants in different parent beta grains
[tot_Cpq,hom_Cpq,het_Cpq,eigen,gomega,tot_Cpqrs,Spq]= greenmatrix_stiffnessmatrix_eigenstrain(grainBs,variants,phiB,psi,ang,rot0,e0,nx,ny,nz,kx,ky,kz,c11,c12,c44,k_pl);

%% -- denominator in Allen-Cahn function
denom = 1+ dt* M* kappa* k2;

%% == microstructure evolution		
for istep= 1: nstep

    ttime= ttime + dt;

	phiA2= phiA.^2; sumphiA2= sum(phiA2, [4,5]); phiAk= fft(fft(fft(phiA, [], 1), [], 2), [], 3);
    
    if max(sumphiA2,[],'all')> 1

        sumphiA2= sumphiA2/max(sumphiA2,[],'all'); 

    end    
	
	dfdphiA= a* phiA- b* phiA2+ c* phiA.* sumphiA2;
	dfdphiAk= fft(fft(fft(dfdphiA, [], 1), [], 2), [], 3);
    
    [Sigmay2, e0_pl]= yieldstress(sumphiA2,sigmaA_pl,sigmaB_pl,kappaHP_pl,D_alpha_pl,G_pl,delta_Gm,Vmol,JCa_pl,JCb_pl,JCn_pl,JCm1_pl,JCm2_pl,T,T0,Tr,e1_pl,e2_pl,e3_pl,e4_pl,e5_pl,e6_pl);
	
    [deldphiAk,E_elastic,e1_pl,e2_pl,e3_pl,e4_pl,e5_pl,e6_pl,...
     output_s11,output_s22,output_s33,output_s23,output_s13,output_s12]= solve_elasticity(phiA,StrucphiB,eigen,s0yy,e1_pl,e2_pl,e3_pl,e4_pl,e5_pl,e6_pl,grainBs,variants,...
                                                                                          tot_Cpq,hom_Cpq,het_Cpq,gomega,tot_Cpqrs,Spq,Sigmay2,nx,ny,nz,kx,ky,kz,dt,delta_Gm,Vmol);

    output_Vons= (output_s11.^2+ output_s22.^2+ output_s33.^2- (output_s11.* output_s22+ output_s11.* output_s33+ output_s22.* output_s33)+ 3* (output_s12.^2+ output_s13.^2+ output_s23.^2)).^0.5;

    if (istep<= nucleationstep)

        r1= rand(nx,ny,nz,variants,grainBs); r2= rand(nx,ny,nz,variants,grainBs);
        noise= -2* log(r1).* sin(2*pi*r2).*cos(2*pi*r2)*ampnoise;
        noisek= fft(fft(fft(noise, [], 1), [], 2), [], 3);

        phiAk= (phiAk- dt* M.* (dfdphiAk+ deldphiAk))./denom+ noisek;

    else

        phiAk= (phiAk- dt* M.* (dfdphiAk+ deldphiAk))./denom;

    end
	
	phiA= StrucphiB.* real(ifft(ifft(ifft(phiAk, [], 1), [], 2), [], 3));
	
	inrange= (phiA> 1); phiA(inrange)= 1; inrange= (phiA< 0); phiA(inrange)= 0;

	tmpphiA2= phiA.^2;
    for g= 1: grainBs
        for v= 1: variants
            
            tmpphiA2_vg= tmpphiA2(:,:,:,v,g);
            inrange= (tmpphiA2_vg> 0.5); tmpphiAplot(inrange)= 1;
            phiAplot(:,:,:,v,g)= tmpphiAplot;

            sumphiAplot1= sumphiAplot1+ phiAplot(:,:,:,v,g);
            sumphiAplot2= sumphiAplot2+ v* phiAplot(:,:,:,v,g);

            tmpphiAplot= zeros(nx, ny, nz);

        end
    end
	
	if (mod(istep, nprint)== 0)

        filename1= ['phiA2_ttime_',num2str(ttime,'%6.4f'),'s.mat']; save(filename1, 'tmpphiA2');
        filename2= ['sumphiAplot2_ttime_',num2str(ttime,'%6.4f'),'s.mat']; save(filename2, 'sumphiAplot2');
        filename3= ['VonMisesStress_ttime_',num2str(ttime,'%6.4f'),'s.mat']; save(filename3, 'output_Vons');
        filename4= ['EquivalentPlasticStrain_ttime_',num2str(ttime,'%6.4f'),'s.mat']; save(filename4, 'e0_pl');
        filename5= ['ElasticEnergy_ttime_',num2str(ttime,'%6.4f'),'s.mat']; save(filename5, 'E_elastic');

    end
	
	maxphiA2= max(phiA,[],'all'); maxsumphiAplot= max(sumphiAplot1,[],'all'); 	
	phiAplot= zeros(nx, ny, nz, variants, grainBs); 
    sumphiAplot1= zeros(nx, ny, nz); sumphiAplot2= zeros(nx, ny, nz);	
	
	if (maxsumphiAplot> 1| maxphiA2< 2e-2| isinf(deldphiAk))

        fileID= fopen('breakvalue.txt','w');
        fprintf(fileID, 'running step:  %d\r\n', istep);
        fprintf(fileID, 'maxsumphiAplot: %d;   maxphiA2:  %6.4f', maxsumphiAplot, maxphiA2);
		break				
    end   

end

filename= 'VolumnFraction.mat'; save(filename, 'VolF'); 
fclose(fileID);
