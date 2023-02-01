function [phiB]= ParentGrains(nx, ny, grainBs) 
 
dt= 0.1; nstep= 1e4; M= 1; Kappa= 1; nxy= nx* ny;

seed= randi([10,nx-10],grainBs,2); R= 2;
phiB= zeros(nx, ny, grainBs);
for g= 1: grainBs
    phiB(seed(g,1)-R:seed(g,1)+R, seed(g,2)-R:seed(g,2)+R, g)=1;
end

kx= 2*pi*[0: nx/2 -nx/2+1: -1]/nx; ky= kx;
[kx, ky]= ndgrid(kx, ky); k2= kx.^2+ ky.^2;

phiBk= zeros(nx,ny,grainBs); phiplot= zeros(nx,ny);
dfdphiB= zeros(nx,ny,grainBs); dfdphiBk= zeros(nx,ny,grainBs);

cc= hsv(grainBs+1); cc(1,1:3)=[0 0 0];

for istep= 1: nstep

    phiB2 = phiB.^2; sumphiB2= sum(phiB2, 3);
    
    for g= 1: grainBs
		
		phiBk(:,:,g)= fftn(phiB(:,:,g));
        dfdphiB(:,:,g)= -phiB(:,:,g)+ phiB(:,:,g).^3+ 3* phiB(:,:,g).*(sumphiB2- phiB2(:,:,g));
        dfdphiBk(:,:,g)= fftn(dfdphiB(:,:,g));

        phiBk(:,:,g)= (phiBk(:,:,g)- dt* M* dfdphiBk(:,:,g))./(1+ dt* M* Kappa* k2);
        phiB(:,:,g)= real(ifftn(phiBk(:,:,g)));
        
        tmpphi= phiB(:,:,g);
        inrange= (tmpphi>0.9999); tmpphi(inrange)=1;
        inrange= (tmpphi<0.0001); tmpphi(inrange)=0;
        
        phiB(:,:,g)= tmpphi;
    end
    
    totalphiB= sum(phiB, 3);  minphi= min(totalphiB,[],'all');
    
    if (minphi>=0.9995)
	
        break
		
    end

end

phiBlong= reshape(phiB,[],grainBs);
for gp= 1: grainBs
    
    inrange= (phiBlong(:,gp)== 1);
    iflag=(1: grainBs); iflag(iflag== gp)= [];   
    phiBlong(inrange,iflag)= 0;
    
end

totalphiBlong= sum(phiBlong, 2); 
diffphiBlong(:,2)= 1- totalphiBlong; diffphiBlong(:,1)= 1: nxy;
id= diffphiBlong(:,2)== 0| diffphiBlong(:,2)== 1; diffphiBlong(id,:)= [];
len= size(diffphiBlong,1); 
for r= 1: len
    
    loc= diffphiBlong(r,1);
    [~,c]= max(phiBlong(loc,:));
    phiBlong(loc,c)= phiBlong(loc,c)+ diffphiBlong(r,2);
    
end

inrange= (phiBlong> 0.9999); phiBlong(inrange)= 1;
inrange= (phiBlong< 0.0001); phiBlong(inrange)= 0;
phiB = reshape(phiBlong, nx, ny, grainBs);

for g= 1: grainBs
    
    inrange= (phiB(:,:,g)>0.5); phiplot(inrange)= g;
    
end

%figure(10);imagesc(phiplot); grid off; axis equal; axis([1 nx 1 ny]); colormap(cc); colorbar
f= figure('visible','off');imagesc(phiplot); grid off; axis equal; axis([1 nx 1 ny]); colormap(cc); colorbar
filename= sprintf('parent_grains.fig'); savefig(f,filename); f;clf 

filename1= ['BetaMap_16p1028grid.mat']; save(filename1, 'phiB');

end % end function