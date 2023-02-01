function [tot_Cpq,hom_Cpq,het_Cpq,eigen,gomega,tot_Cpqrs,Spq] = greenmatrix_stiffnessmatrix_eigenstrain(grainBs,variants,phiB,psi,ang,e0,nx,ny,kx,ky,c11,c12,c44,k_pl)

format long;

%% == position-dependent elastic constant
Cijkl= zeros(2,2,2,2); 
Cijkl(1,1,1,1)= c11; Cijkl(2,2,2,2)= c11;
Cijkl(1,1,2,2)= c12; Cijkl(2,2,1,1)= c12;

Cijkl(1,2,1,2)= c44; Cijkl(1,2,2,1)= c44;
Cijkl(2,1,1,2)= c44; Cijkl(2,1,2,1)= c44;

rot= zeros(2,2,grainBs);
tot_Cpqrs= zeros(nx,ny,2,2,2,2);
hom_Cpqrs= zeros(2,2,2,2);
het_Cpqrs= zeros(nx,ny,2,2,2,2);

for p= 1: 2
    for q= 1: 2
        for r= 1: 2
            for s= 1: 2 

                for g= 1: grainBs

                    rot(:,:,g)= [cos(ang(g)) sin(ang(g)); -sin(ang(g)) cos(ang(g))];
                    
                    for i= 1: 2
                        for j= 1: 2
                            for k= 1: 2
                                for l= 1: 2                                   

                                    tot_Cpqrs(:,:,p,q,r,s)= tot_Cpqrs(:,:,p,q,r,s)+ phiB(:,:,g)* rot(p,i,g)* rot(q,j,g)* rot(r,k,g)* rot(s,l,g)* Cijkl(i,j,k,l);

                                end
                            end
                        end
                    end

                end

                hom_Cpqrs(p,q,r,s)= (max(tot_Cpqrs(:,:,p,q,r,s),[],'all')+ min(tot_Cpqrs(:,:,p,q,r,s),[],'all'))/2;
                het_Cpqrs(:,:,p,q,r,s)= tot_Cpqrs(:,:,p,q,r,s)- hom_Cpqrs(p,q,r,s);
                
            end
        end
    end
end
    
tot_Cpq(:,:,1,1)= tot_Cpqrs(:,:,1,1,1,1); tot_Cpq(:,:,1,2)= tot_Cpqrs(:,:,1,1,2,2); tot_Cpq(:,:,1,3)= tot_Cpqrs(:,:,1,1,1,2);
tot_Cpq(:,:,2,1)= tot_Cpqrs(:,:,2,2,1,1); tot_Cpq(:,:,2,2)= tot_Cpqrs(:,:,2,2,2,2); tot_Cpq(:,:,2,3)= tot_Cpqrs(:,:,2,2,1,2);
tot_Cpq(:,:,3,1)= tot_Cpqrs(:,:,1,2,1,1); tot_Cpq(:,:,3,2)= tot_Cpqrs(:,:,1,2,2,2); tot_Cpq(:,:,3,3)= tot_Cpqrs(:,:,1,2,1,2);

hom_Cpq(1,1)= hom_Cpqrs(1,1,1,1); hom_Cpq(1,2)= hom_Cpqrs(1,1,2,2); hom_Cpq(1,3)= hom_Cpqrs(1,1,1,2);
hom_Cpq(2,1)= hom_Cpqrs(2,2,1,1); hom_Cpq(2,2)= hom_Cpqrs(2,2,2,2); hom_Cpq(2,3)= hom_Cpqrs(2,2,1,2);
hom_Cpq(3,1)= hom_Cpqrs(1,2,1,1); hom_Cpq(3,2)= hom_Cpqrs(1,2,2,2); hom_Cpq(3,3)= hom_Cpqrs(1,2,1,2);

het_Cpq(:,:,1,1)= het_Cpqrs(:,:,1,1,1,1); het_Cpq(:,:,1,2)= het_Cpqrs(:,:,1,1,2,2); het_Cpq(:,:,1,3)= het_Cpqrs(:,:,1,1,1,2);
het_Cpq(:,:,2,1)= het_Cpqrs(:,:,2,2,1,1); het_Cpq(:,:,2,2)= het_Cpqrs(:,:,2,2,2,2); het_Cpq(:,:,2,3)= het_Cpqrs(:,:,2,2,1,2);
het_Cpq(:,:,3,1)= het_Cpqrs(:,:,1,2,1,1); het_Cpq(:,:,3,2)= het_Cpqrs(:,:,1,2,2,2); het_Cpq(:,:,3,3)= het_Cpqrs(:,:,1,2,1,2);

%% == eigen strain for different variants in different parent grains: eigen
eigen= zeros(nx,ny,2,2,variants,grainBs);

for g= 1: grainBs
    for v= 1: variants
        for ii= 1: 2
            for jj= 1: 2
                for kk= 1: 2
                    for ll= 1: 2

                        eigen(:,:,ii,jj,v,g)= eigen(:,:,ii,jj,v,g)+ psi.* rot(ii,kk,g)* rot(jj,ll,g)* e0(kk,ll,v);

                    end
                end
            end
        end
    end
end

%% == green's tensor
gomega= zeros(nx,ny,3);
% Spqrs= zeros(nx,ny,2,2,2,2);                    % compliance tensor Spqrs= inv(Cpqrs);
Spq= zeros(nx,ny,3,3);

for ix= 1: nx
	for iy= 1: ny

        Spq(ix,iy,:,:)= inv(squeeze(tot_Cpq(ix,iy,:,:)));

		n=[kx(ix,iy), ky(ix,iy)];
		for ii= 1: 2
			for jj= 1: 2
				
				iomega(ii,jj)= hom_Cpqrs(1,ii,jj,1)* n(1)* n(1)+ hom_Cpqrs(1,ii,jj,2)* n(1)* n(2)+... 
							   hom_Cpqrs(2,ii,jj,1)* n(2)* n(1)+ hom_Cpqrs(2,ii,jj,2)* n(2)* n(2);
						   
			end
		end
		omega= inv(iomega);

		gomega(ix,iy,1)= omega(1,1); 
		gomega(ix,iy,2)= omega(2,2); 
		gomega(ix,iy,3)= omega(1,2);  
		
	end
end

gomega(1,1,:)=0;

Spq= Spq/k_pl;

end  % end function      