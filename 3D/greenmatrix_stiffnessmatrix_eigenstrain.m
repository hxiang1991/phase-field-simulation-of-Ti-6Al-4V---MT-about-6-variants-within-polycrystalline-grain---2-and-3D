function [tot_Cpq,hom_Cpq,het_Cpq,eigen,gomega,tot_Cpqrs,Spq] = greenmatrix_stiffnessmatrix_eigenstrain(grainBs,variants,phiB,psi,ang,rot0,e0,nx,ny,nz,kx,ky,kz,c11,c12,c44,k_pl)

format long;

%% == position-dependent elastic constant
Cijkl= zeros(3,3,3,3); 

Cijkl(1,1,1,1)= c11; Cijkl(2,2,2,2)= c11; Cijkl(3,3,3,3)= c11;
Cijkl(1,1,2,2)= c12; Cijkl(2,2,3,3)= c12; Cijkl(3,3,1,1)= c12;
Cijkl(2,2,1,1)= c12; Cijkl(3,3,2,2)= c12; Cijkl(1,1,3,3)= c12;

Cijkl(1,2,1,2)= c44; Cijkl(2,3,2,3)= c44; Cijkl(3,1,3,1)= c44;
Cijkl(2,1,2,1)= c44; Cijkl(3,2,3,2)= c44; Cijkl(1,3,1,3)= c44;
Cijkl(1,2,2,1)= c44; Cijkl(3,2,2,3)= c44; Cijkl(1,3,3,1)= c44;
Cijkl(2,1,1,2)= c44; Cijkl(2,3,3,2)= c44; Cijkl(3,1,1,3)= c44;

rot_ang= zeros(3,3,grainBs); rot= zeros(3,3,grainBs);
tot_Cpqrs= zeros(nx,ny,nz,3,3,3,3);
hom_Cpqrs= zeros(3,3,3,3);
het_Cpqrs= zeros(nx,ny,nz,3,3,3,3);

for p= 1: 3
    for q= 1: 3
        for r= 1: 3
            for s= 1: 3 

                for g= 1: grainBs

                    sa= sin(ang(g,1)); sb= sin(ang(g,2)); sg= sin(ang(g,3));   % sa: sin(alpha); sb: sin(beta); sg: sin(gamma);
					ca= cos(ang(g,1)); cb= cos(ang(g,2)); cg= cos(ang(g,3));   % ca: cos(alpha); cb: cos(beta); cg: cos(gamma);
					
					rot_ang(:,:,g)= [ ca* cb- sa* sb* cg,  sa* cb+ ca* sb* cg, sb* sg;
						             -ca* sb- sa* cb* cg, -sa* sb+ ca* cb* cg, cb* sg;
						                          sa* sg,             -ca* sg,     cg];

                    rot(:,:,g)= rot0/rot_ang(:,:,g);
                    
                    for i= 1: 3
                        for j= 1: 3
                            for k= 1: 3
                                for l= 1: 3                                   

                                    tot_Cpqrs(:,:,:,p,q,r,s)= tot_Cpqrs(:,:,:,p,q,r,s)+ phiB(:,:,:,g)* rot(p,i,g)* rot(q,j,g)* rot(r,k,g)* rot(s,l,g)* Cijkl(i,j,k,l);

                                end
                            end
                        end
                    end

                end

                hom_Cpqrs(p,q,r,s)= (max(tot_Cpqrs(:,:,:,p,q,r,s),[],'all')+ min(tot_Cpqrs(:,:,:,p,q,r,s),[],'all'))/2;
                het_Cpqrs(:,:,:,p,q,r,s)= tot_Cpqrs(:,:,:,p,q,r,s)- hom_Cpqrs(p,q,r,s);
                
            end
        end
    end
end
    
tot_Cpq(:,:,:,1,1)= tot_Cpqrs(:,:,:,1,1,1,1); tot_Cpq(:,:,:,1,2)= tot_Cpqrs(:,:,:,1,1,2,2); tot_Cpq(:,:,:,1,3)= tot_Cpqrs(:,:,:,1,1,3,3);
tot_Cpq(:,:,:,1,4)= tot_Cpqrs(:,:,:,1,1,2,3); tot_Cpq(:,:,:,1,5)= tot_Cpqrs(:,:,:,1,1,1,3); tot_Cpq(:,:,:,1,6)= tot_Cpqrs(:,:,:,1,1,1,2);

tot_Cpq(:,:,:,2,1)= tot_Cpqrs(:,:,:,2,2,1,1); tot_Cpq(:,:,:,2,2)= tot_Cpqrs(:,:,:,2,2,2,2); tot_Cpq(:,:,:,2,3)= tot_Cpqrs(:,:,:,2,2,3,3);
tot_Cpq(:,:,:,2,4)= tot_Cpqrs(:,:,:,2,2,2,3); tot_Cpq(:,:,:,2,5)= tot_Cpqrs(:,:,:,2,2,1,3); tot_Cpq(:,:,:,2,6)= tot_Cpqrs(:,:,:,2,2,1,2);

tot_Cpq(:,:,:,3,1)= tot_Cpqrs(:,:,:,3,3,1,1); tot_Cpq(:,:,:,3,2)= tot_Cpqrs(:,:,:,3,3,2,2); tot_Cpq(:,:,:,3,3)= tot_Cpqrs(:,:,:,3,3,3,3);
tot_Cpq(:,:,:,3,4)= tot_Cpqrs(:,:,:,3,3,2,3); tot_Cpq(:,:,:,3,5)= tot_Cpqrs(:,:,:,3,3,1,3); tot_Cpq(:,:,:,3,6)= tot_Cpqrs(:,:,:,3,3,1,2);

tot_Cpq(:,:,:,4,1)= tot_Cpqrs(:,:,:,2,3,1,1); tot_Cpq(:,:,:,4,2)= tot_Cpqrs(:,:,:,2,3,2,2); tot_Cpq(:,:,:,4,3)= tot_Cpqrs(:,:,:,2,3,3,3);
tot_Cpq(:,:,:,4,4)= tot_Cpqrs(:,:,:,2,3,2,3); tot_Cpq(:,:,:,4,5)= tot_Cpqrs(:,:,:,2,3,1,3); tot_Cpq(:,:,:,4,6)= tot_Cpqrs(:,:,:,2,3,1,2);

tot_Cpq(:,:,:,5,1)= tot_Cpqrs(:,:,:,1,3,1,1); tot_Cpq(:,:,:,5,2)= tot_Cpqrs(:,:,:,1,3,2,2); tot_Cpq(:,:,:,5,3)= tot_Cpqrs(:,:,:,1,3,3,3);
tot_Cpq(:,:,:,5,4)= tot_Cpqrs(:,:,:,1,3,2,3); tot_Cpq(:,:,:,5,5)= tot_Cpqrs(:,:,:,1,3,1,3); tot_Cpq(:,:,:,5,6)= tot_Cpqrs(:,:,:,1,3,1,2);

tot_Cpq(:,:,:,6,1)= tot_Cpqrs(:,:,:,1,2,1,1); tot_Cpq(:,:,:,6,2)= tot_Cpqrs(:,:,:,1,2,2,2); tot_Cpq(:,:,:,6,3)= tot_Cpqrs(:,:,:,1,2,3,3);
tot_Cpq(:,:,:,6,4)= tot_Cpqrs(:,:,:,1,2,2,3); tot_Cpq(:,:,:,6,5)= tot_Cpqrs(:,:,:,1,2,1,3); tot_Cpq(:,:,:,6,6)= tot_Cpqrs(:,:,:,1,2,1,2);

hom_Cpq(1,1)= hom_Cpqrs(1,1,1,1); hom_Cpq(1,2)= hom_Cpqrs(1,1,2,2); hom_Cpq(1,3)= hom_Cpqrs(1,1,3,3);
hom_Cpq(1,4)= hom_Cpqrs(1,1,2,3); hom_Cpq(1,5)= hom_Cpqrs(1,1,1,3); hom_Cpq(1,6)= hom_Cpqrs(1,1,1,2);

hom_Cpq(2,1)= hom_Cpqrs(2,2,1,1); hom_Cpq(2,2)= hom_Cpqrs(2,2,2,2); hom_Cpq(2,3)= hom_Cpqrs(2,2,3,3);
hom_Cpq(2,4)= hom_Cpqrs(2,2,2,3); hom_Cpq(2,5)= hom_Cpqrs(2,2,1,3); hom_Cpq(2,6)= hom_Cpqrs(2,2,1,2);

hom_Cpq(3,1)= hom_Cpqrs(3,3,1,1); hom_Cpq(3,2)= hom_Cpqrs(3,3,2,2); hom_Cpq(3,3)= hom_Cpqrs(3,3,3,3);
hom_Cpq(3,4)= hom_Cpqrs(3,3,2,3); hom_Cpq(3,5)= hom_Cpqrs(3,3,1,3); hom_Cpq(3,6)= hom_Cpqrs(3,3,1,2);

hom_Cpq(4,1)= hom_Cpqrs(2,3,1,1); hom_Cpq(4,2)= hom_Cpqrs(2,3,2,2); hom_Cpq(4,3)= hom_Cpqrs(2,3,3,3);
hom_Cpq(4,4)= hom_Cpqrs(2,3,2,3); hom_Cpq(4,5)= hom_Cpqrs(2,3,1,3); hom_Cpq(4,6)= hom_Cpqrs(2,3,1,2);

hom_Cpq(5,1)= hom_Cpqrs(1,3,1,1); hom_Cpq(5,2)= hom_Cpqrs(1,3,2,2); hom_Cpq(5,3)= hom_Cpqrs(1,3,3,3);
hom_Cpq(5,4)= hom_Cpqrs(1,3,2,3); hom_Cpq(5,5)= hom_Cpqrs(1,3,1,3); hom_Cpq(5,6)= hom_Cpqrs(1,3,1,2);

hom_Cpq(6,1)= hom_Cpqrs(1,2,1,1); hom_Cpq(6,2)= hom_Cpqrs(1,2,2,2); hom_Cpq(6,3)= hom_Cpqrs(1,2,3,3);
hom_Cpq(6,4)= hom_Cpqrs(1,2,2,3); hom_Cpq(6,5)= hom_Cpqrs(1,2,1,3); hom_Cpq(6,6)= hom_Cpqrs(1,2,1,2);

het_Cpq(:,:,:,1,1)= het_Cpqrs(:,:,:,1,1,1,1); het_Cpq(:,:,:,1,2)= het_Cpqrs(:,:,:,1,1,2,2); het_Cpq(:,:,:,1,3)= het_Cpqrs(:,:,:,1,1,3,3);
het_Cpq(:,:,:,1,4)= het_Cpqrs(:,:,:,1,1,2,3); het_Cpq(:,:,:,1,5)= het_Cpqrs(:,:,:,1,1,1,3); het_Cpq(:,:,:,1,6)= het_Cpqrs(:,:,:,1,1,1,2);

het_Cpq(:,:,:,2,1)= het_Cpqrs(:,:,:,2,2,1,1); het_Cpq(:,:,:,2,2)= het_Cpqrs(:,:,:,2,2,2,2); het_Cpq(:,:,:,2,3)= het_Cpqrs(:,:,:,2,2,3,3);
het_Cpq(:,:,:,2,4)= het_Cpqrs(:,:,:,2,2,2,3); het_Cpq(:,:,:,2,5)= het_Cpqrs(:,:,:,2,2,1,3); het_Cpq(:,:,:,2,6)= het_Cpqrs(:,:,:,2,2,1,2);

het_Cpq(:,:,:,3,1)= het_Cpqrs(:,:,:,3,3,1,1); het_Cpq(:,:,:,3,2)= het_Cpqrs(:,:,:,3,3,2,2); het_Cpq(:,:,:,3,3)= het_Cpqrs(:,:,:,3,3,3,3);
het_Cpq(:,:,:,3,4)= het_Cpqrs(:,:,:,3,3,2,3); het_Cpq(:,:,:,3,5)= het_Cpqrs(:,:,:,3,3,1,3); het_Cpq(:,:,:,3,6)= het_Cpqrs(:,:,:,3,3,1,2);

het_Cpq(:,:,:,4,1)= het_Cpqrs(:,:,:,2,3,1,1); het_Cpq(:,:,:,4,2)= het_Cpqrs(:,:,:,2,3,2,2); het_Cpq(:,:,:,4,3)= het_Cpqrs(:,:,:,2,3,3,3);
het_Cpq(:,:,:,4,4)= het_Cpqrs(:,:,:,2,3,2,3); het_Cpq(:,:,:,4,5)= het_Cpqrs(:,:,:,2,3,1,3); het_Cpq(:,:,:,4,6)= het_Cpqrs(:,:,:,2,3,1,2);

het_Cpq(:,:,:,5,1)= het_Cpqrs(:,:,:,1,3,1,1); het_Cpq(:,:,:,5,2)= het_Cpqrs(:,:,:,1,3,2,2); het_Cpq(:,:,:,5,3)= het_Cpqrs(:,:,:,1,3,3,3);
het_Cpq(:,:,:,5,4)= het_Cpqrs(:,:,:,1,3,2,3); het_Cpq(:,:,:,5,5)= het_Cpqrs(:,:,:,1,3,1,3); het_Cpq(:,:,:,5,6)= het_Cpqrs(:,:,:,1,3,1,2);

het_Cpq(:,:,:,6,1)= het_Cpqrs(:,:,:,1,2,1,1); het_Cpq(:,:,:,6,2)= het_Cpqrs(:,:,:,1,2,2,2); het_Cpq(:,:,:,6,3)= het_Cpqrs(:,:,:,1,2,3,3);
het_Cpq(:,:,:,6,4)= het_Cpqrs(:,:,:,1,2,2,3); het_Cpq(:,:,:,6,5)= het_Cpqrs(:,:,:,1,2,1,3); het_Cpq(:,:,:,6,6)= het_Cpqrs(:,:,:,1,2,1,2);

%% == eigen strain for different variants in different parent grains: eigen
eigen= zeros(nx,ny,nz,3,3,variants,grainBs);

for g= 1: grainBs
    for v= 1: variants
        for ii= 1: 3
            for jj= 1: 3
                for kk= 1: 3
                    for ll= 1: 3

                        eigen(:,:,:,ii,jj,v,g)= eigen(:,:,:,ii,jj,v,g)+ psi.* rot(ii,kk,g)* rot(jj,ll,g)* e0(kk,ll,v);

                    end
                end
            end
        end
    end
end

%% == green's tensor
gomega= zeros(nx,ny,nz,6);
Spq= zeros(nx,ny,nz,6,6);

warning('off');
for ix= 1: nx
	for iy= 1: ny
        for iz= 1: nz

            Spq(ix,iy,iz,:,:)= inv(squeeze(tot_Cpq(ix,iy,iz,:,:))); 
    
		    n=[kx(ix,iy,iz), ky(ix,iy,iz), kz(ix,iy,iz)];
		    for ii= 1: 3
			    for jj= 1: 3
				    
				    iomega(ii,jj)= hom_Cpqrs(1,ii,jj,1)* n(1)* n(1)+ hom_Cpqrs(1,ii,jj,2)* n(1)* n(2)+ hom_Cpqrs(1,ii,jj,3)* n(1)* n(3)+... 
							       hom_Cpqrs(2,ii,jj,1)* n(2)* n(1)+ hom_Cpqrs(2,ii,jj,2)* n(2)* n(2)+ hom_Cpqrs(2,ii,jj,3)* n(2)* n(3)+...
                                   hom_Cpqrs(3,ii,jj,1)* n(3)* n(1)+ hom_Cpqrs(3,ii,jj,2)* n(3)* n(2)+ hom_Cpqrs(3,ii,jj,3)* n(3)* n(3);
						       
			    end
		    end
		    omega= inv(iomega);
    
		    gomega(ix,iy,iz,1)= omega(1,1); gomega(ix,iy,iz,2)= omega(2,2); gomega(ix,iy,iz,3)= omega(3,3); 
            gomega(ix,iy,iz,4)= omega(2,3); gomega(ix,iy,iz,5)= omega(1,3); gomega(ix,iy,iz,6)= omega(1,2);  

        end		    
	end
end

gomega(1,1,1,:)=0;

Spq= Spq/k_pl;

warning('on');

end  % end function      