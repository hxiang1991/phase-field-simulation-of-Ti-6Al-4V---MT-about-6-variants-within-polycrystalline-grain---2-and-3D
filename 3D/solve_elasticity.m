function [deldphiAk,E_elastic,e1_pl,e2_pl,e3_pl,e4_pl,e5_pl,e6_pl,...
          output_s11,output_s22,output_s33,output_s23,output_s13,output_s12]= solve_elasticity(phiA,StrucphiB,eigen,s0yy,e1_pl,e2_pl,e3_pl,e4_pl,e5_pl,e6_pl,grainBs,variants,...
                                                                                               tot_Cpq,hom_Cpq,het_Cpq,gomega,tot_Cpqrs,Spq,Sigmay2,nx,ny,nz,kx,ky,kz,dt,delta_Gm,Vmol)

format long;

niter= 1e2; tolerance= 3e-4;

%% == eigenstrains
ei1= zeros(nx,ny,nz); ei2= zeros(nx,ny,nz); ei3= zeros(nx,ny,nz); 
ei4= zeros(nx,ny,nz); ei5= zeros(nx,ny,nz); ei6= zeros(nx,ny,nz);  

for g= 1: grainBs
    for v= 1: variants
        
        ei1= ei1+ StrucphiB(:,:,:,v,g).* eigen(:,:,:,1,1,v,g).* phiA(:,:,:,v,g); 
        ei2= ei2+ StrucphiB(:,:,:,v,g).* eigen(:,:,:,2,2,v,g).* phiA(:,:,:,v,g); 
        ei3= ei3+ StrucphiB(:,:,:,v,g).* eigen(:,:,:,3,3,v,g).* phiA(:,:,:,v,g);
		ei4= ei4+ StrucphiB(:,:,:,v,g).* eigen(:,:,:,2,3,v,g).* phiA(:,:,:,v,g);
        ei5= ei5+ StrucphiB(:,:,:,v,g).* eigen(:,:,:,1,3,v,g).* phiA(:,:,:,v,g); 
        ei6= ei6+ StrucphiB(:,:,:,v,g).* eigen(:,:,:,1,2,v,g).* phiA(:,:,:,v,g);
        
    end
end

ei1= ei1+ e1_pl; ei2= ei2+ e2_pl; ei3= ei3+ e3_pl;
ei4= ei4+ e4_pl; ei5= ei5+ e5_pl; ei6= ei6+ e6_pl;

ei_devia= (ei1+ ei2+ ei3)/3;

%% == elastic strain and stress field
u= zeros(nx,ny,nz,3,niter); uk= zeros(nx,ny,nz,3); 
T0= zeros(nx,ny,nz,6);
T= zeros(nx,ny,nz,6); Tk= zeros(nx,ny,nz,6);
hom_e= zeros(6,1);
het_e= zeros(nx,ny,nz,6); het_ek= zeros(nx,ny,nz,6);
s= zeros(nx,ny,nz,6); sk= zeros(nx,ny,nz,6);

% == zeroth-order iteration
s(:,:,:,1)= hom_Cpq(1,1)* ei1+ hom_Cpq(1,2)* ei2+ hom_Cpq(1,3)* ei3+ 2* (hom_Cpq(1,4)* ei4+ hom_Cpq(1,5)* ei5+ hom_Cpq(1,6)* ei6); sk(:,:,:,1)= fftn(s(:,:,:,1));
s(:,:,:,2)= hom_Cpq(2,1)* ei1+ hom_Cpq(2,2)* ei2+ hom_Cpq(2,3)* ei3+ 2* (hom_Cpq(2,4)* ei4+ hom_Cpq(2,5)* ei5+ hom_Cpq(2,6)* ei6); sk(:,:,:,2)= fftn(s(:,:,:,2));
s(:,:,:,3)= hom_Cpq(3,1)* ei1+ hom_Cpq(3,2)* ei2+ hom_Cpq(3,3)* ei3+ 2* (hom_Cpq(3,4)* ei4+ hom_Cpq(3,5)* ei5+ hom_Cpq(3,6)* ei6); sk(:,:,:,3)= fftn(s(:,:,:,3));
s(:,:,:,4)= hom_Cpq(4,1)* ei1+ hom_Cpq(4,2)* ei2+ hom_Cpq(4,3)* ei3+ 2* (hom_Cpq(4,4)* ei4+ hom_Cpq(4,5)* ei5+ hom_Cpq(4,6)* ei6); sk(:,:,:,4)= fftn(s(:,:,:,4));
s(:,:,:,5)= hom_Cpq(5,1)* ei1+ hom_Cpq(5,2)* ei2+ hom_Cpq(5,3)* ei3+ 2* (hom_Cpq(5,4)* ei4+ hom_Cpq(5,5)* ei5+ hom_Cpq(5,6)* ei6); sk(:,:,:,5)= fftn(s(:,:,:,5));
s(:,:,:,6)= hom_Cpq(6,1)* ei1+ hom_Cpq(6,2)* ei2+ hom_Cpq(6,3)* ei3+ 2* (hom_Cpq(6,4)* ei4+ hom_Cpq(6,5)* ei5+ hom_Cpq(6,6)* ei6); sk(:,:,:,6)= fftn(s(:,:,:,6));

uk(:,:,:,1)= -1i*(gomega(:,:,:,1).*(sk(:,:,:,1).*kx+ sk(:,:,:,6).*ky+ sk(:,:,:,5).*kz)+...
				  gomega(:,:,:,6).*(sk(:,:,:,6).*kx+ sk(:,:,:,2).*ky+ sk(:,:,:,4).*kz)+...
				  gomega(:,:,:,5).*(sk(:,:,:,5).*kx+ sk(:,:,:,4).*ky+ sk(:,:,:,3).*kz)); 
uk(:,:,:,2)= -1i*(gomega(:,:,:,6).*(sk(:,:,:,1).*kx+ sk(:,:,:,6).*ky+ sk(:,:,:,5).*kz)+...
				  gomega(:,:,:,2).*(sk(:,:,:,6).*kx+ sk(:,:,:,2).*ky+ sk(:,:,:,4).*kz)+...
				  gomega(:,:,:,4).*(sk(:,:,:,5).*kx+ sk(:,:,:,4).*ky+ sk(:,:,:,3).*kz));
uk(:,:,:,3)= -1i*(gomega(:,:,:,5).*(sk(:,:,:,1).*kx+ sk(:,:,:,6).*ky+ sk(:,:,:,5).*kz)+...
				  gomega(:,:,:,4).*(sk(:,:,:,6).*kx+ sk(:,:,:,2).*ky+ sk(:,:,:,4).*kz)+...
				  gomega(:,:,:,3).*(sk(:,:,:,5).*kx+ sk(:,:,:,4).*ky+ sk(:,:,:,3).*kz));	
					
u(:,:,:,1,1)= real(ifftn(uk(:,:,:,1))); u(:,:,:,2,1)= real(ifftn(uk(:,:,:,2))); u(:,:,:,3,1)= real(ifftn(uk(:,:,:,3)));					

het_ek(:,:,:,1)= 1i* kx.* uk(:,:,:,1); het_ek(:,:,:,2)= 1i* ky.* uk(:,:,:,2); het_ek(:,:,:,3)= 1i* kz.* uk(:,:,:,3); 
het_ek(:,:,:,4)= 0.5*1i* (ky.* uk(:,:,:,3)+ kz.* uk(:,:,:,2)); 
het_ek(:,:,:,5)= 0.5*1i* (kx.* uk(:,:,:,3)+ kz.* uk(:,:,:,1));                       
het_ek(:,:,:,6)= 0.5*1i* (kx.* uk(:,:,:,2)+ ky.* uk(:,:,:,1)); 

het_e(:,:,:,1)= real(ifftn(het_ek(:,:,:,1))); het_e(:,:,:,2)= real(ifftn(het_ek(:,:,:,2))); het_e(:,:,:,3)= real(ifftn(het_ek(:,:,:,3)));
het_e(:,:,:,4)= real(ifftn(het_ek(:,:,:,4))); het_e(:,:,:,5)= real(ifftn(het_ek(:,:,:,5))); het_e(:,:,:,6)= real(ifftn(het_ek(:,:,:,6)));

T0(:,:,:,1)= tot_Cpq(:,:,:,1,1).*(ei1- hom_e(1))+ tot_Cpq(:,:,:,1,2).*(ei2- hom_e(2))+ tot_Cpq(:,:,:,1,3).*(ei3- hom_e(3))+...
		 2* (tot_Cpq(:,:,:,1,4).*(ei4- hom_e(4))+ tot_Cpq(:,:,:,1,5).*(ei5- hom_e(5))+ tot_Cpq(:,:,:,1,6).*(ei6- hom_e(6)));
T0(:,:,:,2)= tot_Cpq(:,:,:,2,1).*(ei1- hom_e(1))+ tot_Cpq(:,:,:,2,2).*(ei2- hom_e(2))+ tot_Cpq(:,:,:,2,3).*(ei3- hom_e(3))+...
		 2* (tot_Cpq(:,:,:,2,4).*(ei4- hom_e(4))+ tot_Cpq(:,:,:,2,5).*(ei5- hom_e(5))+ tot_Cpq(:,:,:,2,6).*(ei6- hom_e(6)));
T0(:,:,:,3)= tot_Cpq(:,:,:,3,1).*(ei1- hom_e(1))+ tot_Cpq(:,:,:,3,2).*(ei2- hom_e(2))+ tot_Cpq(:,:,:,3,3).*(ei3- hom_e(3))+...
		 2* (tot_Cpq(:,:,:,3,4).*(ei4- hom_e(4))+ tot_Cpq(:,:,:,3,5).*(ei5- hom_e(5))+ tot_Cpq(:,:,:,3,6).*(ei6- hom_e(6)));
T0(:,:,:,4)= tot_Cpq(:,:,:,4,1).*(ei1- hom_e(1))+ tot_Cpq(:,:,:,4,2).*(ei2- hom_e(2))+ tot_Cpq(:,:,:,4,3).*(ei3- hom_e(3))+...
		 2* (tot_Cpq(:,:,:,4,4).*(ei4- hom_e(4))+ tot_Cpq(:,:,:,4,5).*(ei5- hom_e(5))+ tot_Cpq(:,:,:,4,6).*(ei6- hom_e(6)));
T0(:,:,:,5)= tot_Cpq(:,:,:,5,1).*(ei1- hom_e(1))+ tot_Cpq(:,:,:,5,2).*(ei2- hom_e(2))+ tot_Cpq(:,:,:,5,3).*(ei3- hom_e(3))+...
		 2* (tot_Cpq(:,:,:,5,4).*(ei4- hom_e(4))+ tot_Cpq(:,:,:,5,5).*(ei5- hom_e(5))+ tot_Cpq(:,:,:,5,6).*(ei6- hom_e(6)));
T0(:,:,:,6)= tot_Cpq(:,:,:,6,1).*(ei1- hom_e(1))+ tot_Cpq(:,:,:,6,2).*(ei2- hom_e(2))+ tot_Cpq(:,:,:,6,3).*(ei3- hom_e(3))+...
		 2* (tot_Cpq(:,:,:,6,4).*(ei4- hom_e(4))+ tot_Cpq(:,:,:,6,5).*(ei5- hom_e(5))+ tot_Cpq(:,:,:,6,6).*(ei6- hom_e(6)));

for iter= 1: niter
    
    T(:,:,:,1)= T0(:,:,:,1)- het_Cpq(:,:,:,1,1).* het_e(:,:,:,1)- het_Cpq(:,:,:,1,2).* het_e(:,:,:,2)- het_Cpq(:,:,:,1,3).* het_e(:,:,:,3)-...
	                      2*(het_Cpq(:,:,:,1,4).* het_e(:,:,:,4)+ het_Cpq(:,:,:,1,5).* het_e(:,:,:,5)+ het_Cpq(:,:,:,1,6).* het_e(:,:,:,6));
							   
    T(:,:,:,2)= T0(:,:,:,2)- het_Cpq(:,:,:,2,1).* het_e(:,:,:,1)- het_Cpq(:,:,:,2,2).* het_e(:,:,:,2)- het_Cpq(:,:,:,2,3).* het_e(:,:,:,3)-...
	                      2*(het_Cpq(:,:,:,2,4).* het_e(:,:,:,4)+ het_Cpq(:,:,:,2,5).* het_e(:,:,:,5)+ het_Cpq(:,:,:,2,6).* het_e(:,:,:,6));
    
	T(:,:,:,3)= T0(:,:,:,3)- het_Cpq(:,:,:,3,1).* het_e(:,:,:,1)- het_Cpq(:,:,:,3,2).* het_e(:,:,:,2)- het_Cpq(:,:,:,3,3).* het_e(:,:,:,3)-...
	                      2*(het_Cpq(:,:,:,3,4).* het_e(:,:,:,4)+ het_Cpq(:,:,:,3,5).* het_e(:,:,:,5)+ het_Cpq(:,:,:,3,6).* het_e(:,:,:,6));
							   
    T(:,:,:,4)= T0(:,:,:,4)- het_Cpq(:,:,:,4,1).* het_e(:,:,:,1)- het_Cpq(:,:,:,4,2).* het_e(:,:,:,2)- het_Cpq(:,:,:,4,3).* het_e(:,:,:,3)-...
	                      2*(het_Cpq(:,:,:,4,4).* het_e(:,:,:,4)+ het_Cpq(:,:,:,4,5).* het_e(:,:,:,5)+ het_Cpq(:,:,:,4,6).* het_e(:,:,:,6));
							   
    T(:,:,:,5)= T0(:,:,:,5)- het_Cpq(:,:,:,5,1).* het_e(:,:,:,1)- het_Cpq(:,:,:,5,2).* het_e(:,:,:,2)- het_Cpq(:,:,:,5,3).* het_e(:,:,:,3)-...
	                      2*(het_Cpq(:,:,:,5,4).* het_e(:,:,:,4)+ het_Cpq(:,:,:,5,5).* het_e(:,:,:,5)+ het_Cpq(:,:,:,5,6).* het_e(:,:,:,6));
    
	T(:,:,:,6)= T0(:,:,:,6)- het_Cpq(:,:,:,6,1).* het_e(:,:,:,1)- het_Cpq(:,:,:,6,2).* het_e(:,:,:,2)- het_Cpq(:,:,:,6,3).* het_e(:,:,:,3)-...
	                      2*(het_Cpq(:,:,:,6,4).* het_e(:,:,:,4)+ het_Cpq(:,:,:,6,5).* het_e(:,:,:,5)+ het_Cpq(:,:,:,6,6).* het_e(:,:,:,6));							   
              
    Tk(:,:,:,1)= fftn(T(:,:,:,1)); Tk(:,:,:,2)= fftn(T(:,:,:,2)); Tk(:,:,:,3)= fftn(T(:,:,:,3));
	Tk(:,:,:,4)= fftn(T(:,:,:,4)); Tk(:,:,:,5)= fftn(T(:,:,:,5)); Tk(:,:,:,6)= fftn(T(:,:,:,6));
	
	uk(:,:,:,1)= -1i*(gomega(:,:,:,1).*(Tk(:,:,:,1).*kx+ Tk(:,:,:,6).*ky+ Tk(:,:,:,5).*kz)+...
					  gomega(:,:,:,6).*(Tk(:,:,:,6).*kx+ Tk(:,:,:,2).*ky+ Tk(:,:,:,4).*kz)+...
					  gomega(:,:,:,5).*(Tk(:,:,:,5).*kx+ Tk(:,:,:,4).*ky+ Tk(:,:,:,3).*kz)); 
	uk(:,:,:,2)= -1i*(gomega(:,:,:,6).*(Tk(:,:,:,1).*kx+ Tk(:,:,:,6).*ky+ Tk(:,:,:,5).*kz)+...
					  gomega(:,:,:,2).*(Tk(:,:,:,6).*kx+ Tk(:,:,:,2).*ky+ Tk(:,:,:,4).*kz)+...
				      gomega(:,:,:,4).*(Tk(:,:,:,5).*kx+ Tk(:,:,:,4).*ky+ Tk(:,:,:,3).*kz));
	uk(:,:,:,3)= -1i*(gomega(:,:,:,5).*(Tk(:,:,:,1).*kx+ Tk(:,:,:,6).*ky+ Tk(:,:,:,5).*kz)+...
					  gomega(:,:,:,4).*(Tk(:,:,:,6).*kx+ Tk(:,:,:,2).*ky+ Tk(:,:,:,4).*kz)+...
					  gomega(:,:,:,3).*(Tk(:,:,:,5).*kx+ Tk(:,:,:,4).*ky+ Tk(:,:,:,3).*kz));	
    
    u(:,:,:,1,iter+1)= real(ifftn(uk(:,:,:,1))); u(:,:,:,2,iter+1)= real(ifftn(uk(:,:,:,2))); u(:,:,:,3,iter+1)= real(ifftn(uk(:,:,:,3)));

	het_ek(:,:,:,1)= 1i* kx.* uk(:,:,:,1); het_ek(:,:,:,2)= 1i* ky.* uk(:,:,:,2); het_ek(:,:,:,3)= 1i* kz.* uk(:,:,:,3); 
	het_ek(:,:,:,4)= 0.5*1i* (ky.* uk(:,:,:,3)+ kz.* uk(:,:,:,2)); 
	het_ek(:,:,:,5)= 0.5*1i* (kx.* uk(:,:,:,3)+ kz.* uk(:,:,:,1));                       
	het_ek(:,:,:,6)= 0.5*1i* (kx.* uk(:,:,:,2)+ ky.* uk(:,:,:,1)); 

	het_e(:,:,:,1)= real(ifftn(het_ek(:,:,:,1))); het_e(:,:,:,2)= real(ifftn(het_ek(:,:,:,2))); het_e(:,:,:,3)= real(ifftn(het_ek(:,:,:,3)));
	het_e(:,:,:,4)= real(ifftn(het_ek(:,:,:,4))); het_e(:,:,:,5)= real(ifftn(het_ek(:,:,:,5))); het_e(:,:,:,6)= real(ifftn(het_ek(:,:,:,6)));

    % == check convergence
    diffu2= (u(:,:,:,1,iter+1)- u(:,:,:,1,iter)).^2+ (u(:,:,:,2,iter+1)- u(:,:,:,2,iter)).^2+ (u(:,:,:,3,iter+1)- u(:,:,:,3,iter)).^2;
    conver= trapz((trapz(trapz(diffu2))))^0.5;
    
    if(conver< tolerance|| conver> 1e2)
       break
    end
        
end  %iter

e11= hom_e(1)+ het_e(:,:,:,1)- ei1; e22= hom_e(2)+ het_e(:,:,:,2)- ei2; e33= hom_e(3)+ het_e(:,:,:,3)- ei3;
e23= hom_e(4)+ het_e(:,:,:,4)- ei4; e13= hom_e(5)+ het_e(:,:,:,5)- ei5; e12= hom_e(6)+ het_e(:,:,:,6)- ei6;

e_devia= (het_e(:,:,:,1)+ het_e(:,:,:,2)+ het_e(:,:,:,3))/3;

% ==  stress field and strain energy
s11= tot_Cpq(:,:,:,1,1).* e11+ tot_Cpq(:,:,:,1,2).* e22+ tot_Cpq(:,:,:,1,3).* e33+...
 2* (tot_Cpq(:,:,:,1,4).* e23+ tot_Cpq(:,:,:,1,5).* e13+ tot_Cpq(:,:,:,1,6).* e12);   output_s11= s11* delta_Gm/Vmol/1e06;   % unit :MPa
s22= tot_Cpq(:,:,:,2,1).* e11+ tot_Cpq(:,:,:,2,2).* e22+ tot_Cpq(:,:,:,2,3).* e33+...
 2* (tot_Cpq(:,:,:,2,4).* e23+ tot_Cpq(:,:,:,2,5).* e13+ tot_Cpq(:,:,:,2,6).* e12);   output_s22= s22* delta_Gm/Vmol/1e06;   % unit :MPa
s33= tot_Cpq(:,:,:,3,1).* e11+ tot_Cpq(:,:,:,3,2).* e22+ tot_Cpq(:,:,:,3,3).* e33+...
 2* (tot_Cpq(:,:,:,3,4).* e23+ tot_Cpq(:,:,:,3,5).* e13+ tot_Cpq(:,:,:,3,6).* e12);   output_s33= s33* delta_Gm/Vmol/1e06;   % unit :MPa
s23= tot_Cpq(:,:,:,4,1).* e11+ tot_Cpq(:,:,:,4,2).* e22+ tot_Cpq(:,:,:,4,3).* e33+...
 2* (tot_Cpq(:,:,:,4,4).* e23+ tot_Cpq(:,:,:,4,5).* e13+ tot_Cpq(:,:,:,4,6).* e12);   output_s23= s23* delta_Gm/Vmol/1e06;   % unit :MPa
s13= tot_Cpq(:,:,:,5,1).* e11+ tot_Cpq(:,:,:,5,2).* e22+ tot_Cpq(:,:,:,5,3).* e33+...
 2* (tot_Cpq(:,:,:,5,4).* e23+ tot_Cpq(:,:,:,5,5).* e13+ tot_Cpq(:,:,:,5,6).* e12);   output_s13= s13* delta_Gm/Vmol/1e06;   % unit :MPa
s12= tot_Cpq(:,:,:,6,1).* e11+ tot_Cpq(:,:,:,6,2).* e22+ tot_Cpq(:,:,:,6,3).* e33+...
 2* (tot_Cpq(:,:,:,6,4).* e23+ tot_Cpq(:,:,:,6,5).* e13+ tot_Cpq(:,:,:,6,6).* e12);   output_s12= s12* delta_Gm/Vmol/1e06;   % unit :MPa

deldphiA = zeros(nx, ny, nz, variants, grainBs);
deldphiAk= zeros(nx, ny, nz, variants, grainBs);

E_elastic= (s11.* e11+ s22.* e22+ s33.* e33+ 2* (s23.* e23+ s13.* e13+ s12.* e12))* delta_Gm/Vmol; 

for g= 1: grainBs
    for v= 1: variants
        
        deldphiA(:,:,:,v,g)= -StrucphiB(:,:,:,v,g).*(s11.* eigen(:,:,:,1,1,v,g)+ s22.* eigen(:,:,:,2,2,v,g)+ s33.* eigen(:,:,:,3,3,v,g)+...
                                                 2* (s23.* eigen(:,:,:,2,3,v,g)+ s13.* eigen(:,:,:,1,3,v,g)+ s12.* eigen(:,:,:,1,2,v,g))+...
                                                     s0yy.* eigen(:,:,:,2,2,v,g));
        deldphiAk(:,:,:,v,g)= fftn(deldphiA(:,:,:,v,g));
				
    end
end

% == judge whether plastic deformation begins
pl_stat= zeros(nx,ny,nz);                                                     % plastic status in different position, pl_stat= 1 if plastic deformation begin
ei_shear= zeros(nx,ny,nz,3,3); tot_e_shear= zeros(nx,ny,nz,3,3);
s_shear= zeros(nx,ny,nz,3,3);
krone_delta= [1 0 0; 0 1 0; 0 0 1];

s_Von2= s11.^2+ s22.^2+ s33.^2- s11.* s22- s11.* s33- s22.* s33+ 3* (s12.^2+ s13.^2+ s23.^2);
inrange= (s_Von2>= Sigmay2); pl_stat(inrange)= 1;

ei_shear(:,:,:,1,1)= ei1- ei_devia;
ei_shear(:,:,:,2,2)= ei2- ei_devia;
ei_shear(:,:,:,3,3)= ei3- ei_devia;
ei_shear(:,:,:,2,3)= ei4; ei_shear(:,:,:,3,2)= ei4;
ei_shear(:,:,:,1,3)= ei5; ei_shear(:,:,:,3,1)= ei5;
ei_shear(:,:,:,1,2)= ei6; ei_shear(:,:,:,2,1)= ei6;

tot_e_shear(:,:,:,1,1)= het_e(:,:,:,1)- e_devia; 
tot_e_shear(:,:,:,2,2)= het_e(:,:,:,2)- e_devia;
tot_e_shear(:,:,:,3,3)= het_e(:,:,:,3)- e_devia;
tot_e_shear(:,:,:,2,3)= het_e(:,:,:,4); tot_e_shear(:,:,:,3,2)= het_e(:,:,:,4);
tot_e_shear(:,:,:,1,3)= het_e(:,:,:,5); tot_e_shear(:,:,:,3,1)= het_e(:,:,:,5);
tot_e_shear(:,:,:,1,2)= het_e(:,:,:,6); tot_e_shear(:,:,:,2,1)= het_e(:,:,:,6);

e_shear= tot_e_shear- ei_shear;

for ii= 1: 3
    for jj= 1: 3
        for kk= 1: 3
            for ll= 1: 3

                s_shear(:,:,:,ii,jj)= s_shear(:,:,:,ii,jj)+ tot_Cpqrs(:,:,:,ii,jj,kk,ll).* e_shear(:,:,:,kk,ll);

            end
        end
    end
end

e1_pl= e1_pl- dt* pl_stat.* (Spq(:,:,:,1,1).* s_shear(:,:,:,1,1)* (krone_delta(1,1)/3-1)+ Spq(:,:,:,1,2).* s_shear(:,:,:,2,2)* (krone_delta(2,2)/3-1)+ Spq(:,:,:,1,3).* s_shear(:,:,:,3,3)* (krone_delta(3,3)/3-1)+...
                             Spq(:,:,:,1,4).* s_shear(:,:,:,2,3)* (krone_delta(2,3)/3-1)+ Spq(:,:,:,1,5).* s_shear(:,:,:,1,3)* (krone_delta(1,3)/3-1)+ Spq(:,:,:,1,6).* s_shear(:,:,:,1,2)* (krone_delta(1,2)/3-1)); 
e2_pl= e2_pl- dt* pl_stat.* (Spq(:,:,:,2,1).* s_shear(:,:,:,1,1)* (krone_delta(1,1)/3-1)+ Spq(:,:,:,2,2).* s_shear(:,:,:,2,2)* (krone_delta(2,2)/3-1)+ Spq(:,:,:,2,3).* s_shear(:,:,:,3,3)* (krone_delta(3,3)/3-1)+...
                             Spq(:,:,:,2,4).* s_shear(:,:,:,2,3)* (krone_delta(2,3)/3-1)+ Spq(:,:,:,2,5).* s_shear(:,:,:,1,3)* (krone_delta(1,3)/3-1)+ Spq(:,:,:,2,6).* s_shear(:,:,:,1,2)* (krone_delta(1,2)/3-1)); 
e3_pl= e3_pl- dt* pl_stat.* (Spq(:,:,:,3,1).* s_shear(:,:,:,1,1)* (krone_delta(1,1)/3-1)+ Spq(:,:,:,3,2).* s_shear(:,:,:,2,2)* (krone_delta(2,2)/3-1)+ Spq(:,:,:,3,3).* s_shear(:,:,:,3,3)* (krone_delta(3,3)/3-1)+...
                             Spq(:,:,:,3,4).* s_shear(:,:,:,2,3)* (krone_delta(2,3)/3-1)+ Spq(:,:,:,3,5).* s_shear(:,:,:,1,3)* (krone_delta(1,3)/3-1)+ Spq(:,:,:,3,6).* s_shear(:,:,:,1,2)* (krone_delta(1,2)/3-1)); 
e4_pl= e4_pl- dt* pl_stat.* (Spq(:,:,:,4,1).* s_shear(:,:,:,1,1)* (krone_delta(1,1)/3-1)+ Spq(:,:,:,4,2).* s_shear(:,:,:,2,2)* (krone_delta(2,2)/3-1)+ Spq(:,:,:,4,3).* s_shear(:,:,:,3,3)* (krone_delta(3,3)/3-1)+...
                             Spq(:,:,:,4,4).* s_shear(:,:,:,2,3)* (krone_delta(2,3)/3-1)+ Spq(:,:,:,4,5).* s_shear(:,:,:,1,3)* (krone_delta(1,3)/3-1)+ Spq(:,:,:,4,6).* s_shear(:,:,:,1,2)* (krone_delta(1,2)/3-1))/2; 
e5_pl= e5_pl- dt* pl_stat.* (Spq(:,:,:,5,1).* s_shear(:,:,:,1,1)* (krone_delta(1,1)/3-1)+ Spq(:,:,:,5,2).* s_shear(:,:,:,2,2)* (krone_delta(2,2)/3-1)+ Spq(:,:,:,5,3).* s_shear(:,:,:,3,3)* (krone_delta(3,3)/3-1)+...
                             Spq(:,:,:,5,4).* s_shear(:,:,:,2,3)* (krone_delta(2,3)/3-1)+ Spq(:,:,:,5,5).* s_shear(:,:,:,1,3)* (krone_delta(1,3)/3-1)+ Spq(:,:,:,5,6).* s_shear(:,:,:,1,2)* (krone_delta(1,2)/3-1))/2; 
e6_pl= e6_pl- dt* pl_stat.* (Spq(:,:,:,6,1).* s_shear(:,:,:,1,1)* (krone_delta(1,1)/3-1)+ Spq(:,:,:,6,2).* s_shear(:,:,:,2,2)* (krone_delta(2,2)/3-1)+ Spq(:,:,:,6,3).* s_shear(:,:,:,3,3)* (krone_delta(3,3)/3-1)+...
                             Spq(:,:,:,6,4).* s_shear(:,:,:,2,3)* (krone_delta(2,3)/3-1)+ Spq(:,:,:,6,5).* s_shear(:,:,:,1,3)* (krone_delta(1,3)/3-1)+ Spq(:,:,:,6,6).* s_shear(:,:,:,1,2)* (krone_delta(1,2)/3-1))/2;
 
if (conver>= tolerance)

    deldphiAk= inf(nx, ny, nz, variants, grainBs);
	
end

end  % end function

   

