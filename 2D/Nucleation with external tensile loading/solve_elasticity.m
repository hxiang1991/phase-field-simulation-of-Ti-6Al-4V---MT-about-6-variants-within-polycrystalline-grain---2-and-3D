function [deldphiAk, E_elastic, e1_pl,e2_pl,e3_pl, output_s11, output_s22, output_s12] = solve_elasticity(phiA,StrucphiB,eigen,e1_pl,e2_pl,e3_pl,grainBs,variants,...
                                                                                                          tot_Cpq,hom_Cpq,het_Cpq,gomega,tot_Cpqrs,Spq,Sigmay2,nx,ny,kx,ky,dt,delta_Gm,Vmol)

format long;

niter= 1e2; tolerance= 1e-4;

%% == eigenstrains
ei1= zeros(nx, ny); ei2= zeros(nx, ny); ei3= zeros(nx, ny); 

% e1_pl= 0; e2_pl= 0; e3_pl= 0;

for g= 1: grainBs
    for v= 1: variants
        
        ei1= ei1+ StrucphiB(:,:,v,g).* eigen(:,:,1,1,v,g).* phiA(:,:,v,g); 
        ei2= ei2+ StrucphiB(:,:,v,g).* eigen(:,:,2,2,v,g).* phiA(:,:,v,g); 
        ei3= ei3+ StrucphiB(:,:,v,g).* eigen(:,:,1,2,v,g).* phiA(:,:,v,g); 
        
    end
end

ei1= ei1+ e1_pl; ei2= ei2+ e2_pl; ei3= ei3+ e3_pl;
ei_devia= (ei1+ ei2)/2;

%% == elastic strain and stress field
u= zeros(nx,ny,2,niter+1); uk= zeros(nx,ny,2); 
T0= zeros(nx,ny,3);
T= zeros(nx,ny,3); Tk= zeros(nx,ny,3);
hom_e= zeros(3,1);
het_e= zeros(nx,ny,3,niter+1); het_ek= zeros(nx,ny,3);
s= zeros(nx,ny,3); sk= zeros(nx,ny,3);

% == zeroth-order iteration
s(:,:,1)= hom_Cpq(1,1)* ei1+ hom_Cpq(1,2)* ei2+ 2* hom_Cpq(1,3)* ei3; sk(:,:,1)= fftn(s(:,:,1));
s(:,:,2)= hom_Cpq(2,1)* ei1+ hom_Cpq(2,2)* ei2+ 2* hom_Cpq(2,3)* ei3; sk(:,:,2)= fftn(s(:,:,2));
s(:,:,3)= hom_Cpq(3,1)* ei1+ hom_Cpq(3,2)* ei2+ 2* hom_Cpq(3,3)* ei3; sk(:,:,3)= fftn(s(:,:,3));

uk(:,:,1,1)= -1i*(gomega(:,:,1).*(sk(:,:,1).*kx+ sk(:,:,3).*ky)+ gomega(:,:,3).*(sk(:,:,3).*kx+ sk(:,:,2).*ky)); u(:,:,1,1)= real(ifftn(uk(:,:,1,1)));
uk(:,:,2,1)= -1i*(gomega(:,:,3).*(sk(:,:,1).*kx+ sk(:,:,3).*ky)+ gomega(:,:,2).*(sk(:,:,3).*kx+ sk(:,:,2).*ky)); u(:,:,2,1)= real(ifftn(uk(:,:,2,1)));

het_ek(:,:,1,1)= 1i* kx.* uk(:,:,1,1);                       
het_ek(:,:,2,1)= 1i* ky.* uk(:,:,2,1);                       
het_ek(:,:,3,1)= 0.5*1i* (kx.* uk(:,:,2,1)+ ky.* uk(:,:,1,1)); 

het_e(:,:,1,1)= real(ifftn(het_ek(:,:,1,1))); het_e(:,:,1,1)= het_e(:,:,1,1)- trapz(trapz(het_e(:,:,1,1)))/(nx*ny);
het_e(:,:,2,1)= real(ifftn(het_ek(:,:,2,1))); het_e(:,:,2,1)= het_e(:,:,2,1)- trapz(trapz(het_e(:,:,2,1)))/(nx*ny);
het_e(:,:,3,1)= real(ifftn(het_ek(:,:,3,1))); het_e(:,:,3,1)= het_e(:,:,3,1)- trapz(trapz(het_e(:,:,3,1)))/(nx*ny);

T0(:,:,1)= tot_Cpq(:,:,1,1).*(ei1- hom_e(1))+ tot_Cpq(:,:,1,2).*(ei2- hom_e(2))+ 2* tot_Cpq(:,:,1,3).*(ei3- hom_e(3));
T0(:,:,2)= tot_Cpq(:,:,2,1).*(ei1- hom_e(1))+ tot_Cpq(:,:,2,2).*(ei2- hom_e(2))+ 2* tot_Cpq(:,:,2,3).*(ei3- hom_e(3));
T0(:,:,3)= tot_Cpq(:,:,3,1).*(ei1- hom_e(1))+ tot_Cpq(:,:,3,2).*(ei2- hom_e(2))+ 2* tot_Cpq(:,:,3,3).*(ei3- hom_e(3));

for iter= 1: niter
    
    T(:,:,1)= T0(:,:,1)- het_Cpq(:,:,1,1).* het_e(:,:,1,iter)- het_Cpq(:,:,1,2).* het_e(:,:,2,iter)- 2* het_Cpq(:,:,1,3).* het_e(:,:,3,iter);
              
    T(:,:,2)= T0(:,:,2)- het_Cpq(:,:,2,1).* het_e(:,:,1,iter)- het_Cpq(:,:,2,2).* het_e(:,:,2,iter)- 2* het_Cpq(:,:,2,3).* het_e(:,:,3,iter);
              
    T(:,:,3)= T0(:,:,3)- het_Cpq(:,:,3,1).* het_e(:,:,1,iter)- het_Cpq(:,:,3,2).* het_e(:,:,2,iter)- 2* het_Cpq(:,:,3,3).* het_e(:,:,3,iter);
              
    Tk(:,:,1)= fftn(T(:,:,1)); Tk(:,:,2)= fftn(T(:,:,2)); Tk(:,:,3)= fftn(T(:,:,3));
    
    uk(:,:,1)= -1i*(gomega(:,:,1).*(Tk(:,:,1).*kx+ Tk(:,:,3).*ky)+ gomega(:,:,3).*(Tk(:,:,3).*kx+ Tk(:,:,2).*ky));
    uk(:,:,2)= -1i*(gomega(:,:,3).*(Tk(:,:,1).*kx+ Tk(:,:,3).*ky)+ gomega(:,:,2).*(Tk(:,:,3).*kx+ Tk(:,:,2).*ky));
    
    u(:,:,1,iter+1)= real(ifftn(uk(:,:,1))); u(:,:,2,iter+1)= real(ifftn(uk(:,:,2))); 
    
    het_ek(:,:,1)= 1i* kx.* uk(:,:,1);                              
    het_ek(:,:,2)= 1i* ky.* uk(:,:,2);                              
    het_ek(:,:,3)= 0.5*1i* (kx.* uk(:,:,2)+ ky.* uk(:,:,1)); 
    
    het_e(:,:,1,iter+1)= real(ifftn(het_ek(:,:,1))); het_e(:,:,1,iter+1)= het_e(:,:,1,iter+1)- trapz(trapz(het_e(:,:,1,iter+1)))/(nx*ny);
    het_e(:,:,2,iter+1)= real(ifftn(het_ek(:,:,2))); het_e(:,:,2,iter+1)= het_e(:,:,2,iter+1)- trapz(trapz(het_e(:,:,2,iter+1)))/(nx*ny);
    het_e(:,:,3,iter+1)= real(ifftn(het_ek(:,:,3))); het_e(:,:,3,iter+1)= het_e(:,:,3,iter+1)- trapz(trapz(het_e(:,:,3,iter+1)))/(nx*ny);
    
    % == check convergence
    diffu2= (u(:,:,1,iter+1)- u(:,:,1,iter)).^2+ (u(:,:,2,iter+1)- u(:,:,2,iter)).^2;
    conver= (trapz(trapz(diffu2)))^0.5;
    
    if(conver< tolerance)
       break
    end
        
end  %iter

e11= hom_e(1)+ het_e(:,:,1,iter+1)- ei1; 
e22= hom_e(2)+ het_e(:,:,2,iter+1)- ei2;
e12= hom_e(3)+ het_e(:,:,3,iter+1)- ei3;

e_devia= (het_e(:,:,1,iter+1)+ het_e(:,:,2,iter+1))/2;

% ==  stress field and strain energy
s11= tot_Cpq(:,:,1,1).* e11+ tot_Cpq(:,:,1,2).* e22+ 2* tot_Cpq(:,:,1,3).* e12; output_s11= s11* delta_Gm/Vmol/1e06;   % unit :MPa
s22= tot_Cpq(:,:,2,1).* e11+ tot_Cpq(:,:,2,2).* e22+ 2* tot_Cpq(:,:,2,3).* e12; output_s22= s22* delta_Gm/Vmol/1e06;   % unit :MPa
s12= tot_Cpq(:,:,3,1).* e11+ tot_Cpq(:,:,3,2).* e22+ 2* tot_Cpq(:,:,3,3).* e12; output_s12= s12* delta_Gm/Vmol/1e06;   % unit :MPa

deldphiA = zeros(nx, ny, variants, grainBs);
deldphiAk= zeros(nx, ny, variants, grainBs);

E_elastic= (s11.* e11+ s22.* e22+ 2* s12.* e12)* delta_Gm/Vmol; 

for g= 1: grainBs
    for v= 1: variants
        
        deldphiA(:,:,v,g)= -StrucphiB(:,:,v,g).*(s11.* eigen(:,:,1,1,v,g)+ s22.* eigen(:,:,2,2,v,g)+ 2* s12.* eigen(:,:,1,2,v,g));
        deldphiAk(:,:,v,g)= fftn(deldphiA(:,:,v,g));
				
    end
end

% == judge whether plastic deformation begins
pl_stat= zeros(nx,ny);                                                     % plastic status in different position, pl_stat= 1 if plastic deformation begin
ei_shear= zeros(nx,ny,2,2); tot_e_shear= zeros(nx,ny,2,2);
s_shear= zeros(nx,ny,2,2);
krone_delta= [1 0; 0 1];

s_Von2= s11.^2- s11.* s22+ s22.^2+ 3* s12.^2;
inrange= (s_Von2>= Sigmay2); pl_stat(inrange)= 1;

ei_shear(:,:,1,1)= ei1- ei_devia;
ei_shear(:,:,2,2)= ei2- ei_devia;
ei_shear(:,:,1,2)= ei3; 
ei_shear(:,:,2,1)= ei3;

tot_e_shear(:,:,1,1)= het_e(:,:,1,iter+1)- e_devia; 
tot_e_shear(:,:,2,2)= het_e(:,:,2,iter+1)- e_devia;
tot_e_shear(:,:,1,2)= het_e(:,:,3,iter+1);
tot_e_shear(:,:,2,1)= het_e(:,:,3,iter+1);

e_shear= tot_e_shear- ei_shear;

for ii= 1: 2
    for jj= 1: 2
        for kk= 1: 2
            for ll= 1: 2

                s_shear(:,:,ii,jj)= s_shear(:,:,ii,jj)+ tot_Cpqrs(:,:,ii,jj,kk,ll).* e_shear(:,:,kk,ll);

            end
        end
    end
end

e1_pl= e1_pl- dt* pl_stat.* (Spq(:,:,1,1).* s_shear(:,:,1,1)*(krone_delta(1,1)/2-1)+ Spq(:,:,1,2).* s_shear(:,:,2,2)*(krone_delta(2,2)/2-1)+ Spq(:,:,1,3).* s_shear(:,:,1,2)*(krone_delta(1,2)/2-1)); 
e2_pl= e2_pl- dt* pl_stat.* (Spq(:,:,2,1).* s_shear(:,:,1,1)*(krone_delta(1,1)/2-1)+ Spq(:,:,2,2).* s_shear(:,:,2,2)*(krone_delta(2,2)/2-1)+ Spq(:,:,2,3).* s_shear(:,:,1,2)*(krone_delta(1,2)/2-1));  
e3_pl= e3_pl- dt* pl_stat.* (Spq(:,:,3,1).* s_shear(:,:,1,1)*(krone_delta(1,1)/2-1)+ Spq(:,:,3,2).* s_shear(:,:,2,2)*(krone_delta(2,2)/2-1)+ Spq(:,:,3,3).* s_shear(:,:,1,2)*(krone_delta(1,2)/2-1))/2; 
 
if (conver>= tolerance)

    deldphiAk= inf(nx, ny, variants, grainBs);
	
end

end  % end function

   

