function [N_p,Qy_p,Qz_p,T_p,My_p,Mz_p,uint_plates] = internalforces_plates2(xnodes,Tplates,Td_p,u,Kel_p,Re_p)
%INTERNALFORCES_BEAMS Summary of this function goes here
% uint_beams=zeros(ndof*n_nod_b,Beam_elements)
for e=1:size(Tplates,1)
   x1=xnodes(Tplates(e,1),1); x2=xnodes(Tplates(e,2),1);x3=xnodes(Tplates(e,3),1); x4=xnodes(Tplates(e,4),1);  %x-pos of the nodes of each bar
   y1=xnodes(Tplates(e,1),2); y2=xnodes(Tplates(e,2),2);y3=xnodes(Tplates(e,3),2);y4=xnodes(Tplates(e,4),2); %y-pos of the nodes of each bar
   z1=xnodes(Tplates(e,1),3); z2=xnodes(Tplates(e,2),3);z3=xnodes(Tplates(e,3),3);z4=xnodes(Tplates(e,4),3); %z-pos of the nodes of each bar   
 %% ELEMENT'S DISPLACEMENTS IN GLOBAL COORDINATES
 for r=1:size(Td_p,2)
        I = Td_p(e,r);
        u_e(r,1) = u(I,1);
 end
Fint_beams(:,e)=Re_p(:,:,e)*Kel_p(:,:,e)*u_e;    
uint_plates(:,e)=Re_p(:,:,e)*u_e;
  
end
    N_p = Fint_beams([1,7],:);
    Qy_p = Fint_beams([2,8],:);
    Qz_p = Fint_beams([3,9],:);
    T_p = Fint_beams([4,10],:);
    My_p = Fint_beams([5,11],:);
    Mz_p = Fint_beams([6,12],:); 
end