function [N_b,Qy_b,Qz_b,T_b,My_b,Mz_b,uint_beams] = internalforces_beams(xnodes,Tbeams,Td_b,u,Kel_b,Re_b)
%INTERNALFORCES_BEAMS Summary of this function goes here
% uint_beams=zeros(ndof*n_nod_b,Beam_elements)
for e=1:size(Tbeams,1)
   x1=xnodes(Tbeams(e,1),1); x2=xnodes(Tbeams(e,2),1); %x-pos of the nodes of each bar
   y1=xnodes(Tbeams(e,1),2); y2=xnodes(Tbeams(e,2),2); %y-pos of the nodes of each bar
   z1=xnodes(Tbeams(e,1),3); z2=xnodes(Tbeams(e,2),3); %z-pos of the nodes of each bar
 %% ELEMENT'S DISPLACEMENTS IN GLOBAL COORDINATES
    for r=1:size(Td_b,2)
        I = Td_b(e,r);
        u_e(r,1) = u(I,1);
    end
    %% INTERNAL FORCES, DISPLACEMENTS AND ROTATIONS IN LOCAL COORDINATES
    Fint_beams(:,e)=Re_b(:,:,e)*Kel_b(:,:,e)*u_e;
    uint_beams(:,e)=Re_b(:,:,e)*u_e;
end
    N_b = Fint_beams([1,7],:);
    Qy_b = Fint_beams([2,8],:);
    Qz_b = Fint_beams([3,9],:);
    T_b = Fint_beams([4,10],:);
    My_b = Fint_beams([5,11],:);
    Mz_b = Fint_beams([6,12],:); 
end

