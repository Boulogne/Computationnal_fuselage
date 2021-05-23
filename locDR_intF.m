function [Fint_e,Uint_e,N, Qy,Qz,T,My,Mz]=locDR_intF(x,Tnod,Tdof,u,Kel,Re)
%% COMPUTE ROTATION MATRIX
for e=1:size(Tnod,1)
    x1=x(Tnod(e,1),1); x2=x(Tnod(e,2),1); %x-pos of the nodes of each bar
    y1=x(Tnod(e,1),2); y2=x(Tnod(e,2),2); %y-pos of the nodes of each bar
    z1=x(Tnod(e,1),3); z2=x(Tnod(e,2),3); %z-pos of the nodes of each bar
    
    %% ELEMENT'S DISPLACEMENTS IN GLOBAL COORDINATES
    for r=1:size(Tdof,2)
        I = Tdof(e,r);
        u_e(r,1) = u(I,1);
    end
    
    %% INTERNAL FORCES, DISPLACEMENTS AND ROTATIONS IN LOCAL COORDINATES
    Fint_e(:,e)=Re(:,:,e)*Kel(:,:,e)*u_e;
    N(e) = -Fint_e(1);
    Qy(e) = -Fint_e(2);
    Qz(e) = -Fint_e(3);
    T(e) = -Fint_e(4);
    My(:,e) = transpose([-Fint_e(5), -Fint_e(11)]);
    Mz(:,e) = transpose([-Fint_e(6), -Fint_e(12)]);
    Uint_e(:,e)=Re(:,:,e)*u_e;
    
    %Axial deformation as:
    %epsilon(e) = 1/elem_params.L(e)*(Uint_e(7,e)-Uint_e(1,e));
    
end

end