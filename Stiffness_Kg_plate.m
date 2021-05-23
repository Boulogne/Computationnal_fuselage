function [Kg_p,Re_p,Kel_p] = Stiffness_Kg_plate(n_nod_p,ndof,n_nod,Plate_elements,plates_parameters,Td_p,dat_plates);
%STIFFNESS_KG_BEAM Summary of this function goes here

R_fun = @(alpha, beta, gamma) [cos(beta)*cos(gamma),   cos(beta)*sin(gamma),                            sin(beta)       ; 
                          -sin(alpha)*sin(beta)*cos(gamma) - cos(alpha)*sin(gamma), -sin(alpha)*sin(beta)*sin(gamma) + cos(alpha)*cos(gamma), sin(alpha)*cos(beta);
                          -cos(alpha)*sin(beta)*cos(gamma) + sin(alpha)*sin(gamma), -cos(alpha)*sin(beta)*sin(gamma) - sin(alpha)*cos(gamma), cos(alpha)*cos(beta)];


Ke_p = zeros(n_nod_p*ndof); % initialize global stiffness matrix
z = zeros(3,3);
%Re = zeros(10,12)
for e=1:Plate_elements % for e!ach element (1 to 79)
    % compute the local stiffness matrix for each beam element data
   Kl =  local_Ke_plates(plates_parameters.E(e),plates_parameters.nu(e),plates_parameters.a(e),plates_parameters.b(e),plates_parameters.h(e));%(beam_params.E(e),beam_params.A(e), beam_params.L(e),...
%                   beam_params.Iy(e),beam_params.Iz(e),...
%                   beam_params.nu(e),beam_params.J(e)); % local stiffness in local coordinates
    R_p  = R_fun(dat_plates(e,1),dat_plates(e,2),dat_plates(e,3));
    Re_p(1:3,1:3,e)=R_p; Re_p(4:6,4:6,e)=R_p;
    Re_p(7:9,7:9,e)=R_p; Re_p(10:12,10:12,e)=R_p;
    Re_p(13:15,13:15,e)=R_p; Re_p(16:18,16:18,e)=R_p;
    Re_p(19:21,19:21,e)=R_p; Re_p(22:24,22:24,e)=R_p;
    
    Ke_p = transpose(Re_p(:,:,e))*Kl*Re_p(:,:,e); % local stiff. matrix in local coordinates.
    for r=1:n_nod_p*ndof
        for s=1:n_nod_p*ndof
            Kel_p(r,s,e) = Ke_p(r,s); %bar (element) stiffness matrices
        end
    end
end

% Kg Assembly process:
Kg_p=zeros(ndof*n_nod,ndof*n_nod);  % Global Stiffness matrix Initialization
for e = 1:Plate_elements
   for i = 1:n_nod_p*ndof
       I = Td_p(e, i);            %corresponding global degree of freedom (rows)
       for j = 1:n_nod_p*ndof
           J = Td_p(e, j);        %corresponding global degree of freedom (columns)
           Kg_p(I,J) = Kg_p(I,J)+Kel_p(i, j, e);
       end
   end
end


end

