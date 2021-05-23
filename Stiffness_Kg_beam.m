function [Kg_b,Re_b,Kel_b] = Stiffness_Kg_beam(n_nod_b,ndof,n_nod,Beam_elements,beam_parameters,Td_b,dat_beams)
%STIFFNESS_KG_BEAM Summary of this function goes here

R_fun = @(alpha, beta, gamma) [cos(beta)*cos(gamma),   cos(beta)*sin(gamma),                            sin(beta)       ; 
                          -sin(alpha)*sin(beta)*cos(gamma) - cos(alpha)*sin(gamma), -sin(alpha)*sin(beta)*sin(gamma) + cos(alpha)*cos(gamma), sin(alpha)*cos(beta);
                          -cos(alpha)*sin(beta)*cos(gamma) + sin(alpha)*sin(gamma), -cos(alpha)*sin(beta)*sin(gamma) - sin(alpha)*cos(gamma), cos(alpha)*cos(beta)];
Ke_b = zeros(n_nod_b*ndof); % initialize global stiffness matrix
z = zeros(3,3);
%Re = zeros(10,12)
for e=1:Beam_elements % for each Beam elements 
    % compute the local stiffness matrix for each beam element data
   Kl =  local_Ke_beams(beam_parameters.E(e),beam_parameters.A(e),beam_parameters.L(e),beam_parameters.Iy(e),beam_parameters.Iz(e),...
       beam_parameters.J(e),beam_parameters.nu(e));%(beam_params.E(e),beam_params.A(e), beam_params.L(e),...
%                   beam_params.Iy(e),beam_params.Iz(e),...
%                   beam_params.nu(e),beam_params.J(e)); % local stiffness in local coordinates
    R_b  = R_fun(dat_beams(e,1),dat_beams(e,2),dat_beams(e,3));
    Re_b(1:3,1:3,e)=R_b; Re_b(4:6,4:6,e)=R_b;
    Re_b(7:9,7:9,e)=R_b; Re_b(10:12,10:12,e)=R_b;
  
    Ke_b = transpose(Re_b(:,:,e))*Kl*Re_b(:,:,e); % local stiff. matrix in local coordinates.
    for r=1:n_nod_b*ndof
        for s=1:n_nod_b*ndof
            Kel_b(r,s,e) = Ke_b(r,s); %bar (element) stiffness matrices
        end
    end
end

% Kg Assembly process:
Kg_b=zeros(ndof*n_nod,ndof*n_nod);  % Global Stiffness matrix Initialization
for e = 1:Beam_elements
   for i = 1:n_nod_b*ndof
       I = Td_b(e, i);            %corresponding global degree of freedom (rows)
       for j = 1:n_nod_b*ndof
           J = Td_b(e, j);        %corresponding global degree of freedom (columns)
           Kg_b(I,J) = Kg_b(I,J) + Kel_b(i, j, e);
       end
   end
end


end

