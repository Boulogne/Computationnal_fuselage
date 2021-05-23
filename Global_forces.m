function [Fg_b,Fg_p] = Global_forces(ndof,n_nod,n_nod_b,Td_b,Beam_elements,beam_parameters,dat_beams,n_nod_p,Td_p,Plate_elements,plates_parameters,dat_plates,qE_b,qE_p,Re_b,Re_p)
%GLOBAL_FORCES Summary of this function goes here
Fg_b=zeros(ndof*n_nod,1);
Fg_p=zeros(ndof*n_nod,1);
% 
fe_p = zeros(n_nod_p*ndof, Plate_elements);
fe_b = zeros(n_nod_b*ndof, Beam_elements);

R_fun = @(alpha, beta, gamma) [cos(beta)*cos(gamma),   cos(beta)*sin(gamma),                            sin(beta)       ; 
                          -sin(alpha)*sin(beta)*cos(gamma) - cos(alpha)*sin(gamma), -sin(alpha)*sin(beta)*sin(gamma) + cos(alpha)*cos(gamma), sin(alpha)*cos(beta);
                          -cos(alpha)*sin(beta)*cos(gamma) + sin(alpha)*sin(gamma), -cos(alpha)*sin(beta)*sin(gamma) - sin(alpha)*cos(gamma), cos(alpha)*cos(beta)];

                      

for e = 1:Beam_elements
    L = beam_parameters.L(e);
    R  = R_fun(dat_beams(e,1),dat_beams(e,2),dat_beams(e,3));
    qE = R*transpose(qE_b(e,:));
    
    fe_b_pr([1,7]) = 0.5*qE(1)*L*[1; 1];
    fe_b_pr([2,6,8,12]) = 0.5*qE(2)*L*[1; L/6; 1; -L/6];
    fe_b_pr([3,5,9,11]) = 0.5*qE(3)*L*[1; -L/6; 1; L/6];
    
    fe_b(:,e) = transpose(Re_b(:,:,e))*fe_b_pr';
end

for e = 1:Plate_elements
    a=plates_parameters.a(e);
    b=plates_parameters.b(e);
    R  = R_fun(dat_plates(e,1),dat_plates(e,2),dat_plates(e,3));
    qE = R*transpose(qE_p(e,:));
    fe_p_pr=zeros(1,ndof*n_nod_p); % size of fe_p_pr seems to not be full with pdf
    fe_p_pr([1,2,7,8,13,14,19,20])=qE(1)*a*b*[1;0;1;0;1;0;1;0]+qE(3)*a*b*[0;1;0;1;0;1;0;1];
    fe_p_pr([3,4,5,9,10,11,15,16,17,21,22,23])=a*b*qE(3)*[1;b/3;a/3;1;-b/3;a/3;1;-b/3;-a/3;1;b/3;-a/3];
    
    fe_p(:,e) = transpose(Re_p(:,:,e))*fe_p_pr';
end

for e = 1:Beam_elements
    for i=1:n_nod_b*ndof
        I=Td_b(e,i); %corresponding global degree of freedom
        Fg_b(I)=Fg_b(I)+fe_b(i,e);
    end
end

for e = 1:Plate_elements
    for i=1:n_nod_p*ndof
        I=Td_p(e,i); %corresponding global degree of freedom
        Fg_p(I)=Fg_p(I)+fe_p(i,e);
    end
end

end

