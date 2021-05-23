function [W_beams,W_plates] = structural_weight(Beam_elements,Plate_elements,Tmat_beams,Tmat_plates,mat_beams,mat_plates,beam_parameters,plates_parameters)
%STRUCTURAL_WEIGHT Summary of this function goes here
%   calculation of the structual weight for beams and plates 

% Plates thickness
hs = 4e-3; %Outer skin
hf = 8e-3;% Cabin floor
Ms=22900; %[kg]
g=9.81;

W_beams=zeros(Beam_elements,1);
rho_hat=(Ms-beam_parameters.M-plates_parameters.M)/(beam_parameters.V+plates_parameters.V);
for e=1:Beam_elements
    if Tmat_beams(e)==1
        W_beams(e)=(mat_beams(1,1)+rho_hat)*beam_parameters.Aa*g;
  
    elseif Tmat_beams(e)==2
        W_beams(e)=(mat_beams(2,1)+rho_hat)*beam_parameters.Ab*g;
    
    elseif Tmat_beams(e)==3
        W_beams(e)=(mat_beams(3,1)+rho_hat)*beam_parameters.Ac*g;
    end
end
% Force on Z axes for beams  % W_beams N/m


W_plates=zeros(Plate_elements,1);
for e=1:Plate_elements
    if Tmat_plates(e)==1
        W_plates(e)=(mat_plates(1,1)+rho_hat)*hs*g;
   
    elseif Tmat_plates(e)==2
        W_plates(e)=(mat_plates(2,1)+rho_hat)*hf*g;
    end
end

end

