clear
close all
clc
% Load mesh data
input_fuselage

% Material properties
% Not the same as in the pdf
rho1 =   2650; % kg/m3
E1   = 70.1e9; % Pa
nu1  =   0.30; % Poisson's ratio
rho2 =   2910; % kg/m3
E2   =   76.3e9; % Pa
nu2  =   0.35; % Poisson's ratio

% Section A parameters
h1 = 74e-3;
t1 = 2.0e-3;
a  = 20e-3;
% Get Aa, Iya, Iza, Ja ...

% Section B parameters
h2 = 24e-3;
t2 = 2.5e-3;
b  = 18e-3;
c  = 20e-3;
% Get Ab, Iyb, Izb, Jb ...

% Section C parameters
h3 = 48e-3;
t3 = 2.2e-3;
d  = 22e-3;

% Get Ac, Iyc, Izc, Jc ...
data_mat = [
%  Density  Young   Poisson     h     t     a     b
      rho1,    E1,      nu1,   h1,   t1,   a,   0; % Frames (1)
      rho1,    E1,      nu1,   h2,   t2,   b,   c; % Stringers (2)
      rho2,    E2,      nu2,   h3,   t3,   d,   0; % Reinforcements (3)
];
Beam_elements = size(Tbeams,1);
[beam_parameters] = beam_parameters(xnodes, Tbeams, Tmat_beams, data_mat);

beam_parameters.L;

% Beam properties
mat_beams = [
%  Density  Young   Poisson     A     Iy     Iz     J
      rho1,    E1,      nu1,   beam_parameters.Aa,   beam_parameters.Iya,   beam_parameters.Iza,   beam_parameters.Ja; % Frames (1)
      rho1,    E1,      nu1,   beam_parameters.Ab,   beam_parameters.Iyb,   beam_parameters.Izb,   beam_parameters.Jb; % Stringers (2)
      rho2,    E2,      nu2,   beam_parameters.Ac,   beam_parameters.Iyc,   beam_parameters.Izc,   beam_parameters.Jc; % Reinforcements (3)
];

% Plates thickness
hs = 4e-3; %Outer skin
hf = 8e-3;% Cabin floor

% Plate properties
mat_plates = [
%  Density  Young   Poisson     h
      rho2,    E2,      nu2,   hs; % Skin (1)
      rho2,    E2,      nu2,   hf; % Floor (2)
];
Plate_elements = size(Tplates,1);
 [plates_parameters] = plates_parameters(xnodes, Tplates, Tmat_plates, mat_plates);
 
Ms=22900; %[kg]
Mp=13500; %[kg]
g=9.81;


ndof = 6;                   % DOF per node (3D: 3 displacements and 3 rotations)
n_nod_b = 2;% # nodes per element.
n_nod_p = 4;
n_nod=size(xnodes,1);

Td_b=compute_Tdof(Beam_elements, ndof, n_nod_b, Tbeams);
Td_p=compute_Tdof(Plate_elements, ndof, n_nod_p, Tplates);
%% Stiffness matrix computation


%% Write code for beam elements ...
% Define the rotation matrix as afunction
R_fun = @(alpha, beta, gamma) [cos(beta)*cos(gamma),   cos(beta)*sin(gamma),                            sin(beta)       ; 
                          -sin(alpha)*sin(beta)*cos(gamma) - cos(alpha)*sin(gamma), -sin(alpha)*sin(beta)*sin(gamma) + cos(alpha)*cos(gamma), sin(alpha)*cos(beta);
                          -cos(alpha)*sin(beta)*cos(gamma) + sin(alpha)*sin(gamma), -cos(alpha)*sin(beta)*sin(gamma) - sin(alpha)*cos(gamma), cos(alpha)*cos(beta)];

[Kg_b,Re_b,Kel_b] = Stiffness_Kg_beam(n_nod_b,ndof,n_nod,Beam_elements,beam_parameters,Td_b,dat_beams);
%% Write code for plate elements ...
% Define the rotation matrix as afunction
[Kg_p,Re_p,Kel_p] = Stiffness_Kg_plate(n_nod_p,ndof,n_nod,Plate_elements,plates_parameters,Td_p,dat_plates);
Kg=Kg_b+Kg_p;

%% A) Structural weight
[W_beams,W_plates] = structural_weight(Beam_elements,Plate_elements,Tmat_beams,Tmat_plates,mat_beams,mat_plates,beam_parameters,plates_parameters);
%% B) Weight of the cabin passengers
 [W_p] = passengers_weight(Plate_elements,plates_parameters,Tfloor);
%% C) Loads transmitted by the wing
[q_front_x,q_front_y,q_front_z,q_rear_x,q_rear_y,q_rear_z,Meq, M_r, M_f,H,L_wing,W] = loadstrby_wing(xnodes,Tfront,Trear,Beam_elements,beam_parameters);
%% D) Loads transmitted by the nose and the tail cone
[q_nose_x,q_tail_x,q_tail_z,L_nose,L_tail] = load_nosetail(Tnose,Ttail,Beam_elements , beam_parameters);
%% E) Cabin pressure
[p_skin_x,p_skin_y,p_skin_z,p_nose_x,p_tail_x] = cabinpressure(Tskin,dat_plates,Plate_elements,Beam_elements,L_nose,L_tail,Tnose,Ttail);
%% Global Forces matrix
Fg_b=zeros(ndof*n_nod,1);
Fg_p=zeros(ndof*n_nod,1);

fe_p = zeros(n_nod_p*ndof, Plate_elements);
fe_b = zeros(n_nod_b*ndof, Beam_elements);

% qE_b Loads distribution [N/m]% can be an issue beam and plates
%Case 1
qE_b = zeros(Beam_elements,3);                    % q = [qxE; qyE; qzE]
qE_b(:,1)=q_front_x + q_rear_x;      % in x -> 
qE_b(:,2)= q_front_y+q_rear_x;                                                  % in y -> 
qE_b(:,3)= q_front_z+q_rear_x;                    % in z -> 

% qE_p Loads distribution [N/m2]
%Case 1 
qE_p = zeros(Plate_elements,3);                    % q = [qxE; qyE; qzE]
qE_p(:,1)=0 ;                                        % in x -> 
qE_p(:,2)=0;                                         % in y -> 
qE_p(:,3)= 0;        % in z -> 

[Fg_b,Fg_p] = Global_forces(ndof,n_nod,n_nod_b,Td_b,Beam_elements,beam_parameters,dat_beams,n_nod_p,Td_p,Plate_elements,plates_parameters,dat_plates,qE_b,qE_p,Re_b,Re_p);

Fg=Fg_b+Fg_p;

Fg((1058-1)*6+1)=Fg((1058-1)*6+1) - M_r(2)/H;
Fg((1058-1)*6+2)=Fg((1058-1)*6+2) - M_r(1)/H;
Fg((1043-1)*6+2)=Fg((1043-1)*6+2) + Meq/L_wing;
Fg((1035-1)*6+1)=Fg((1035-1)*6+1) - M_r(2)/H;
Fg((1035-1)*6+2)=Fg((1035-1)*6+2 )- M_r(1)/H;
Fg((1340-1)*6+1)=Fg((1340-1)*6+1) - M_f(2)/H;
Fg((1340-1)*6+2)=Fg((1340-1)*6+2 )- M_f(1)/H;
Fg((1328-1)*6+2)=Fg((1328-1)*6+2) + Meq/L_wing;
Fg((1323-1)*6+1)=Fg((1323-1)*6+1) - M_f(2)/H;
Fg((1323-1)*6+2)=Fg((1323-1)*6+2) - M_f(1)/H;

%% Prescribed degrees of freedom
 VR=[];
for i = 1:length(Tsym)
    Tsym(i);
%     if xnodes(Tsym(i),3) ==1.9750
%         VR{i} = [ 6*(Tsym(i)-1)+3 , 6*(Tsym(i)-1)+1,6*(Tsym(i)-1)+2];
        
    if Tsym(i) == 1921
         VR{i} = [ 6*(Tsym(i)-1)+1 ,6*(Tsym(i)-1)+3];   
         
    %     if Tsym(i) == 932 % New boundary condition 
%          VR{i} = [ 6*(Tsym(i)-1)+1,6*(Tsym(i)-1)+3]; 
    else
        VR{i} = [6*(Tsym(i)-1)+2 ,  6*(Tsym(i)-1)+4 ,6*(Tsym(i)-1)+6];  
    end
end
VR=cell2mat(VR);
VR=sort(VR);
VL= setdiff(1:ndof*n_nod , VR);

[u,R_R]=Global_System(VL,VR,Kg,Fg);
%% Local displacements, rotations and internal forces distribution    
[N_b,Qy_b,Qz_b,T_b,My_b,Mz_b,uint_beams] = internalforces_beams(xnodes,Tbeams,Td_b,u,Kel_b,Re_b);
[uint_plates] = internalforces_plates(xnodes,Tplates,Td_p,u,Re_p);

%% Postprocess
plotFuselage(xnodes,Tbeams,Tplates,Tmat_beams,Tmat_plates,u,uint_beams,N_b,Qy_b,Qz_b,T_b,My_b,Mz_b,uint_plates,mat_plates)