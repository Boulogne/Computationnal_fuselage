clear
close all
clc
% Load mesh data
input_fuselage

% Material properties
% Not the same as in the pdf
rho1 =   2700; % kg/m3
E1   = 68.9e9; % Pa
nu1  =   0.33; % Poisson's ratio
rho2 =   2810; % kg/m3
E2   =   70e9; % Pa
nu2  =   0.33; % Poisson's ratio

% Section A parameters
h1 = 75e-3;
t1 = 2.5e-3;
a  = 25e-3;
% Get Aa, Iya, Iza, Ja ...

% Section B parameters
h2 = 25e-3;
t2 = 2e-3;
b  = 16e-3;
c  = 21e-3;
% Get Ab, Iyb, Izb, Jb ...

% Section C parameters
h3 = 50e-3;
t3 = 3e-3;
d  = 25e-3;
% Get Ac, Iyc, Izc, Jc ...
data_mat = [
%  Density  Young   Poisson     h     t     a     b
      rho1,    E1,      nu1,   h1,   t1,   a,   0; % Frames (1)
      rho1,    E1,      nu1,   h2,   t2,   b,   c; % Stringers (2)
      rho2,    E2,      nu2,   h3,   t3,   d,   0; % Reinforcements (3)
];
Beam_elements = size(Tbeams,1);
[beam_parameters] = beam_parameters(xnodes, Tbeams, Tmat_beams, data_mat);

% beam_parameters.L;

% Beam properties
mat_beams = [
%  Density  Young   Poisson     A     Iy     Iz     J
      rho1,    E1,      nu1,   beam_parameters.Aa,   beam_parameters.Iya,   beam_parameters.Iza,   beam_parameters.Ja; % Frames (1)
      rho1,    E1,      nu1,   beam_parameters.Ab,   beam_parameters.Iyb,   beam_parameters.Izb,   beam_parameters.Jb; % Stringers (2)
      rho2,    E2,      nu2,   beam_parameters.Ac,   beam_parameters.Iyc,   beam_parameters.Izc,   beam_parameters.Jc; % Reinforcements (3)
];

% Plates thickness
hs = 4e-3; %Outer skin [m]
hf = 8e-3;% Cabin floor [m]

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

R_fun = @(alpha, beta, gamma) [cos(beta)*cos(gamma),   cos(beta)*sin(gamma),                            sin(beta)       ; 
                          -sin(alpha)*sin(beta)*cos(gamma) - cos(alpha)*sin(gamma), -sin(alpha)*sin(beta)*sin(gamma) + cos(alpha)*cos(gamma), sin(alpha)*cos(beta);
                          -cos(alpha)*sin(beta)*cos(gamma) + sin(alpha)*sin(gamma), -cos(alpha)*sin(beta)*sin(gamma) - sin(alpha)*cos(gamma), cos(alpha)*cos(beta)];

%% Stiffness matrix computation


%% Write code for beam elements ...
% Define the rotation matrix as afunction

[Kg_b,Re_b,Kel_b] = Stiffness_Kg_beam(n_nod_b,ndof,n_nod,Beam_elements,beam_parameters,Td_b,dat_beams);

%% Write code for plate elements ...
% Define the rotation matrix as afunction

[Kg_p,Re_p,Kel_p] = Stiffness_Kg_plate(n_nod_p,ndof,n_nod,Plate_elements,plates_parameters,Td_p,dat_plates);

Kg=Kg_b+Kg_p;




%% A) Structural weight

[W_beams,W_plates] = structural_weight(Beam_elements,Plate_elements,Tmat_beams,Tmat_plates,mat_beams,mat_plates,beam_parameters,plates_parameters);

% Solve problem for loading case A ...

%% B) Weight of the cabin passengers

 [W_p] = passengers_weight(Plate_elements,plates_parameters,Tfloor);
% Solve problem for loading case B ...

%% C) Loads transmitted by the wing

[q_front_x,q_front_y,q_front_z,q_rear_x,q_rear_y,q_rear_z,Meq, M_r, M_f,H,L_wing,W] = loadstrby_wing(xnodes,Tfront,Trear,Beam_elements,beam_parameters);

% Solve problem for loading case C ... 

%% D) Loads transmitted by the nose and the tail cone

% Solve problem for loading case D ...
[q_nose_x,q_tail_x,q_tail_z,L_nose,L_tail] = load_nosetail(Tnose,Ttail,Beam_elements , beam_parameters);

%% E) Cabin pressure

[p_skin_x,p_skin_y,p_skin_z,p_nose_x,p_tail_x] = cabinpressure(Tskin,dat_plates,Plate_elements,Beam_elements,L_nose,L_tail,Tnose,Ttail);

% Solve problem for loading case E ...

%% Global Forces matrix


% qE_b Loads distribution [N/m]% can be an issue beam and plates
%Case 1
qE_b = zeros(Beam_elements,3);                    % q = [qxE; qyE; qzE]
qE_b(:,1)= p_tail_x+p_nose_x;%(q_tail_x+p_nose_x+p_tail_x+q_front_x+q_rear_x) ;      % in x -> 
qE_b(:,2)= 0;%(q_front_y+q_rear_y);                                                  % in y -> 
qE_b(:,3)= 0;%(W_beams+q_tail_z+q_front_z+q_rear_z );                    % in z -> 

% qE_p Loads distribution [N/m2]
%Case 1 
qE_p = zeros(Plate_elements,3);                    % q = [qxE; qyE; qzE]
qE_p(:,1)=p_skin_x ;                                        % in x -> 
qE_p(:,2)=p_skin_y;                                         % in y -> 
qE_p(:,3)= p_skin_z;%(-W_plates -W_p + p_skin_z);        % in z -> 


[Fg_b,Fg_p] = Global_forces(ndof,n_nod,n_nod_b,Td_b,Beam_elements,beam_parameters,dat_beams,n_nod_p,Td_p,Plate_elements,plates_parameters,dat_plates,qE_b,qE_p,Re_b,Re_p);

Fg=Fg_b+Fg_p;

% Fg((1058-1)*6+1)=Fg((1058-1)*6+1) - M_r(2)/H;
% Fg((1058-1)*6+2)=Fg((1058-1)*6+2) - M_r(1)/H;
% Fg((1043-1)*6+2)=Fg((1043-1)*6+2) + Meq/L_wing;
% Fg((1035-1)*6+1)=Fg((1035-1)*6+1) - M_r(2)/H;
% Fg((1035-1)*6+2)=Fg((1035-1)*6+2 )- M_r(1)/H;
% Fg((1340-1)*6+1)=Fg((1340-1)*6+1) - M_f(2)/H;
% Fg((1340-1)*6+2)=Fg((1340-1)*6+2 )- M_f(1)/H;
% Fg((1328-1)*6+2)=Fg((1328-1)*6+2) + Meq/L_wing;
% Fg((1323-1)*6+1)=Fg((1323-1)*6+1) - M_f(2)/H;
% Fg((1323-1)*6+2)=Fg((1323-1)*6+2) - M_f(1)/H;

%% Prescribed degrees of freedom

 VR=[];
for i = 1:length(Tsym)
    A=Tsym(i);
    if Tsym(i) ==79
        VR{i} = [ 6*(Tsym(i)-1)+3 , 6*(Tsym(i)-1)+1];
     end
     if Tsym(i) == 1
        VR{i} = [ 6*(Tsym(i)-1)+1 ,6*(Tsym(i)-1)+3];    
    else
        VR{i} = [6*(A-1)+2 ,  6*(A-1)+4 ,6*(A-1)+6];  
    end
end
VR=cell2mat(VR);
VR=sort(VR);
VL= setdiff(1:ndof*n_nod , VR);

 [u,R_R]=Global_System(VL,VR,Kg,Fg);
   
% The following part as issues i think

%% Local displacements, rotations and internal forces distribution

% Tnod is not in fuselage we have Tbeams and Tplates
% 
uint_beams= zeros(ndof*n_nod_b,Beam_elements) ;
Fint_beams= zeros(ndof*n_nod_b,Beam_elements) ;
uint_plates= zeros(ndof*n_nod_p,Plate_elements) ;



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
    N_b = Fint_beams(7,:);
    Qy_b = Fint_beams(8,:);
    Qz_b = Fint_beams(9,:);
    T_b = Fint_beams(10,:);
    My_b = transpose([-Fint_beams(5,:), Fint_beams(11,:)]);
    Mz_b =transpose( [-Fint_beams(6,:), Fint_beams(12,:)]);


for e=1:size(Tplates,1)
   x1=xnodes(Tplates(e,1),1); x2=xnodes(Tplates(e,2),1);x3=xnodes(Tplates(e,3),1); x4=xnodes(Tplates(e,4),1);  %x-pos of the nodes of each bar
   y1=xnodes(Tplates(e,1),2); y2=xnodes(Tplates(e,2),2);y3=xnodes(Tplates(e,3),2);y4=xnodes(Tplates(e,4),2); %y-pos of the nodes of each bar
   z1=xnodes(Tplates(e,1),3); z2=xnodes(Tplates(e,2),3);z3=xnodes(Tplates(e,3),3);z4=xnodes(Tplates(e,4),3); %z-pos of the nodes of each bar
 
   %% ELEMENT'S DISPLACEMENTS IN GLOBAL COORDINATES
    for r=1:size(Td_p,2)
        I = Td_p(e,r);
        u_e(r,1) = u(I,1);
    end
    %% INTERNAL FORCES, DISPLACEMENTS AND ROTATIONS IN LOCAL COORDINATES
    uint_plates(:,e)=Re_p(:,:,e)*u_e;

end

%% Postprocess

plotFuselage(xnodes,Tbeams,Tplates,Tmat_beams,Tmat_plates,u,uint_beams,N_b,Qy_b,Qz_b,T_b,My_b,Mz_b,uint_plates,mat_plates)


