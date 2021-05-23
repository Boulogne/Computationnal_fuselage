function [q_front_x,q_front_y,q_front_z,q_rear_x,q_rear_y,q_rear_z,Meq, M_r, M_f,H,L,W] = loadstrby_wing(xnodes,Tfront,Trear,Beam_elements,beam_parameters)
%LOADSTRBY_WING Summary of this function goes here
%   Detailed explanation goes here
F_f=[78.94 ; 72.22 ; -71.07]*1e3;% N 
F_r=[17.18 ; -72.22 ; -25.47]*1e3;% N
M_f=[-349.08 ; 194.87 ; -86.29]*1e3;% N.m
M_r=[-203.02 ; 83.85; -91.93]*1e3;% N.m

L_front=sum(beam_parameters.L(Tfront));% m
L_rear=sum(beam_parameters.L(Trear));% Same but logical

q_front_x=zeros(Beam_elements,1);
q_front_y=zeros(Beam_elements,1);
q_front_z=zeros(Beam_elements,1);

q_rear_x=zeros(Beam_elements,1);
q_rear_y=zeros(Beam_elements,1);
q_rear_z=zeros(Beam_elements,1);

q_front_x(Tfront)=-F_f(1)/L_front; % N/m
q_front_y(Tfront)=-F_f(2)/L_front; % N/m
q_front_z(Tfront)=-F_f(3)/L_front; % N/m

q_rear_x(Trear)=-F_r(1)/L_rear; % N/m
q_rear_y(Trear)=-F_r(2)/L_rear; % N/m
q_rear_z(Trear)=-F_r(3)/L_rear; % N/m


H=xnodes(1323,3)-xnodes(1340,3);
L=xnodes(1043,1)-xnodes(1328,3);
W=xnodes(1323,2)-xnodes(1340,2);

Meq_z=zeros(Beam_elements,1);

Meq=-(M_f(3)+M_r(3))-(M_f(2)+M_r(2))*W/H; % N.m on z axes

end

