function [q_nose_x,q_tail_x,q_tail_z,L_nose,L_tail] = load_nosetail(Tnose,Ttail,Beam_elements , beam_parameters)
%LOAD_NOSETAIL Summary of this function goes here
%   Detailed explanation goes here

L_nose=2*sum(beam_parameters.L(Tnose));
L_tail=2*sum(beam_parameters.L(Ttail));
rho=1.225;
V=210;
Cd=0.39;
S=12.84;

q_nose_x=zeros(Beam_elements,1);

D_nose=(rho*(V^2)*S*Cd)/2;
q_nose_x(Tnose)=(1/L_nose)*D_nose;% N/m % Tnose

q_tail_x=zeros(Beam_elements,1);
q_tail_z=zeros(Beam_elements,1);

Dl_tail=164.01*1e3; % N
D_tail=56.98*1e3; % N %56.98e3 new old 17.58e3

q_tail_x(Ttail)=(1/L_tail)*D_tail; % N/m % Ttail
q_tail_z(Ttail)=(1/L_tail)*Dl_tail; % N/m % Ttail
end

