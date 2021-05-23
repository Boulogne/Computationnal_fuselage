function [p_skin_x,p_skin_y,p_skin_z,p_nose_x,p_tail_x] = cabinpressure(Tskin,dat_plates,Plate_elements,Beam_elements,L_nose,L_tail,Tnose,Ttail)
%CABINPRESSURE Summary of this function goes here
%   Detailed explanation goes here
pin=78191.21;
pout=22632.06;

Alpha=dat_plates(Tskin,1);
Beta=dat_plates(Tskin,2);
Gamma=dat_plates(Tskin,3);
ne =zeros(3,size(Tskin,1));
p_skin=zeros(3,size(Tskin,1));

for i=1:size(Tskin,1)
ne(:,i)=[-cos(Alpha(i))*sin(Beta(i))*cos(Gamma(i))+sin(Alpha(i))*sin(Gamma(i)) ; ....
    -cos(Alpha(i))*sin(Beta(i))*sin(Gamma(i))-sin(Alpha(i))*cos(Gamma(i)) ; ...
    cos(Alpha(i))*cos(Beta(i))];
    p_skin(:,i)=(pin-pout)*ne(:,i);
end
p_skin_x=zeros(Plate_elements,1);
p_skin_y=zeros(Plate_elements,1);
p_skin_z=zeros(Plate_elements,1);

p_skin_x(Tskin)=p_skin(1,:).'; % N/m2 plates on x
p_skin_y(Tskin)=p_skin(2,:).';% N/m2 plates on y
p_skin_z(Tskin)=p_skin(3,:).';% N/m2 plates on z

S=12.84;

p_nose_x=zeros(Beam_elements,1);
p_tail_x=zeros(Beam_elements,1);
p_nose_x(Tnose)=(1/L_nose)*(-(pin-pout)*S); % nose on x axes
p_tail_x(Ttail)=(1/L_tail)*((pin-pout)*S); % Tail on x axes

end

