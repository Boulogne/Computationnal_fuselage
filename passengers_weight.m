function [W_p] = passengers_weight(Plate_elements,plates_parameters,Tfloor)
%PASSENGERS_WEIGHT Summary of this function goes here

Mp=13500; %[kg]
g=9.81;

W_p=zeros(Plate_elements,1);

W_p(Tfloor)=Mp*g/plates_parameters.Sfloor; % N/m2

end

