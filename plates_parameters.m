function  [plates_parameters] = plates_parameters(xnodes, Tplates, Tmat_plates, mat_plates)

Plate_elements = size(Tplates,1);
a_plate = zeros(Plate_elements,1); 
b_plate= zeros(Plate_elements,1); 
h_plate=zeros(Plate_elements,1); 
nu_plate=zeros(Plate_elements,1); 
E_plate=zeros(Plate_elements,1); 
rho_p=zeros(Plate_elements,1);
Mplate=zeros(Plate_elements,1); 
Vplate=zeros(Plate_elements,1); 
S_f=0;
for e=1:Plate_elements
    a_plate(e)=sqrt((xnodes(Tplates(e,2),1) - xnodes(Tplates(e,1),1))^2 + ...  % length 
                     (xnodes(Tplates(e,2),2) - xnodes(Tplates(e,1),2))^2 + ...
                     (xnodes(Tplates(e,2),3) - xnodes(Tplates(e,1),3))^2)/2;
    b_plate(e)=sqrt((xnodes(Tplates(e,4),1) - xnodes(Tplates(e,1),1))^2 + ...  % length 
                     (xnodes(Tplates(e,4),2) - xnodes(Tplates(e,1),2))^2 + ...
                     (xnodes(Tplates(e,4),3) - xnodes(Tplates(e,1),3))^2)/2;
    if Tmat_plates(e)==1% Skin
        h_plate(e)=mat_plates(1,4);  
        nu_plate(e)=mat_plates(1,3);
        E_plate(e)=mat_plates(1,2);
    end
    if Tmat_plates(e)==2% floor
         h_plate(e)=mat_plates(2,4);
         nu_plate(e)=mat_plates(2,3);
         E_plate(e)=mat_plates(2,2);
         S_f=S_f+4*a_plate(e)*b_plate(e);
    end
   Vplate(e)=4*a_plate(e)*b_plate(e)*h_plate(e);
   Mplate(e)=mat_plates(1,1)*Vplate(e);             
end
V_p=2*sum(Vplate);
M_p=2*sum(Mplate);
S_floor=2*S_f;

plates_parameters=struct('a',a_plate,'b',b_plate,'h',h_plate, 'nu',nu_plate, 'E',E_plate, 'V',V_p,'M',M_p,'Sfloor',S_floor); 
                     
                     
end


