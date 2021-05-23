function [uint_plates] = internalforces_plates(xnodes,Tplates,Td_p,u,Re_p)
%INTERNALFORCES_PLATES Summary of this function goes here
for e=1:size(Tplates,1)
   x1=xnodes(Tplates(e,1),1); x2=xnodes(Tplates(e,2),1);x3=xnodes(Tplates(e,3),1); x4=xnodes(Tplates(e,4),1);  %x-pos of the nodes of each bar
   y1=xnodes(Tplates(e,1),2); y2=xnodes(Tplates(e,2),2);y3=xnodes(Tplates(e,3),2);y4=xnodes(Tplates(e,4),2); %y-pos of the nodes of each bar
   z1=xnodes(Tplates(e,1),3); z2=xnodes(Tplates(e,2),3);z3=xnodes(Tplates(e,3),3);z4=xnodes(Tplates(e,4),3); %z-pos of the nodes of each bar   
   for r=1:size(Td_p,2)
        I = Td_p(e,r);
        u_e(r,1) = u(I,1);
    end
uint_plates(:,e)=Re_p(:,:,e)*u_e;
    
end

