%function [u, RR, F] = Global_System(Fext, KG, x)
function [u,R_R]=Global_System(vL,vR,K_G,F_G)
    K_LL=K_G(vL,vL);
    K_LR=K_G(vL,vR);
    K_RL=K_G(vR,vL);
    K_RR=K_G(vR,vR);
    F_L=F_G(vL,1);
    F_R=F_G(vR,1);

    u_R=zeros(size(vR,2),1); %imposed displacement vector
    u_L=K_LL\(F_L-K_LR*u_R); %free displacement vector
    
    R_R=K_RR*u_R+K_RL*u_L-F_R; %reactions vector

    %generalized displacement vector
    u(vL,1)=u_L;
    u(vR,1)=u_R;

end