function  [beam_parameters] = beam_parameters(xnodes, Tbeams, Tmat_beams, data_mat)

h1=data_mat( 1,4);
t1=data_mat( 1,5);
a=data_mat(1,6);

h2=data_mat( 2,4);
t2=data_mat( 2,5);
b=data_mat(2,6);
c=data_mat(2,7);

h3=data_mat(3,4);
t3=data_mat( 3,5);
d=data_mat(3,6);

Beam_elements = size(Tbeams,1);

    % Initialize vectors
    L_beam = zeros(Beam_elements,1);     % (a) Associated length
    A_beam = zeros(Beam_elements,1);     % (b) Associated section area
    % 
    Iy_beam = zeros(Beam_elements,1);    % (c-1) Associated inertia (y)
    Iz_beam = zeros(Beam_elements,1);    % (c-2) Associated inertia (z) 
    J_beam= zeros(Beam_elements,1);     % (d) Associated torsional constant (y)
    %
    V_beam= zeros(Beam_elements,1);     % (e) Beam element volume
    m_beam = zeros(Beam_elements,1);     % (f) Beam element mass
    rho_beam = zeros(Beam_elements,1);   % (g) Beam elements effective density
    E_beam = zeros(Beam_elements,1);     % (f) Beam element mass
    nu_beam = zeros(Beam_elements,1);   % (g) Beam elements effective density
    Mb = zeros(Beam_elements,1); 
for e=1:Beam_elements
       L_beam(e) = sqrt((xnodes(Tbeams(e,2),1) - xnodes(Tbeams(e,1),1))^2 + ...  % length of the beam
                     (xnodes(Tbeams(e,2),2) - xnodes(Tbeams(e,1),2))^2 + ...
                     (xnodes(Tbeams(e,2),3) - xnodes(Tbeams(e,1),3))^2);  
    if Tmat_beams(e)==1
               %Section 1 
               a1=[2*data_mat(1,6),data_mat(1,6)+t1/2,t1];
               b1=[t1,t1,h1-t1];
               for i = 1:length(a1)
                   Ai(i) = a1(i)*b1(i);
                   Iyi(i) = a1(i)*b1(i)^3/12;
                   Izi(i) = b1(i)*a1(i)^3/12;
                   if i == 1
                       zi(i) = (h1-t1)/2 + t1/2;
                       yi(i) = 0;
                   elseif i == 2
                       zi(i) = -(h1-t1)/2 - t1/2;
                       yi(i) = -a/2;
                   else
                       zi(i) = 0;
                       yi(i) = 0;
                   end
               end
               % Compute the total area of the section
               Aa = sum(Ai);
               % Compute the centroid of the section
               yG = sum(yi.*Ai)/Aa;
               zG = sum(zi.*Ai)/Aa;
               % Apply the parallel axis theorem to compute the inertia of the section:
               Iya = sum(Iyi + Ai.*(zi-zG).*(zi-zG));
               Iza = sum(Izi + Ai.*(yi-yG).*(yi-yG));               
              Ja = 4*sum(min(Iyi, Izi));               
              V_beam(e) = Aa*L_beam(e);
               rho_beam(e)=data_mat(Tmat_beams(e),1);
               Mb(e)=rho_beam(e)*V_beam(e);
               Iy_beam(e)=Iya;
               Iz_beam(e)=Iza;
               J_beam(e)=Ja;
               E_beam(e)=data_mat(Tmat_beams(e),2);
               nu_beam(e)=data_mat(Tmat_beams(e),3);
               A_beam(e)=Aa;
    end
    if Tmat_beams(e)==2
              %Section 2 
               a2=[c+t2/2,b+t2/2,t2];
               b2=[t2,t2,h2-t2];
               for i = 1:length(a2)
                   Ai(i) = a2(i)*b2(i);
                   Iyi(i) = a2(i)*b2(i)^3/12;
                   Izi(i) = b2(i)*a2(i)^3/12;
                   if i == 1
                       zi(i) = (h2-t2)/2 + t2/2;
                       yi(i) = 0;
                   elseif i == 2
                       zi(i) = -(h2-t2)/2 - t2/2;
                       yi(i) = -a/2;
                   else
                       zi(i) = 0;
                       yi(i) = 0;
                   end
               end
               % Compute the total area of the section
               Ab = sum(Ai);
               % Compute the centroid of the section
               yG = sum(yi.*Ai)/Ab;
               zG = sum(zi.*Ai)/Ab;
               % Apply the parallel axis theorem to compute the inertia of the section:
               Iyb = sum(Iyi + Ai.*(zi-zG).*(zi-zG));
               Izb = sum(Izi + Ai.*(yi-yG).*(yi-yG));
               Jb = 4*sum(min(Iyi, Izi));                
                V_beam(e) = Ab*L_beam(e);
                rho_beam(e)=data_mat(Tmat_beams(e),1);
               Mb(e)=rho_beam(e)*V_beam(e);
               Iy_beam(e)=Iyb;
               Iz_beam(e)=Izb;
               J_beam(e)=Jb;
               E_beam(e)=data_mat(Tmat_beams(e),2);
               nu_beam(e)=data_mat(Tmat_beams(e),3);
               A_beam(e)=Ab;
    end
    
    if Tmat_beams(e)==3
              %Section 3
              a3=[2*d,0,t3];
              b3=[t3,t3,h3-t3];
               for i = 1:length(a3)
                   Ai(i) = a3(i)*b3(i);
                   Iyi(i) = a3(i)*b3(i)^3/12;
                   Izi(i) = b3(i)*a3(i)^3/12;
                   if i == 1
                       zi(i) = (h3-t3)/2 + t3/2;
                       yi(i) = 0;
                   elseif i == 2
                       zi(i) = -(h1-t1)/2 - t1/2;
                       yi(i) = -a/2;
                   else
                       zi(i) = 0;
                       yi(i) = 0;
                   end
               end
               % Compute the total area of the section
               Ac = sum(Ai);
               % Compute the centroid of the section
               yG = sum(yi.*Ai)/Ac;
               zG = sum(zi.*Ai)/Ac;
               % Apply the parallel axis theorem to compute the inertia of the section:
               Iyc = sum(Iyi + Ai.*(zi-zG).*(zi-zG));
               Izc = sum(Izi + Ai.*(yi-yG).*(yi-yG));
               Jc = 4*sum(min(Iyi, Izi));
               V_beam(e) = Ac*L_beam(e);
               rho_beam(e)=data_mat(3,1);
               Mb(e)=rho_beam(e)*V_beam(e);
               Iy_beam(e)=Iyc;
               Iz_beam(e)=Izc;
               J_beam(e)=Jc;
               E_beam(e)=data_mat(3,2);
               nu_beam(e)=data_mat(3,3);
               A_beam(e)=Ac;
    end    

end
M_beams=2*sum(Mb);
V_beams =2*sum(V_beam);

beam_parameters = struct('L',L_beam, 'A',A_beam, 'Iy',Iy_beam, 'Iz',Iz_beam,...
                             'J',J_beam, 'V1b',V_beam, 'm',m_beam, 'rho',rho_beam,...
                             'nu', nu_beam, 'E', E_beam,'M',M_beams,'V',V_beams,'Aa',Aa,...
                         'Iya',Iya,'Iza', Iza,  'Ja', Ja,'Ab',Ab,  'Iyb',Iyb,'Izb', Izb,  'Jb', Jb,'Ac',Ac,  'Iyc',Iyc,'Izc', Izc,  'Jc', Jc);

end

