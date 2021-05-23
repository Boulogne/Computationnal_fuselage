function Ke = local_Ke_beams(E,A,L,Iy,Iz,J,v)

K1 = E*A/L;           % Axial load
K2 = 2*E*Iz/L^3;      % Shear load in y-direction and bending moment in z-direction
K3 = 2*E*Iy/L^3;      % Shear load in z-direction and bending moment in y-direction
K4 = E*J/(2*(1+v)*L); % Torsion moment in x-direction

Ke = zeros(12,12);    % initialize local stiffness matrix

% Axial load
Ke([1,7],[1,7]) = K1*[ 1  -1;
                      -1   1];
% Shear load in y-direction and bending moment in z-direction
Ke([2,6,8,12],[2,6,8,12]) = K2*[6      3*L     -6      3*L;
                                -3*L   2*L^2   -3*L    L^2;
                                -6     -3*L    6       -3*L^2;
                                3*L    L^2     -3*L    2*L^2];

% Shear load in z-direction and bending moment in y-direction
Ke([3,5,9,11],[3,5,9,11]) = K3*[6       -3*L    -6     -3*L;
                                -3*L    2*L^2   3*L    L^2;
                                -6      3*L     6      3*L;
                                -3*L    L^2     3*L    2*L^2];
% Torsion moment in x-direction
Ke([4,10],[4,10]) = K4*[ 1  -1;
                        -1   1];


end