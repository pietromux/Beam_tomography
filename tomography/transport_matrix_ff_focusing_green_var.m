function [ R , R_all, z_all ,pos] = transport_matrix_ff_focusing_green_var(quads,gammaBeta,dbs)
% This function calculates the transport matrix from screen4 to the
% microscope objective screen as a function of quad4 and quad5 currents
% taking into account the real field maps (input by tanh function)
bRho=gammaBeta*1.70451e-3 ;
leff=0.0768;
qgrad4=0.45;
qgrad5=0.45;
qgrad6=0.45;

dz = 0.001;
R = eye(4);
b = 135;
Ioffset = 0.0;

pos_screen4 = 3.191; % position screen4
ds4 = 0.104; % screen4 to Q4
d45 = 0.086; % Center Q4 to center Q5 0.086
d56 = 0.085; % Center Q5 to center Q6 0.085
%dbs = 1.044;  % Center Q6 to Yag after box
%dbs = 0.476;  % Center Q6 to Yag in box
total_length = ds4 + d45 + d56 + dbs ;

%quad strengths. hysteresis is responsible for a constant gradient offset. 
k = [...
    (quads(1)+Ioffset)*qgrad4/bRho ...
    (quads(2)-Ioffset)*qgrad5/bRho ...
    (quads(3)+Ioffset)*qgrad6/bRho ...
    ];    

%quad positions
pos = pos_screen4 +[...
    ds4 ...
    ds4 + d45 ...
    ds4 + d45 + d56 ...
    ];

z_all = pos_screen4 + (dz:dz:total_length);
R_all = zeros(4,4,length(z_all)); i=1;
for z = z_all
%    gradquad = [...
%        k(1)/2*(tanh(b/2*(leff/2-(z-pos(1))))+tanh(b/2*(leff/2+(z-pos(1))))) ...
 %       k(2)/2*(tanh(b/2*(leff/2-(z-pos(2))))+tanh(b/2*(leff/2+(z-pos(2))))) ...
 %      ];
%    grad = sum(gradquad);
grad = k(1)/2*(tanh(b/2*(leff/2-(z-pos(1))))+tanh(b/2*(leff/2+(z-pos(1)))));
grad = grad + k(2)/2*(tanh(b/2*(leff/2-(z-pos(2))))+tanh(b/2*(leff/2+(z-pos(2)))));
grad = grad + k(3)/2*(tanh(b/2*(leff/2-(z-pos(3))))+tanh(b/2*(leff/2+(z-pos(3)))));

    if abs(grad) < 1e-4
        Rdelta = [1 dz 0 0;0 1 0 0;0 0 1 dz;0 0 0 1];
    else
        sK=sqrt(grad);
        sKL=sK*dz;
        Rdelta= [cos(sKL) sin(sKL)/sK 0 0;-sK*sin(sKL) cos(sKL) 0 0;0 0 cosh(sKL) sinh(sKL)/sK;0 0 sK*sinh(sKL) cosh(sKL)];
    end
    R = Rdelta*R;   
    R_all(:,:,i) = R;
    i=i+1;
end
end
