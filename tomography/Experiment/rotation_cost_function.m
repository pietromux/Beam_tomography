function [diff] = rotation_cost_function(quads, rx_goal, ry_goal, gammaBeta, dbs)
[ R ,~, ~ ,~] = transport_matrix_ff_focusing_green_var(-quads,gammaBeta,dbs);
S=[1 0 0 0; 
  (R(1,2)*R(2,2)+R(1,1)*R(2,1))/(R(1,1)^2+R(1,2)^2) 1 0 0;
   0 0 1 0;
   0 0 (R(3,4)*R(4,4)+R(3,3)*R(4,3))/(R(3,3)^2+R(3,4)^2) 1];
E=[sqrt(R(1,1)^2+R(1,2)^2) 0 0 0;
    0 1/sqrt(R(1,1)^2+R(1,2)^2) 0 0;
    0 0 sqrt(R(3,3)^2+R(3,4)^2) 0;
    0 0 0 1/sqrt(R(3,3)^2+R(3,4)^2)];
Rot=[R(1,1)/sqrt(R(1,1)^2+R(1,2)^2) R(1,2)/sqrt(R(1,1)^2+R(1,2)^2) 0 0;
    -R(1,2)/sqrt(R(1,1)^2+R(1,2)^2) R(1,1)/sqrt(R(1,1)^2+R(1,2)^2) 0 0; ...
    0 0 R(3,3)/sqrt(R(3,3)^2+R(3,4)^2) R(3,4)/sqrt(R(3,3)^2+R(3,4)^2);
    0 0 -R(3,4)/sqrt(R(3,3)^2+R(3,4)^2) R(3,3)/sqrt(R(3,3)^2+R(3,4)^2)];
e1=E(1,1);
e3=E(3,3);
shearx=S(2,1);
sheary=S(4,3);
rotationsx=acos(Rot(1,1));
rotationsy=acos(Rot(3,3));
diff_rx = rx_goal - rotationsx;
diff_ry = ry_goal - rotationsy;
diff = [diff_rx, diff_ry];