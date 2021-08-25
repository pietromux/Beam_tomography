clear all
close all
%%
G=7; %gamma factor
gammaBeta = sqrt(G^2-1);
dbs = 1.1;
Ip = [1.0, -2.0, 1.0];
current_settings = [{'I3', 'I4', 'I5', 'actual x rotation', 'actual y rotation', 'goal x rotation', 'goal y rotation'}];
%%
for i = 0:0.1:3.14
    for j = 0:0.1:3.14
        rx_goal = i;
        ry_goal = j;
        fun=@(quads)rotation_cost_function(quads, rx_goal, ry_goal, gammaBeta, dbs);
%         options = optimset('Display','iter');
        options = optimoptions ('fsolve', 'Display', 'none','Algorithm', 'levenberg-marquardt');
        options.MaxIter = 50000 ;
        options.MaxFunEvals = 5000000 ;
        options.FunctionTolerance = 1e-50;
        options.OptimalityTolerance = 1e-50;
        options.StepTolerance = 1e-25;
%         options = optimoptions(@lsqnonlin,options);
        [Is,fval]=fsolve(fun,Ip,options);
%         lb = [-6.0, -6.0, -6.0];
%         ub = [6.0, 6.0, 6.0];
%         [Is,fval]=lsqnonlin(fun,Ip,lb,ub,options);
        Is = round(Is,2);
        [shearx, sheary, e1, e3, rotationsx, rotationsy] = SER_decomposition(Is, gammaBeta, dbs);        
        format = '\nguess currents %.4f_%.4f_%.4f: rotX = %.4f, rotY = %.4f. Goal: rx = %.4f, ry = %.4f';
        fprintf(format, Is(1), Is(2), Is(3), rotationsx, rotationsy, rx_goal, ry_goal)
        current_settings = [current_settings; {Is(1), Is(2), Is(3), rotationsx, rotationsy, rx_goal, ry_goal}];
    end
end
%%
writecell(current_settings, 'D:\UCLA\PEGASUS\tomography\current_settings_G7_no_constraint.xls')
%%
current_settings = readtable('current_settings_2.xls');
%%
current_settings_cell_array = current_settings;
current_settings_cell_array(1,:) = [];
current_settings_matrix = table2array(current_settings_cell_array);
% current_settings_matrix = cell2mat(current_settings_cell_array);
%%
figure
scatter3(current_settings_matrix(:,4),current_settings_matrix(:,5),current_settings_matrix(:,1), 8.5, 'filled', 'blue','MarkerFaceAlpha',.5)
hold on
scatter3(current_settings_matrix(:,4),current_settings_matrix(:,5),current_settings_matrix(:,2), 8.5, 'filled', 'red', 'MarkerFaceAlpha',.5)
hold on
scatter3(current_settings_matrix(:,4),current_settings_matrix(:,5),current_settings_matrix(:,3), 8.5, 'filled', 'cyan','MarkerFaceAlpha',.6)
xlim([0 3.2])
ylim([0 3.2])
patch([[-3.2 -3.2] [3.2 3.2]], [[-3.2 3.2] [3.2 -3.2]], [1 1 1 1]*0, 'k', 'FaceAlpha',0.2);
legend('I3', 'I4', 'I5')
ax = gca;
ax.BoxStyle = 'full';
xlabel('X Rotation (rad)');
ylabel('Y Rotation (rad)');
zlabel('Current Value');
title('Current Settings vs. Rotations');
%%
figure
scatter(current_settings_matrix(:,4),current_settings_matrix(:,5), 8.5, 'filled', 'blue','MarkerFaceAlpha',.5)
hold on
scatter(current_settings_matrix(:,6),current_settings_matrix(:,7), 8.5, 'filled', 'red', 'MarkerFaceAlpha',.5)
legend('Actual Rotations', 'Goal Rotations')
xlabel('X Rotation (rad)');
ylabel('Y Rotation (rad)');
ylim([0 3.3])
title('Current Setting Optimization');
%%
figure
scatter3(current_settings_matrix(:,4),current_settings_matrix(:,5),current_settings_matrix(:,6)-current_settings_matrix(:,4), 8.5, 'filled', 'blue','MarkerFaceAlpha',.5)
hold on
scatter3(current_settings_matrix(:,4),current_settings_matrix(:,5), current_settings_matrix(:,7)-current_settings_matrix(:,5),8.5, 'filled', 'red', 'MarkerFaceAlpha',.5)
hold on
xlim([0 3.2])
ylim([0 3.2])
patch([[-3.2 -3.2] [3.2 3.2]], [[-3.2 3.2] [3.2 -3.2]], [1 1 1 1]*0, 'k', 'FaceAlpha',0.3);
legend('X Rotation Goal Difference', 'Y Rotation Goal Difference')
ax = gca;
ax.BoxStyle = 'full';
xlabel('X Rotation (rad)');
ylabel('Y Rotation (rad)');
zlabel('Difference (rad)');
title('Difference between Goal Rotation and Actual Rotation Achieved for Optimized Current Values');