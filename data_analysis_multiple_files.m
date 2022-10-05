clc
clear variables
close all

format long

set(groot,'DefaultAxesFontSize',14)
set(groot,'DefaultLineLineWidth',1.5)

%% data reading

MyFolderInfo = dir('usable Data full');

force = struct; torque = struct;    % definition of force and torque related quantities' containers
exp_value  = struct;

exp_value.f_avg = zeros(length(MyFolderInfo), 3);
exp_value.f_std = zeros(length(MyFolderInfo), 3);
exp_value.t_avg = zeros(length(MyFolderInfo), 3);
exp_value.t_std = zeros(length(MyFolderInfo), 3);
exp_value.aoa = zeros(length(MyFolderInfo), 1);
exp_value.vel = zeros(length(MyFolderInfo), 1);
exp_value.inflation = zeros(length(MyFolderInfo), 1);

force.avg = zeros(length(MyFolderInfo), 3);
force.std = zeros(length(MyFolderInfo), 3);
force.aoa = zeros(length(MyFolderInfo), 1);
force.vel = zeros(length(MyFolderInfo), 1);
force.inflation = zeros(length(MyFolderInfo), 1);

torque.avg = zeros(length(MyFolderInfo), 3);
torque.std = zeros(length(MyFolderInfo), 3);
torque.aoa = zeros(length(MyFolderInfo), 1);
torque.vel = zeros(length(MyFolderInfo), 1);
torque.inflation = zeros(length(MyFolderInfo), 1);


for k = 1:length(MyFolderInfo) 

    if MyFolderInfo(k).name == '.' 
        continue
    end

    if MyFolderInfo(k).name == ".DS_Store" 
        continue
    end

    exp_table = readtable("usable Data full/" + MyFolderInfo(k).name, 'Delimiter', ', ', "Range", "D:I");

    exp_value.f_avg(k, :) = mean(exp_table{1:end, 1:3});  % average force vector for all of wing's config.
    exp_value.f_std(k, :) = std(exp_table{1:end, 1:3});   % standard on for each force component of every wing config.
    
    exp_value.t_avg(k, :) = mean(exp_table{1:end, 4:6}); % average torque vector for all of wing's config.
    exp_value.t_std(k, :) = std(exp_table{1:end, 4:6});  % standard deviation for each torque component of every wing config.

    exp_value.aoa(k) = str2double(MyFolderInfo(k).name(1:2));
    exp_value.vel(k, 1) = str2double(MyFolderInfo(k).name(4:5));
    exp_value.inflation(k) = str2double(MyFolderInfo(k).name(7));

end

%% data sorting

T = struct2table(exp_value); % convert the struct array to a table
sortedT = sortrows(T, 'aoa'); % sort the table by 'aoa'
exp_value = table2struct(sortedT); % convert the table back to the struct array

%% data processing

L = 1;

% given a reference frame, Sf(x_s, y_s, z_s) represents the force sensing point,
% whereas St(x_t, y_t, z_t) represents the torque sensing point.

% given a reference frame, T(x_t, y_t, z_t) represents the transposition
% point of forces and momenta

% d(x, y, z) = (dx; dy; dz) vector defining distance between force and moment sensing
% points and transposition point (one vector because force and momenta sensing
% points coincide)

% CL = L / (dynamic pressure * wing_surface)
% CD = D / (dynamic pressure * wing_surface)


x_s = 1;
y_s = 0;
z_s = L / 2;

x_t = 0;
y_t = 0;
z_t = 0;

dx = x_t - x_s;
dy = y_t - y_s;
dz = z_t - z_s;

d = [dx; dy; dz];

S = 1; % m^2, wing's surface
rho = 1000; % kg / m^3 density of water

tor_transposed = zeros(length(MyFolderInfo), 3);
tor_transposed(1:end, 1:3) = exp_value.t_avg(1:end, 1:3) + exp_value.f_avg(1:end, 1:3) * d;

%% data visualization

close all

% plots against AoA

sel_speed = [10, 15, 20, 30, 40, 50];
sel_inflation = [0, 1, 2, 3, 4];

% CD / CL: presenting one plot per selected speed and selected inflation, varying
% angle of attack

% for j = 1:length(sel_speed)
%     dyn_pressure = 0.5 * rho * sel_speed(j) ^ 2; % calculation of dynamic pressure
% 
%     for m = 1:length(sel_inflation)
%         figure
%         title(['Flow Speed: ', num2str(sel_speed(j)), ';   Inflation: ', num2str(sel_inflation(m)), ';'])
%         hold on
%         grid on
%         xlabel('AoA [deg]')
%         ylabel('CL/CD')
%         xlim([-10 35])
%         ylim([-10 10])
% 
%         for k = 1:length(MyFolderInfo)
%             if (force.vel(k) == sel_speed(j)) && (force.inflation(k) == sel_inflation(m))
%                 scatter(force.aoa(k), force.avg(k, 2) / force.avg(k, 1), 'r', 'filled');
%             end
%         end
% 
%         hold off
%     end
% end



% CD / CL: presenting one plot per selected speed and all inflations, varying
% angle of attack

% for j = 1:length(sel_speed)
%     dyn_pressure = 0.5 * rho * sel_speed(j) ^ 2; % calculation of dynamic pressure
% 
%     figure
%     title(['CD/CL plot # ', num2str(j), '; Flow Speed: ', num2str(sel_speed(j))])
%     hold on
%     grid on
%     xlabel('AoA [deg]')
%     ylabel('CL / CD')
%     xlim([-10 35])
%     ylim([-10 10])
% 
% 
%     for k = 1:length(MyFolderInfo)
%          if (force.vel(k) == sel_speed(j)) && (force.inflation(k) == sel_inflation(1))
%             scatter(force.aoa(k), force.avg(k, 2) / force.avg(k, 1), 'r', 'filled');
% %             plot(force.aoa(k), force.avg(k, 2) / force.avg(k, 1), '-r');
%          elseif (force.vel(k) == sel_speed(j)) && (force.inflation(k) == sel_inflation(2))
%              scatter(force.aoa(k), force.avg(k, 2) / force.avg(k, 1), 'o', 'filled');
%          elseif (force.vel(k) == sel_speed(j)) && (force.inflation(k) == sel_inflation(3))
%              scatter(force.aoa(k), force.avg(k, 2) / force.avg(k, 1), 'y', 'filled');
%          elseif (force.vel(k) == sel_speed(j)) && (force.inflation(k) == sel_inflation(4))
%              scatter(force.aoa(k), force.avg(k, 2) / force.avg(k, 1), 'b', 'filled');
%          elseif (force.vel(k) == sel_speed(j)) && (force.inflation(k) == sel_inflation(5))
%              scatter(force.aoa(k), force.avg(k, 2) / force.avg(k, 1), 'k', 'filled');
%          end
%     end
%     
%     legend({'inf. = 0 mL', 'inf. = 60 mL', 'inf. = 90 mL', 'inf. = 120 mL', 'inf. = 30 mL'}, ... 
%      'Location','northwest','Orientation','horizontal')
%     hold off
% 
% end



% CL: presenting one plot per selected speed and all inflations, varying
% angle of attack

% for j = 1:length(sel_speed)
%     dyn_pressure = 0.5 * rho * sel_speed(j) ^ 2; % calculation of dynamic pressure
% 
%     figure
%     title(['CL plot # ', num2str(j), '; Flow Speed: ', num2str(sel_speed(j))])
%     hold on
%     grid on
%     xlabel('AoA [deg]')
%     ylabel('CL / CD')
%     xlim([-10 35])
%     ylim([-10 10])
% 
%     for k = 1:length(MyFolderInfo)
%          if (force.vel(k) == sel_speed(j)) && (force.inflation(k) == sel_inflation(1))
%             scatter(force.aoa(k), force.avg(k, 2), 'r', 'filled');
% %             plot(force.aoa(k), force.avg(k, 2) / force.avg(k, 1), '-r');
%          elseif (force.vel(k) == sel_speed(j)) && (force.inflation(k) == sel_inflation(2))
%              scatter(force.aoa(k), force.avg(k, 2), 'o', 'filled');
%          elseif (force.vel(k) == sel_speed(j)) && (force.inflation(k) == sel_inflation(3))
%              scatter(force.aoa(k), force.avg(k, 2), 'y', 'filled');
%          elseif (force.vel(k) == sel_speed(j)) && (force.inflation(k) == sel_inflation(4))
%              scatter(force.aoa(k), force.avg(k, 2), 'b', 'filled');
%          elseif (force.vel(k) == sel_speed(j)) && (force.inflation(k) == sel_inflation(5))
%              scatter(force.aoa(k), force.avg(k, 2), 'k', 'filled');
%          end
%     end
%     
%     legend({'inf. = 0 mL', 'inf. = 60 mL', 'inf. = 90 mL', 'inf. = 120 mL', 'inf. = 30 mL'}, ... 
%      'Location','northwest','Orientation','horizontal')
%     hold off
% 
% end
% 


% CD: presenting one plot per selected speed and all inflations, varying
% angle of attack


for j = 1:length(sel_speed)
    dyn_pressure = 0.5 * rho * sel_speed(j) ^ 2; % calculation of dynamic pressure
    clear k1

    figure
    title(['CD plot # ', num2str(j), '; Flow Speed: ', num2str(sel_speed(j))])
    hold on
    grid on
    xlabel('AoA [deg]')
    ylabel('CL / CD')
    xlim([-10 35])
    ylim([-10 10])

    for k = 1:length(MyFolderInfo)
        if k == 87 
            continue
        end
         if (force.vel(k) == sel_speed(j)) && (force.inflation(k) == sel_inflation(1))
            if exist('k1','var') == 1
                scatter(force.aoa(k), force.avg(k, 1), 'r', 'filled');
                x_vec = [force.aoa(k1), force.aoa(k)];
                y_vec = [force.avg(k1, 1), force.avg(k, 1)];
                plot(x_vec, y_vec, '-or')
            else
                scatter(force.aoa(k), force.avg(k, 1), 'r', 'filled');
            end
            k1 = k;
         elseif (force.vel(k) == sel_speed(j)) && (force.inflation(k) == sel_inflation(2))
             scatter(force.aoa(k), force.avg(k, 1), 'o', 'filled');
         elseif (force.vel(k) == sel_speed(j)) && (force.inflation(k) == sel_inflation(3))
             scatter(force.aoa(k), force.avg(k, 1), 'y', 'filled');
         elseif (force.vel(k) == sel_speed(j)) && (force.inflation(k) == sel_inflation(4))
             scatter(force.aoa(k), force.avg(k, 1), 'b', 'filled');
         elseif (force.vel(k) == sel_speed(j)) && (force.inflation(k) == sel_inflation(5))
             scatter(force.aoa(k), force.avg(k, 1), 'k', 'filled');
         end
    end
    
    legend({'inf. = 0 mL', 'inf. = 60 mL', 'inf. = 90 mL', 'inf. = 120 mL', 'inf. = 30 mL'}, ... 
     'Location','northwest','Orientation','horizontal')
    hold off

end



% CM pitching: presenting one plot per selected speed and all inflations, varying
% angle of attack

% for j = 1:length(sel_speed)
%     dyn_pressure = 0.5 * rho * sel_speed(j) ^ 2; % calculation of dynamic pressure
% 
%     figure
%     title(['CM pitch plot # ', num2str(j), '; Flow Speed: ', num2str(sel_speed(j))])
%     hold on
%     grid on
%     xlabel('AoA [deg]')
%     ylabel('CL / CD')
%     xlim([-10 35])
%     ylim([-10 10])
% 
%     for k = 1:length(MyFolderInfo)
%          if (torque.vel(k) == sel_speed(j)) && (torque.inflation(k) == sel_inflation(1))
%             scatter(torque.aoa(k), torque.avg(k, 2), 'r', 'filled');
%          elseif (torque.vel(k) == sel_speed(j)) && (torque.inflation(k) == sel_inflation(2))
%              scatter(torque.aoa(k), torque.avg(k, 2), 'o', 'filled');
%          elseif (torque.vel(k) == sel_speed(j)) && (torque.inflation(k) == sel_inflation(3))
%              scatter(torque.aoa(k), torque.avg(k, 2), 'y', 'filled');
%          elseif (torque.vel(k) == sel_speed(j)) && (torque.inflation(k) == sel_inflation(4))
%              scatter(torque.aoa(k), torque.avg(k, 2), 'b', 'filled');
%          elseif (torque.vel(k) == sel_speed(j)) && (torque.inflation(k) == sel_inflation(5))
%              scatter(torque.aoa(k), torque.avg(k, 2), 'k', 'filled');
%          end
%     end
%     
%     legend({'inf. = 0 mL', 'inf. = 60 mL', 'inf. = 90 mL', 'inf. = 120 mL', 'inf. = 30 mL'}, ... 
%      'Location','northwest','Orientation','horizontal')
%     hold off
% 
% end
% 


