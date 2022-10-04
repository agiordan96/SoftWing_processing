clc
clear variables
close all

format long

set(groot,'DefaultAxesFontSize',14)
set(groot,'DefaultLineLineWidth',1.5)

%% data read

MyFolderInfo = dir('usable Data');

force = struct; torque = struct;    % definition of force and torque related quantities' containers

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

    force_table = readtable("usable Data/" + MyFolderInfo(k).name, 'Delimiter', ', ', "Range", "D:F");
    torque_table = readtable("usable Data/" + MyFolderInfo(k).name, 'Delimiter', ', ', "Range", "G:I");   

    force.avg(k, :) = mean(force_table{1:end, :});  % average force vector for all of wing's config.
    torque.avg(k, 1:3) = mean(torque_table{1:end, :}); % average torque vector for all of wing's config.

    force.std(k, :) = std(force_table{1:end, :});   % standard on for each force component of every wing config.
    torque.std(k, :) = std(torque_table{1:end, :});  % standard deviation for each torque component of every wing config.

    if length(MyFolderInfo(k).name) == 11 || length(MyFolderInfo(k).name) > 12

        force.aoa(k) = str2double(MyFolderInfo(k).name(1:2));
        force.vel(k, 1) = str2double(MyFolderInfo(k).name(4:5));
        force.inflation(k) = str2double(MyFolderInfo(k).name(7));

        torque.aoa(k) = str2double(MyFolderInfo(k).name(1:2));
        torque.vel(k) = str2double(MyFolderInfo(k).name(4:5));
        torque.inflation(k) = str2double(MyFolderInfo(k).name(7));

    elseif length(MyFolderInfo(k).name) == 12

        force.aoa(k) = str2double(MyFolderInfo(k).name(2:3));
        force.vel(k, 1) = str2double(MyFolderInfo(k).name(5:6));
        force.inflation(k) = str2double(MyFolderInfo(k).name(8));

        torque.aoa(k) = str2double(MyFolderInfo(k).name(2:3));
        torque.vel(k) = str2double(MyFolderInfo(k).name(5:6));
        torque.inflation(k) = str2double(MyFolderInfo(k).name(8));

    end

end

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
tor_transposed(1:end, 1:3) = torque.avg(1:end, 1:3) + force.avg(1:end, 1:3) * d;

%% data visualization

sel_speed = [10, 15, 20, 30, 40, 50];
sel_inflation = [0, 1, 2, 3, 4];

% presenting one plot per selected speed and selected inflation, varying
% angle of attack

% for j = 1:length(sel_speed)
%     dyn_pressure = 0.5 * rho * sel_speed(j) ^ 2; % calculation of dynamic pressure
% 
%     for m = 1:length(sel_inflation)
%         figure
%         title(['Flow Speed: ', num2str(sel_speed(j)), ';   Inflation: ', num2str(sel_inflation(m)), ';'])
%         hold on
%         grid on
%         xlabel('AoA')
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
% 
% close all

% presenting one plot per selected speed and all inflations, varying
% angle of attack

for j = 1:length(sel_speed)
    dyn_pressure = 0.5 * rho * sel_speed(j) ^ 2; % calculation of dynamic pressure

    figure
    title(['Flow Speed: ', num2str(sel_speed(j)), ';   Inflation: ', num2str(sel_inflation(m)), ';'])
    hold on
    grid on
    xlabel('AoA')
    ylabel('CL / CD')
    xlim([-10 35])
    ylim([-10 10])


    for k = 1:length(MyFolderInfo)
         if (force.vel(k) == sel_speed(j)) && (force.inflation(k) == sel_inflation(1))
            scatter(force.aoa(k), force.avg(k, 2) / force.avg(k, 1), 'r', 'filled');
         elseif (force.vel(k) == sel_speed(j)) && (force.inflation(k) == sel_inflation(2))
             scatter(force.aoa(k), force.avg(k, 2) / force.avg(k, 1), 'o', 'filled');
         elseif (force.vel(k) == sel_speed(j)) && (force.inflation(k) == sel_inflation(3))
             scatter(force.aoa(k), force.avg(k, 2) / force.avg(k, 1), 'y', 'filled');
         elseif (force.vel(k) == sel_speed(j)) && (force.inflation(k) == sel_inflation(4))
             scatter(force.aoa(k), force.avg(k, 2) / force.avg(k, 1), 'b', 'filled');
         elseif (force.vel(k) == sel_speed(j)) && (force.inflation(k) == sel_inflation(5))
             scatter(force.aoa(k), force.avg(k, 2) / force.avg(k, 1), 'k', 'filled');
         end
    end
    
    legend({'inf. = 0 mL', 'inf. = 60 mL', 'inf. = 90 mL', 'inf. = 120 mL', 'inf. = 30 mL'}, ... 
     'Location','northwest','Orientation','horizontal')
    hold off

end


























