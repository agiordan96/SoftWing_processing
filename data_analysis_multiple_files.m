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
force.aoa = zeros(length(MyFolderInfo));
force.vel = zeros(length(MyFolderInfo));
force.inflation = zeros(length(MyFolderInfo));

torque.avg = zeros(length(MyFolderInfo), 3);
torque.std = zeros(length(MyFolderInfo), 3);
torque.aoa = zeros(length(MyFolderInfo));
torque.vel = zeros(length(MyFolderInfo));
torque.inflation = zeros(length(MyFolderInfo));

for k = 1:length(MyFolderInfo) 

    if MyFolderInfo(k).name == '.' 
        continue
    end

    if MyFolderInfo(k).name == ".DS_Store" 
        continue
    end

    force_table = readtable("usable Data/" + MyFolderInfo(k).name, 'Delimiter', ', ', "Range", "D:F");
    torque_table = readtable("usable Data/" + MyFolderInfo(k).name, 'Delimiter', ', ', "Range", "G:I");   

    avg_F(k, 1:3) = mean(force_table{1:end, :});  % average force vector for all of wing's config.
    avg_T(k, 1:3) = mean(torque_table{1:end, :}); % average torque vector for all of wing's config.

    std_F(k, :) = std(force_table{1:end, :});   % standard on for each force component of every wing config.
    std_T(k, :) = std(torque_table{1:end, :});  % standard deviation for each torque component of every wing config.

    if length(MyFolderInfo(k).name) == 11

        avg_F(k, 4) = str2double(MyFolderInfo(k).name(1:2));
        avg_T(k, 4) = str2double(MyFolderInfo(k).name(1:2));

        avg_F(k, 5) = str2double(MyFolderInfo(k).name(1:2));
        avg_T(k, 5) = str2double(MyFolderInfo(k).name(1:2));

        avg_F(k, 6) = str2double(MyFolderInfo(k).name(1:2));
        avg_T(k, 6) = str2double(MyFolderInfo(k).name(1:2));
    end

end

%% data processing

L = 1;

% given a reference frame, Sf(x_s, y_s, z_s) represents the force sensing point,
% whereas St(x_t, y_t, z_t) represents the torque sensing point.

% given a reference frame, T(x_t, y_t, z_t) represents the transposition
% point of forces and momenta

% d(x, y, z) = (dx; dy; dz) vector defining distance between force and moment sensing
% point and transposition point

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

tor_transposed = zeros(length(MyFolderInfo), 3);
tor_transposed(1:end, 1:3) = avg_T(1:end, 1:3) + avg_F(1:end, 1:3) * d;

%flow_speed = string_flow_speed(MyFolderInfo, length(MyFolderInfo));



%% data visualization
%
% figure
% 
% grid on
% hold on
% 
% plot(1:length(Fx), Fx, 'b');
% 
% plot(1:length(Fx), mean(Fx), '-r');