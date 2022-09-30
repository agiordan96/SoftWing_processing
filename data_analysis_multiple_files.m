clc
clear all
close all

format long

set(groot,'DefaultAxesFontSize',14)
set(groot,'DefaultLineLineWidth',1.5)

%% data read

MyFolderInfo = dir('usable Data');

avg_F = zeros(length(MyFolderInfo), 3);
avg_T = zeros(length(MyFolderInfo), 3);

std_F = zeros(length(MyFolderInfo), 3);
std_T = zeros(length(MyFolderInfo), 3);

for k = 1:length(MyFolderInfo) 

    if MyFolderInfo(k).name == '.' 
        continue
    end

    if MyFolderInfo(k).name == ".DS_Store" 
        continue
    end

    force_table = readtable("usable Data/" + MyFolderInfo(k).name, 'Delimiter', ', ', "Range", "D:F");
    torque_table = readtable("usable Data/" + MyFolderInfo(k).name, 'Delimiter', ', ', "Range", "G:I");   

    avg_F(k, :) = mean(force_table{1:end, :});  % average force vector for all of wing's config.
    avg_T(k, :) = mean(torque_table{1:end, :}); % average torque vector for all of wing's config.

    std_F(k, :) = std(force_table{1:end, :});   % standard on for each force component of every wing config.
    std_T(k, :) = std(torque_table{1:end, :});  % standard deviation for each torque component of every wing config.

end

%% data processing

L = 0.5; 

% given a reference frame, Sf(x_s, y_s, z_s) represents the force sensing point,
% whereas St(x_t, y_t, z_t) represents the torque sensing point.

% given a reference frame, T(x_t, y_t, z_t) represents the transposition
% point of forces and momenta

% d(x, y, z) = (dx; dy; dz) vector defining distance between force and moment sensing
% point and transposition point

dx = 1;
dy = 0;
dz = L / 2;

d = [dx; dy; dz];

tor_transposed = zeros(length(MyFolderInfo), 3);
tor_transposed(1:end, 1:3) = avg_T(1:end, 1:3) + avg_F(1:end, 1:3) * d;


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


