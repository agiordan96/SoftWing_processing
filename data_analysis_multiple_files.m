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

for k = 1:40 % for k = 1:length(MyFolderInfo)

    force_table = readtable("usable Data/" + MyFolderInfo(k).name, 'Delimiter', ', ', "Range", "D:F");
    torque_table = readtable("usable Data/" + MyFolderInfo(k).name, 'Delimiter', ', ', "Range", "G:I");

    avg_F(k, :) = mean(force_table{:, :});  % average force vector for all of wing's config.
    avg_T(k, :) = mean(torque_table{:, :}); % average torque vector for all of wing's config.

    std_F(k, :) = std(force_table{:, :});   % standard deviation for each force component of every wing config.
    std_T(k, :) = std(torque_table{:, :});  % standard deviation for each torque component of every wing config.

end

%% data processing

avg_F = mean(M_corrected(1:end, :, 4:6), 2); 
avg_T = mean(M_corrected(:, :, 7:9), 2);  

%% data visualization

figure 

grid on
hold on

plot(1:length(Fx), Fx, 'b');

plot(1:length(Fx), mean(Fx), '-r');


