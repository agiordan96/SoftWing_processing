clc
clear all
close all

format long

set(groot,'DefaultAxesFontSize',14)
set(groot,'DefaultLineLineWidth',1.5)

%% data read

MyFolderInfo = dir('usable Data');

avg_F = zeros(3, length(MyFolderInfo));
avg_T = zeros(3, length(MyFolderInfo));

% for k = 1:length(MyFolderInfo)
for k = 1:40

    force_table = readtable("usable Data/" + MyFolderInfo(k).name, 'Delimiter', ', ', "Range", "D:F");
    torque_table = readtable("usable Data/" + MyFolderInfo(k).name, 'Delimiter', ', ', "Range", "G:I");

    avg_F(:, k) = mean(force_table{:, :}, 1);
    avg_T(:, k) = mean(torque_table{:, :}, 1);

end

%% data processing

avg_F = mean(M_corrected(1:end, :, 4:6), 2); % average force vector for all of wing's config.
avg_T = mean(M_corrected(:, :, 7:9), 2); % average torque vector for all of wing's config. 

%% data visualization

figure 

grid on
hold on

plot(1:length(Fx), Fx, 'b');

plot(1:length(Fx), mean(Fx), '-r');


