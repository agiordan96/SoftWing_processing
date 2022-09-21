clc
clear all
close all

format long

set(groot,'DefaultAxesFontSize',14)
set(groot,'DefaultLineLineWidth',1.5)

%% data read

MyFolderInfo = dir('usable Data');
for j = 1:length(MyFolderInfo)
    M(j, :, :) = readmatrix("usable Data/" + MyFolderInfo(j).name, "NumHeaderLines", 7, "Range", "A:I");
end

%% data correction

M_corrected(:, 1 : size(M, 2)-7, 1 : size(M, 3) - 1) = M(:, 8 : end, 1 : 9); % removed header and last column (time, NaN)

%% data processing

avg_F = mean(M_corrected(1:end, :, 4:6), 2); % average force vector for all of wing's config.
avg_T = mean(M_corrected(:, :, 7:9), 2); % average torque vector for all of wing's config. 

%% data visualization

figure 

grid on
hold on

plot(1:length(Fx), Fx, 'b');

plot(1:length(Fx), mean(Fx), '-r');


