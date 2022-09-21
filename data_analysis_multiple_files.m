clc
clear all
close all

set(groot,'DefaultAxesFontSize',14)
set(groot,'DefaultLineLineWidth',1.5)

MyFolderInfo = dir('usable Data');
for j = 1:length(MyFolderInfo)
    M(j, :, :) = readmatrix("usable Data/" + MyFolderInfo(1).name, "NumHeaderLines", 7, "Range", "A:I");
end

%%
M_corrected(:, 1 : size(M, 2)-7, 1 : size(M, 3) - 1) = M(:, 8 : end, 1 : 9); %removed header and last column







%% data read



%% data processing

avg_Fx = mean(Fx);
avg_Fy = mean(Fy);
avg_Fz = mean(Fz);
avg_Tx = mean(Tx);
avg_Ty = mean(Ty);
avg_Tz = mean(Tz);

%% data visualization

figure 

grid on
hold on

plot(1:length(Fx), Fx, 'b');


plot(1:length(Fx), mean(Fx), '-r');


