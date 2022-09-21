clc
clear all
close all

set(groot,'DefaultAxesFontSize',14)
set(groot,'DefaultLineLineWidth',1.5)

%% data read

T = readtable("usable Data/05-40-4.csv");



status = table2array(T(:, 1));
RDT_seq = table2array(T(:, 2));
FT_seq = table2array(T(:, 3));
Fx = table2array(T(:, 4));
Fy = table2array(T(:, 5));
Fz = table2array(T(:, 6));
Tx = table2array(T(:, 7));
Ty = table2array(T(:, 8));
Tz = table2array(T(:, 9));
rec_time = table2array(T(:, 10));

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


