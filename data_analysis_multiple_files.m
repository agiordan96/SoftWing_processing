clc
clear variables
close all

format long

set(groot,'DefaultAxesFontSize',14)
set(groot,'DefaultLineLineWidth',1.5)

%% data reading

MyFolderInfo = dir('usable Data');

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

    if length(MyFolderInfo(k).name) > 11    % skipping "submerged" files 
        continue
    end

    exp_table = readtable("usable Data/" + MyFolderInfo(k).name, 'Delimiter', ', ', "Range", "D:I");

    exp_value.f_avg(k, :) = mean(exp_table{1:end, 1:3});  % average force vector for all of wing's config.
    exp_value.f_std(k, :) = std(exp_table{1:end, 1:3});   % standard on for each force component of every wing config.
    
    exp_value.t_avg(k, :) = mean(exp_table{1:end, 4:6}); % average torque vector for all of wing's config.
    exp_value.t_std(k, :) = std(exp_table{1:end, 4:6});  % standard deviation for each torque component of every wing config.

    exp_value.aoa(k) = str2double(MyFolderInfo(k).name(1:2));
    exp_value.vel(k, 1) = str2double(MyFolderInfo(k).name(4:5));
    exp_value.inflation(k) = str2double(MyFolderInfo(k).name(7));

end

remove = (exp_value.f_avg(:, 1) == 0);
fields = fieldnames(exp_value);
for k = 1:numel(fields)
    exp_value.(fields{k})(remove, :) = [];
end

clear exp_table

%% data sorting

T = struct2table(exp_value); % convert the struct array to a table
sortedT = sortrows(T, 'aoa'); % sort the table by 'aoa'
exp_value = table2struct(sortedT,'ToScalar',true); % convert the table back to the struct array

clear T sortedT

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
chord = 0.5; % m
rho = 1000; % kg / m^3 density of water
dyn_viscosity = 10^(-3); % Pa*s
kin_viscosity = dyn_viscosity / rho;

tor_transposed = zeros(length(exp_value.t_avg), 3);
tor_transposed(1:end, 1:3) = exp_value.t_avg(1:end, 1:3) + exp_value.f_avg(1:end, 1:3) * d;

%% data visualization: plots against AoA

close all
format short

sel_speed = [10, 15, 20, 25, 30, 40, 50];
sel_inflation = [0, 4, 1, 2, 3];

% CL / CD: presenting one plot per selected speed and all inflations, varying
% angle of attack

% mkdir('CL_over_CD_plot');
% 
% for j = 1:length(sel_speed)
% 
%     [status, msg, msgID] = mkdir(sprintf('../data_analysis/CL_over_CD_plot/ flow_speed_%d', sel_speed(j)));
%     rehash;  
%     clear k1 k2 k3 k4 k5
% 
%     figure('Position', [200, 200, 1000, 1000])
% 
%     title(['CL / CD plot # ', num2str(j), '; Flow Speed: ', num2str(sel_speed(j))],'fontweight','bold','fontsize', 24)
%     hold on
%     grid on
%     xlabel('AoA [deg]','fontweight','bold','fontsize', 20);
%     ylabel('CL [ ]','fontweight','bold','fontsize', 20);
%     xlim([-10 35])
%     ylim([-8.5*10^(-06) 6*10^(-06)])

%     for k = 1:length(exp_value.f_avg)
%         if k == 87 
%             continue
%         end
%          if (exp_value.vel(k) == sel_speed(j)) && (exp_value.inflation(k) == sel_inflation(1))
%             scatter(exp_value.aoa(k), exp_value.f_avg(k, 2) / exp_value.f_avg(k, 1), 'or', 'filled', 'LineWidth',5)
%             if exist('k1','var') == 1
%                 x_vec = [exp_value.aoa(k1), exp_value.aoa(k)];
%                 y_vec = [exp_value.f_avg(k1, 2) / exp_value.f_avg(k1, 1), exp_value.f_avg(k, 2) / exp_value.f_avg(k, 1)];
%                 plot(x_vec, y_vec, '--r')
%             end
%             k1 = k;
%          elseif (exp_value.vel(k) == sel_speed(j)) && (exp_value.inflation(k) == sel_inflation(2))
%              scatter(exp_value.aoa(k), exp_value.f_avg(k, 2) / exp_value.f_avg(k, 1), 'ok', 'filled')
%              if exist('k2','var') == 1
%                 x_vec = [exp_value.aoa(k2), exp_value.aoa(k)];
%                 y_vec = [exp_value.f_avg(k2, 2) / exp_value.f_avg(k2, 1), exp_value.f_avg(k, 2) / exp_value.f_avg(k, 1)];
%                 plot(x_vec, y_vec, '--k')
%             end
%             k2 = k;
%          elseif (exp_value.vel(k) == sel_speed(j)) && (exp_value.inflation(k) == sel_inflation(3))
%              scatter(exp_value.aoa(k), exp_value.f_avg(k, 2) / exp_value.f_avg(k, 1), 'om', 'filled')
%              if exist('k3','var') == 1
%                 x_vec = [exp_value.aoa(k3), exp_value.aoa(k)];
%                 y_vec = [exp_value.f_avg(k3, 2) / exp_value.f_avg(k3, 1), exp_value.f_avg(k, 2) / exp_value.f_avg(k, 1)];
%                 plot(x_vec, y_vec, '--m')
%             end
%             k3 = k;
%          elseif (exp_value.vel(k) == sel_speed(j)) && (exp_value.inflation(k) == sel_inflation(4))
%              scatter(exp_value.aoa(k), exp_value.f_avg(k, 2) / exp_value.f_avg(k, 1), 'ob', 'filled')
%              if exist('k4','var') == 1
%                 x_vec = [exp_value.aoa(k4), exp_value.aoa(k)];
%                 y_vec = [exp_value.f_avg(k4, 2) / exp_value.f_avg(k4, 1), exp_value.f_avg(k, 2) / exp_value.f_avg(k, 1)];
%                 plot(x_vec, y_vec, '--b')
%             end
%             k4 = k;
%          elseif (exp_value.vel(k) == sel_speed(j)) && (exp_value.inflation(k) == sel_inflation(5))
%             scatter(exp_value.aoa(k), exp_value.f_avg(k, 2) / exp_value.f_avg(k, 1), 'og', 'filled')
%             if exist('k5','var') == 1
%                 x_vec = [exp_value.aoa(k5), exp_value.aoa(k)];
%                 y_vec = [exp_value.f_avg(k5, 2) / exp_value.f_avg(k5, 1), exp_value.f_avg(k, 2) / exp_value.f_avg(k, 1)];
%                 plot(x_vec, y_vec, '--g')
%             end
%             k5 = k;
%          end
%     end
%     
%     legend({'inf. = 0 mL', 'inf. = 60 mL', 'inf. = 90 mL', 'inf. = 120 mL', 'inf. = 30 mL'}, ... 
%      'Location','north','Orientation','horizontal','fontsize', 16)
%     hold off
%     %saveas(gcf, ['../data_analysis/CL_over_CD_plot/flow_speed_', num2str(sel_speed(j)),'/CL_over_CD_plot_#', num2str(j), 'flow_speed', num2str(sel_speed(j))], 'svg');
% end



% CL: presenting one plot per selected speed and all inflations, varying
% angle of attack

for j = 1:length(sel_speed)

    [status, msg, msgID] = mkdir(sprintf('../pic/CL_plot/flow_speed_%d', sel_speed(j)));
    dyn_pressure = 0.5 * rho * sel_speed(j) ^ 2; % calculation of dynamic pressure
    div = dyn_pressure * S;
    Re = sel_speed(j) * chord / kin_viscosity; % Reynolds number
    
    clear k1 k2 k3 k4 k5

    figure('Position', [200, 200, 1000, 1000])

    title(['CL plot # ', num2str(j), '; Flow Speed: ', num2str(sel_speed(j))],'fontweight','bold','fontsize', 24)
    legend('Location','north','Orientation','horizontal','fontsize', 16)
    hold on
    grid on
    xlabel('AoA [deg]','fontweight','bold','fontsize', 20);
    ylabel('CL [E-06]','fontweight','bold','fontsize', 20);
    xlim([-10 35])
    ylim([-8.5 6])

    for k = 1:length(exp_value.f_avg)
        if k == 87 
            continue
        end
         if (exp_value.vel(k) == sel_speed(j)) && (exp_value.inflation(k) == sel_inflation(1))
            if exist('k1','var') == 0
                errorbar(exp_value.aoa(k), (exp_value.f_avg(k, 2) / div) * 10^(06), exp_value.f_std(k), 'or', 'DisplayName', 'inf. = 0 mL', 'CapSize', 18, 'MarkerFaceColor', 'r', 'LineWidth', 1, MarkerEdgeColor = 'red')
            else
                errorbar(exp_value.aoa(k), (exp_value.f_avg(k, 2) / div) * 10^(06), exp_value.f_std(k), 'or', 'CapSize', 18, 'MarkerFaceColor', 'r', 'LineWidth', 1, MarkerEdgeColor = 'red')
                x_vec = [exp_value.aoa(k1), exp_value.aoa(k)];
                y_vec = [(exp_value.f_avg(k1, 2) / div) * 10^(06), (exp_value.f_avg(k, 2) / div) * 10^(06)];
                plot(x_vec, y_vec, '--r')
            end
            k1 = k;
         elseif (exp_value.vel(k) == sel_speed(j)) && (exp_value.inflation(k) == sel_inflation(2))
             if exist('k2','var') == 0
                errorbar(exp_value.aoa(k), (exp_value.f_avg(k, 2) / div) * 10^(06), exp_value.f_std(k), 'ok', 'DisplayName', 'inf. = 30 mL', 'CapSize', 18, 'MarkerFaceColor', 'k', 'LineWidth', 1, MarkerEdgeColor = 'black')
             else
                errorbar(exp_value.aoa(k), (exp_value.f_avg(k, 2) / div) * 10^(06), exp_value.f_std(k), 'ok', 'CapSize', 18, 'MarkerFaceColor', 'k', 'LineWidth', 1, MarkerEdgeColor = 'black')
                x_vec = [exp_value.aoa(k2), exp_value.aoa(k)];
                y_vec = [(exp_value.f_avg(k2, 2) / div) * 10^(06), (exp_value.f_avg(k, 2) / div) * 10^(06)];
                plot(x_vec, y_vec, '--k')
            end
            k2 = k;
         elseif (exp_value.vel(k) == sel_speed(j)) && (exp_value.inflation(k) == sel_inflation(3))
             if exist('k3','var') == 0
                errorbar(exp_value.aoa(k), (exp_value.f_avg(k, 2) / div) * 10^(06), exp_value.f_std(k), 'om', 'DisplayName', 'inf. = 60 mL', 'CapSize', 18, 'MarkerFaceColor', 'm', 'LineWidth', 1, MarkerEdgeColor = 'magenta')
             else
                errorbar(exp_value.aoa(k), (exp_value.f_avg(k, 2) / div) * 10^(06), exp_value.f_std(k), 'om', 'CapSize', 18, 'MarkerFaceColor', 'm', 'LineWidth', 1, MarkerEdgeColor = 'magenta')
                x_vec = [exp_value.aoa(k3), exp_value.aoa(k)];
                y_vec = [(exp_value.f_avg(k3, 2) / div) * 10^(06), (exp_value.f_avg(k, 2) / div) * 10^(06)];
                plot(x_vec, y_vec, '--m')
            end
            k3 = k;
         elseif (exp_value.vel(k) == sel_speed(j)) && (exp_value.inflation(k) == sel_inflation(4))
             if exist('k4','var') == 0
                errorbar(exp_value.aoa(k), (exp_value.f_avg(k, 2) / div) * 10^(06), exp_value.f_std(k), 'ob', 'DisplayName', 'inf. = 90 mL', 'CapSize', 18, 'MarkerFaceColor', 'b', 'LineWidth', 1, MarkerEdgeColor = 'blue')
             else
                errorbar(exp_value.aoa(k), (exp_value.f_avg(k, 2) / div) * 10^(06), exp_value.f_std(k), 'ob', 'CapSize', 18, 'MarkerFaceColor', 'b', 'LineWidth', 1, MarkerEdgeColor = 'blue')
                x_vec = [exp_value.aoa(k4), exp_value.aoa(k)];
                y_vec = [(exp_value.f_avg(k4, 2) / div) * 10^(06), (exp_value.f_avg(k, 2) / div) * 10^(06)];
                plot(x_vec, y_vec, '--b')
            end
            k4 = k;
         elseif (exp_value.vel(k) == sel_speed(j)) && (exp_value.inflation(k) == sel_inflation(5))
            if exist('k5','var') == 0 
             errorbar(exp_value.aoa(k), (exp_value.f_avg(k, 2) / div) * 10^(06), exp_value.f_std(k), 'og', 'DisplayName', 'inf. = 120 mL', 'CapSize', 18, 'MarkerFaceColor', 'g', 'LineWidth', 1, MarkerEdgeColor = 'green')
            else
                errorbar(exp_value.aoa(k), (exp_value.f_avg(k, 2) / div) * 10^(06), exp_value.f_std(k), 'og', 'CapSize', 18, 'MarkerFaceColor', 'g', 'LineWidth', 1, MarkerEdgeColor = 'green')
                x_vec = [exp_value.aoa(k5), exp_value.aoa(k)];
                y_vec = [(exp_value.f_avg(k5, 2) / div) * 10^(06), (exp_value.f_avg(k, 2) / div) * 10^(06)];
                plot(x_vec, y_vec, '--g')
            end
            k5 = k;
         end
    end
    
    str_annotation = sprintf('Re = %.3e', Re);
    annotation('textbox', [0.696 0.77 0.1 0.1], 'String', str_annotation, ...
           'BackgroundColor','white','LineStyle','-','Fontsize', 16, 'Interpreter','latex' )
    hold off
    saveas(gcf, ['../data_analysis/CL_plot/', 'CL_plot_#', num2str(j), 'flow_speed', num2str(sel_speed(j))], 'svg');
    if ~isfolder('..')
        error('Corrupt or very very old file system, missing .. directory entry')
    elseif ~isfolder('../pic')
        error('No folder ../data_analysis')
    elseif ~isfolder('../pic/CL_plot')
        error('No folder ../pic/CL_plot')
    else
        fprintf('folder path ../pic/CL_plot/ is okay \n')
    end

end


% CD: presenting one plot per selected speed and all inflations, varying
% angle of attack

for j = 1:length(sel_speed)

    [status, msg, msgID] = mkdir(sprintf('../pic/CD_plot/flow_speed_%d', sel_speed(j)));
    dyn_pressure = 0.5 * rho * sel_speed(j) ^ 2; % calculation of dynamic pressure
    div = dyn_pressure * S;
    Re = sel_speed(j) * chord / kin_viscosity; % Reynolds number
    
    clear k1 k2 k3 k4 k5

    figure('Position', [200, 200, 1000, 1000])

    title(['CD plot # ', num2str(j), '; Flow Speed: ', num2str(sel_speed(j))], 'fontweight','bold','fontsize', 24)
    hold on
    grid on
    xlabel('AoA [deg]','fontweight','bold','fontsize', 20);
    ylabel('CD [E-06]','fontweight','bold','fontsize', 20);
    xlim([-10 35])
    ylim([-7 1])

    for k = 1:length(exp_value.f_avg)
        if k == 87 
            continue
        end
         if (exp_value.vel(k) == sel_speed(j)) && (exp_value.inflation(k) == sel_inflation(1))
            if exist('k1','var') == 0
                errorbar(exp_value.aoa(k), (exp_value.f_avg(k, 1) / div) * 10^(06), exp_value.f_std(k), 'or', 'DisplayName', 'inf. = 0 mL', 'CapSize', 18, 'MarkerFaceColor', 'r', 'LineWidth', 1, MarkerEdgeColor = 'red')
            else
                errorbar(exp_value.aoa(k), (exp_value.f_avg(k, 1) / div) * 10^(06), exp_value.f_std(k), 'or', 'CapSize', 18, 'MarkerFaceColor', 'r', 'LineWidth', 1, MarkerEdgeColor = 'red')
                x_vec = [exp_value.aoa(k1), exp_value.aoa(k)];
                y_vec = [(exp_value.f_avg(k1, 1) / div) * 10^(06), (exp_value.f_avg(k, 1) / div) * 10^(06)];
                plot(x_vec, y_vec, '--r')
            end
            k1 = k;
         elseif (exp_value.vel(k) == sel_speed(j)) && (exp_value.inflation(k) == sel_inflation(2))
             if exist('k2','var') == 0
                errorbar(exp_value.aoa(k), (exp_value.f_avg(k, 1) / div) * 10^(06), exp_value.f_std(k), 'ok', 'DisplayName', 'inf. = 30 mL', 'CapSize', 18, 'MarkerFaceColor', 'k', 'LineWidth', 1, MarkerEdgeColor = 'black')
             else
                errorbar(exp_value.aoa(k), (exp_value.f_avg(k, 1) / div) * 10^(06), exp_value.f_std(k), 'ok', 'CapSize', 18, 'MarkerFaceColor', 'k', 'LineWidth', 1, MarkerEdgeColor = 'black')
                x_vec = [exp_value.aoa(k2), exp_value.aoa(k)];
                y_vec = [(exp_value.f_avg(k2, 1)  / div) * 10^(06), (exp_value.f_avg(k, 1) / div) * 10^(06)];
                plot(x_vec, y_vec, '--k')
            end
            k2 = k;
         elseif (exp_value.vel(k) == sel_speed(j)) && (exp_value.inflation(k) == sel_inflation(3))
             if exist('k3','var') == 0
                errorbar(exp_value.aoa(k), (exp_value.f_avg(k, 1) / div) * 10^(06), exp_value.f_std(k), 'om', 'DisplayName', 'inf. = 60 mL', 'CapSize', 18, 'MarkerFaceColor', 'm', 'LineWidth', 1, MarkerEdgeColor = 'magenta')
             else
                 errorbar(exp_value.aoa(k), (exp_value.f_avg(k, 1) / div) * 10^(06), exp_value.f_std(k), 'ok', 'CapSize', 18, 'MarkerFaceColor', 'm', 'LineWidth', 1, MarkerEdgeColor = 'magenta')
                x_vec = [exp_value.aoa(k3), exp_value.aoa(k)];
                y_vec = [(exp_value.f_avg(k3, 1) / div) * 10^(06), (exp_value.f_avg(k, 1) / div) * 10^(06)];
                plot(x_vec, y_vec, '--m')
            end
            k3 = k;
         elseif (exp_value.vel(k) == sel_speed(j)) && (exp_value.inflation(k) == sel_inflation(4))
             if exist('k4','var') == 0
                errorbar(exp_value.aoa(k), (exp_value.f_avg(k, 1) / div) * 10^(06), exp_value.f_std(k), 'ob', 'DisplayName', 'inf. = 90 mL', 'CapSize', 18, 'MarkerFaceColor', 'b', 'LineWidth', 1, MarkerEdgeColor = 'blue')
             else
                 errorbar(exp_value.aoa(k), (exp_value.f_avg(k, 1) / div) * 10^(06), exp_value.f_std(k), 'ok', 'CapSize', 18, 'MarkerFaceColor', 'b', 'LineWidth', 1, MarkerEdgeColor = 'blue')
                x_vec = [exp_value.aoa(k4), exp_value.aoa(k)];
                y_vec = [(exp_value.f_avg(k4, 1) / div) * 10^(06), (exp_value.f_avg(k, 1) / div) * 10^(06)];
                plot(x_vec, y_vec, '--b')
            end
            k4 = k;
         elseif (exp_value.vel(k) == sel_speed(j)) && (exp_value.inflation(k) == sel_inflation(5))
             if exist('k5','var') == 0
                errorbar(exp_value.aoa(k), (exp_value.f_avg(k, 1) / div) * 10^(06), exp_value.f_std(k), 'og', 'DisplayName', 'inf. = 120 mL', 'CapSize', 18, 'MarkerFaceColor', 'g', 'LineWidth', 1, MarkerEdgeColor = 'green')
             else
                errorbar(exp_value.aoa(k), (exp_value.f_avg(k, 1) / div) * 10^(06), exp_value.f_std(k), 'og', 'CapSize', 18, 'MarkerFaceColor', 'g', 'LineWidth', 1, MarkerEdgeColor = 'green')
                x_vec = [exp_value.aoa(k5), exp_value.aoa(k)];
                y_vec = [(exp_value.f_avg(k5, 1) / div) * 10^(06), (exp_value.f_avg(k, 1) / div) * 10^(06)];
                plot(x_vec, y_vec, '--g')
            end
            k5 = k;
         end
    end
    
    legend({'inf. = 0 mL', 'inf. = 60 mL', 'inf. = 90 mL', 'inf. = 120 mL', 'inf. = 30 mL'}, ... 
     'Location','north','Orientation','horizontal','fontsize', 16)
    str_annotation = sprintf('Re = %.3e', Re);
    annotation('textbox', [0.696 0.77 0.1 0.1], 'String', str_annotation, ...
           'BackgroundColor','white','LineStyle','-','Fontsize', 16, 'Interpreter','latex' )
    hold off
    saveas(gcf, ['../pic/CD_plot/','/CD_plot_#', num2str(j), 'flow_speed', num2str(sel_speed(j))], 'svg');
    if ~isfolder('..')
        error('Corrupt or very very old file system, missing .. directory entry')
    elseif ~isfolder('../pic')
        error('No folder ../data_analysis')
    elseif ~isfolder('../pic/CD_plot')
        error('No folder ../pic/CD_plot')
    else
        fprintf('folder path ../pic/CD_plot/ is okay \n')
    end

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