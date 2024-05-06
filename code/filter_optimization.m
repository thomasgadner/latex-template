%% *********************************************************************
% Name: filter_optimization.m
% Date: 03-APR-2023
% Version: MATLAB R2021b
% **********************************************************************

%% *********************************************************************
% Clean workspace
% **********************************************************************
close all;
clear;
clc;

%% *********************************************************************
% Definition of constants
% **********************************************************************
T_amb = 273.15+25; % ambient temperature in K
T_max = 273.15+100; % maximum component temperature in K
V_dc = 300; % DC-link voltage in V
f_sw = 300e3; % switching frequency in Hz
I_nom = 33; % nominal DC current in A
I_ripple_pp = 2*0.10*I_nom; % peak-to-peak current ripple in A
I_peak = I_nom + I_ripple_pp/2; % peak current in A
L1_req = 2*V_dc/(4*f_sw*I_ripple_pp)*1e6; % req. filter inductance in uH
L1_req = L1_req/2; % due to filter symmetry reasons L1_req
C_up = 1; % penalty factor for rounding up of values
C_down = 1.5; % penalty factor for rounding down of values
num_series = 10; % number of acceptable series connections
L1_req = 78.56/2;

%% *********************************************************************
% Import of inductor database
% **********************************************************************
L_data = import_L_data('input_inductor_database.xlsx',...
    'input_inductor_database');
L_data = table2struct(L_data);
% preallocate additional fields
[L_data(:).L_90] = deal(NaN); % inductance value @90% of L_nom
[L_data(:).L_eff] = deal(NaN); % effective inductance (parallel conn)
[L_data(:).n_p] = deal(NaN); % number of parallel inductors
[L_data(:).C] = deal(NaN); % penalty factor for individual inductor
[L_data(:).Area] = deal(NaN); % inductor area
[L_data(:).Volume] = deal(NaN); % inductor volume
[L_data(:).T_op] = deal(NaN); % inductor operating temperature
% first iteration: determines number of paralleled inductors to satisfy
% current requirement; assigns corresponding part temperature
for idx = 1:size(L_data,1)
    % calculate inductance value 90% of L_nom
    L_data(idx).L_90 = 0.9*L_data(idx).L_nom;
    % calculate number of required parallel inductors
    n_p_tmp =  I_peak/L_data(idx).I_90;
    if n_p_tmp < 0.5
        L_data(idx).n_p = 1; % no parallel connection
        L_data(idx).L_eff = L_data(idx).L_90;
        L_data(idx).C = 1+C_up*abs(L_data(idx).n_p - n_p_tmp);
    else
        L_data(idx).n_p = round(n_p_tmp+0.25); % round to nearest int
        L_data(idx).L_eff = L_data(idx).L_90/L_data(idx).n_p;
        if(L_data(idx).n_p < n_p_tmp) % rounded down
            L_data(idx).C = 1+C_down*abs(L_data(idx).n_p - n_p_tmp);
        elseif(L_data(idx).n_p > n_p_tmp) % rounded up
            L_data(idx).C = 1+C_up*abs(L_data(idx).n_p - n_p_tmp);
        else
            L_data(idx).C = 0; % no panelty
        end
    end
    % assign corresponding part temperature
    if L_data(idx).n_p == 1
        L_data(idx).T_op = T_amb + L_data(idx).delta_T_TOT;
    else
        L_data(idx).T_op = T_amb + L_data(idx).delta_T_TOT_50;
    end
end
% second iteration: checks part temperature against maximum limit; 
% if yes, in case of single inductor the number of parallel connection
% is incremented; in case of already paralleled inductor, it is deleted
for idx = size(L_data,1):-1:1
    if L_data(idx).T_op > T_max
        if L_data(idx).n_p == 1
            % increment number of paralled inductors
            L_data(idx).n_p = 2;
            L_data(idx).L_eff = L_data(idx).L_90/L_data(idx).n_p;
            L_data(idx).C = ...
                1+C_up*abs(L_data(idx).n_p-I_peak/L_data(idx).I_90);
            L_data(idx).T_op = T_amb + L_data(idx).delta_T_TOT_50;
        else
            % delete entry
            L_data(idx) = [];
        end
    end
    % calculate occupied area
    L_data(idx).Area = L_data(idx).n_p*L_data(idx).Axial*...
        L_data(idx).Vertical;
    % calculate occupied volume
    L_data(idx).Volume = L_data(idx).n_p*L_data(idx).Axial*...
        L_data(idx).Vertical*L_data(idx).Height;
end
% optional: delete THT parts
for idx = size(L_data,1):-1:1
    if strcmp(L_data(idx).Technology,"THT")
        L_data(idx) = [];
    end
end
% generate zero element
for idx = 1:size(L_data,1)
    L_data(end+1) = struct('PartNumber',"",'Manufacturer',"",...
        'L_0',0,'L_nom',0,'I_90',0,'I_50',0,'Axial',0,'Vertical',0,...
        'Height',0,'T_op_max',Inf,'delta_T_TOT_50',0,...
        'delta_T_TOT',0,'delta_T_DC_50',0,'delta_T_DC',0,...
        'Technology',"",'L_90',0,'L_eff',0,'n_p',0,'C',0,'Area',0,...
        'Volume',0,'T_op',0); %#ok
end
clearvars n_p_tmp;

%% *********************************************************************
% Evaluation of results by series connection of same inductors
% **********************************************************************
for idx = 1:size(L_data,1)
    % part number
    simple_series(idx).PartNumber = L_data(idx).PartNumber; %#ok
    % number of parallel inductors
    simple_series(idx).n_p = L_data(idx).n_p; %#ok
    % number of series connected inductors
    simple_series(idx).n_s = round(L1_req/L_data(idx).L_eff); %#ok
    % number of required inductors
    simple_series(idx).n_req = simple_series(idx).n_p*...
        simple_series(idx).n_s; %#ok
    % obtained inductance
    simple_series(idx).L_res = ...
        simple_series(idx).n_s*L_data(idx).L_eff; %#ok
    % occupied area on PCB
    simple_series(idx).Area = ...
        imple_series(idx).n_s*L_data(idx).Area; %#ok
    % height of inductors
    simple_series(idx).Height = L_data(idx).Height; %#ok
    % temperature of inductors
    simple_series(idx).T_op = L_data(idx).T_op-273.15; %#ok
end
% optional: sort struct by height
simple_series = struct2table(simple_series); % convert to table
simple_series(strcmp(simple_series.PartNumber,""),:) = []; % delete 0
simple_series = sortrows(simple_series, 'Height'); % sort by 'Height'
simple_series = table2struct(simple_series); % change back to struct

%% *********************************************************************
% Optimization method using genetic algorithm (ga)
% **********************************************************************
% create an optimization problem with default properties
prob = optimproblem('ObjectiveSense','minimize');
% create optimization variable vector
L_idx = optimvar('L_idx',1,num_series,'Type','integer',...
    'LowerBound',1,'UpperBound',size(L_data,1));
for idx = 1:num_series
    L_s(idx) = fcn2optimexpr(@(x)L_data(x).L_eff,L_idx(idx)); %#ok
    C_s(idx) = fcn2optimexpr(@(x)L_data(x).C,L_idx(idx)); %#ok
    A_s(idx) = fcn2optimexpr(@(x)L_data(x).Area,L_idx(idx)); %#ok
    V_s(idx) = fcn2optimexpr(@(x)L_data(x).Volume,L_idx(idx)); %#ok
end
% create objective function
prob.Objective = sum(A_s.*C_s);
%prob.Objective = sum(V_s);
% calculate total inductance
L1_s_tot = sum(L_s);
% tolerance (+/-) of allowed inductance
L_tolerance = 3;
% create constraint
prob.Constraints.L1_ub = L1_s_tot <= (L1_req + L_tolerance);
prob.Constraints.L1_lp = L1_s_tot >= (L1_req - L_tolerance);
% population size
Populations_Size = 4000;
% number of elite childeren in %
Elite_Count = 15;
% number of crossover children and mutation children in %
Crossover_Children = 40;
% definition of options
options = optimoptions('ga',...
    'PlotFcn', @gaplotbestf,...
    'EliteCount',(Elite_Count/100)*Populations_Size,...
    'InitialPopulationRange',...
    [ones(1,num_series);size(L_data,1)*ones(1,num_series)],...
    'CrossoverFraction',(Crossover_Children/100),...
    'PopulationSize',Populations_Size);
% set the random number generator to the default seed
rng('default');
% solve problem
[sol,fval,exitflag,output,lambda] = ...
   solve(prob,'Solver','ga','Options',options);

%% *********************************************************************
% Post-processing
% **********************************************************************
% order solution vector in ascending order
sol_sorted = sort(sol.L_idx);
L1_sol = 0;
Area_sol = 0;
num_parts = 0;
% create solution struct
for idx = 1:size(sol_sorted,2)
    L_data_sol(idx).PartNumber = ...
        L_data(sol_sorted(idx)).PartNumber; %#ok
    L_data_sol(idx).Manufacturer = ...
        L_data(sol_sorted(idx)).Manufacturer; %#ok
    L_data_sol(idx).L_0 = L_data(sol_sorted(idx)).L_0; %#ok
    L_data_sol(idx).L_nom = L_data(sol_sorted(idx)).L_nom; %#ok
    L_data_sol(idx).n_p = L_data(sol_sorted(idx)).n_p; %#ok
    L_data_sol(idx).L_eff = L_data(sol_sorted(idx)).L_eff; %#ok
    L_data_sol(idx).Area = L_data(sol_sorted(idx)).Area; %#ok
    Area_sol = Area_sol + L_data_sol(idx).Area;
    L1_sol = L1_sol + L_data_sol(idx).L_eff;
    num_parts = num_parts + L_data_sol(idx).n_p;
end
disp(['Obtained inductance: ',num2str(L1_sol),' uH']);
disp(['Occupied PCB area: ',num2str(Area_sol),' mm^2']);
disp(['Required parts: ',num2str(num_parts)]);
% EOF