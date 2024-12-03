% Define domain and discretization
L = 300;        % (mm) length of the domain
xh = 7e-1;      % (mm) spatial discretization
xmesh = 0:xh:L; % Discretized spatial domain
rng(234)
NSiml = 20000; % Number of synthetic patients
all_parameters = cell(NSiml, 1);
IW_at1 = zeros(NSiml, 1); 
IW_at2 = zeros(NSiml, 1); 
IW_at3 = zeros(NSiml, 1); 
IW = zeros(NSiml, 1); 
m1_values = zeros(NSiml, 5); % Random samples
m2_values = zeros(NSiml, 5); % Random samples
m1_rand_values = zeros(NSiml, 5);
m2_rand_values = zeros(NSiml, 5);
mass_80 = zeros(NSiml, 1);
mass_16 = zeros(NSiml, 1);
ki67 = zeros(NSiml, 1);
column_headers = {'param_b', 'param_d', 'param_h2', 'param_r', 'param_h3', 'param_Dm', ...
    'param_ts', 'param_S', 'IW', 'IW_after_3_months', ...
    'IW_after_6_months', 'IW_after_12_months', 'ki67', ...
    'm1_1', 'm1_2', 'm1_3', 'm1_4', 'm1_5', ...
    'm2_1', 'm2_2', 'm2_3', 'm2_4', 'm2_5', ...
    'm1_rand1', 'm1_rand2', 'm1_rand3', 'm1_rand4', 'm1_rand5', ...
    'm2_rand1', 'm2_rand2', 'm2_rand3', 'm2_rand4', 'm2_rand5', ...
    'T1Gd_MRI', 'FLAIR_MRI'};
all_data_for_csv = column_headers; % This will be your first row

% Main loop
for k = 1:NSiml
    % Load simulation results before and after resection
    simul_pdesys = load(strcat('G:/Data_npj_Journal/scene5/simulationResults_SimulID_', num2str(k), '.mat')); 
    simul_pdesys2 = load(strcat('G:/Data_npj_Journal/scene5/simulationResults2_SimulID_', num2str(k), '.mat'));
    simulationResults = simul_pdesys.simulationResults;
    results = simulationResults.Results;
    simulationResults2 = simul_pdesys2.simulationResults2;
    results2 = simulationResults2.Results;
    parameters = simulationResults.Parameters;
    % Calculate Ki67
    [ki67(k,:)] = solveForPI(parameters.param_b);
    % Calculate tumor masses
    Xmesh2 = simulationResults2.Parameters.Xmesh2;
    [mass_80(k,:), mass_16(k,:)] = func_Mass(results.p(end,:), xmesh);
    
    % Calculate infiltration widths
    [IW_at1(k,:)] = func_IW(results2.pr(3,:), Xmesh2);
    [IW_at2(k,:)] = func_IW(results2.pr(6,:), Xmesh2);
    [IW_at3(k,:)] = func_IW(results2.pr(12,:), Xmesh2);
    [IW(k,:)] = func_IW(results.p(end,:), xmesh);
    
    % Compute infiltration width values and select specific indices
    [m1_values(k,:), m2_values(k,:), m1_rand_values(k,:), m2_rand_values(k,:), ...
     ] = func_InfiltrationWidth(results.p(end,:), results.m1(end,:), results.m2(end,:), xmesh);
        

    % Combine all data
    param_values = [parameters.param_b, parameters.param_d, parameters.param_h2,...
        parameters.param_r, parameters.param_h3, parameters.param_Dm, parameters.param_ts, parameters.param_S];
    numeric_data_row = num2cell([IW(k,:), IW_at1(k,:), IW_at2(k,:), IW_at3(k,:), ki67(k,:), ...
        m1_values(k,:), m2_values(k,:), m1_rand_values(k,:), m2_rand_values(k,:), ...
        mass_80(k,:), mass_16(k,:)]);
    combined_row = [num2cell(param_values), numeric_data_row];
    all_data_for_csv = [all_data_for_csv; combined_row];
end

% Save results to CSV
writecell(all_data_for_csv, 'dataset_last4.csv');

% Function to calculate infiltration widths and specific samples
function [m1_values, m2_values,m1_rand_values, m2_rand_values] = func_InfiltrationWidth(vals, vals1, vals2, xmesh)
    ind99 = [];
    ind80 = [];
    ind60 = [];
    ind16 = [];
    ind05 = [];
    ind01 = [];

    for k = 1:length(vals)
        if(vals(1,k) >= (0.99 * max(vals)))
            ind99 = k;
        end
        if(vals(1,k) >= (0.80 * max(vals)))
            ind80 = k;
        end
        if(vals(1,k) >= (0.60 * max(vals)))
            ind60 = k;
        end
        if(vals(1,k) >= (0.16 * max(vals)))
            ind16 = k;
        end
        if(vals(1,k) >= (0.05 * max(vals)))
            ind05 = k;
        end
        if(vals(1,k) >= (0.01 * max(vals)))
            ind01 = k;
        end
    end
    inds = [ind99, ind80,ind60,ind16,ind05,ind01];

    random_indices = [];
    random_indices_2 = [];
    rng(123)
    for i = 1:5
       random_indices = [random_indices, randi([inds(i), inds(i+1)], 1)];
    end    

    for i = 1:5
       random_indices_2 = [random_indices_2, randi([inds(1), inds(4)], 1)];
    end    
 
        m1_values = vals1(end,round(xmesh(1,random_indices))); % Random values from m1
        m2_values = vals2(end,round(xmesh(1,random_indices))); % Random values from m2
        m1_rand_values = vals1(end,round(xmesh(1,random_indices_2)));
        m2_rand_values = vals2(end,round(xmesh(1,random_indices_2)));
end

% Function to calculate infiltration width
function [IW_at] = func_IW(vals, xmesh)
    for k = 1:length(vals)
        if(vals(1,k) >= (0.80 * max(vals)))
            ind80 = k;
        end
    end
    for k = length(vals):-1:1
        if(vals(1,k) <= (0.02 * max(vals)))
            ind02 = k;
        end
    end
    IW_at = xmesh(1,ind02) - xmesh(1,ind80);
end

% Function to calculate Ki67 biomarker
function [ki67] = solveForPI(a_val)
    syms a t ta AI PI
    equation = a_val *(1-AI-PI)-((1/t)*(PI+PI^2)-(1/ta)*AI*PI)==0;
    t_val = 24/24;
    ta_val = 8/24;
    AI_val = 0.7/100;
    equation_substituted = subs(equation, [a, t, ta, AI], [a_val, t_val, ta_val, AI_val]);
    solutions = solve(equation_substituted, PI);
    positive_solution = max(double(solutions));
    if positive_solution < 0
        error('No positive solution found.');
    else
        ki67 = positive_solution * 100; % Report as percentage
    end
end

% Function to calculate tumor masses
function [mass_80, mass_16] = func_Mass(p, xmesh)
    max_p = max(p);
    threshold_80 = 0.80 * max_p;
    threshold_16 = 0.16 * max_p;
    for t = 1:length(p)
        if(p(1,t) >= threshold_80)
            ind_80 = t;
        end
    end
    for t = 1:length(p)
        if(p(1,t) >= threshold_16)
            ind_16 = t;
        end
    end
    start_index = find(xmesh >= 0, 1);
    if isempty(start_index)
        start_index = 1; 
    end
    mass_80 = trapz(xmesh(start_index:ind_80), p(start_index:ind_80) / max_p);
    mass_16 = trapz(xmesh(start_index:ind_16), p(start_index:ind_16) / max_p);
end
