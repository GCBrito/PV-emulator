clc; clear; close all; format long
%% Module Parameters

ns      = 60;       % number of series cells
np      = 1;        % number of parallel branches

Vmp_mod_ref  = 30.1;     % voltage at maximum power point (V)
Imp_mod_ref  = 8.30;     % current at maximum power point (A)
Voc_mod_ref  = 37.2;     % open-circuit voltage (V)
Isc_mod_ref  = 8.87;     % short-circuit current (A)

Tref = 25 + 273.15; % Reference temperature (K)
Sref = 1000; % Reference irradiance (W/m²)

alpha   = 0.00065;  % temperature coefficient of Isc (%/K)

%% Operating Conditions

T = 44.5 + 273.15; % Current temperature (K)
S = 765; % Current irradiance (W/m²)

Vmp_cell_ref = Vmp_mod_ref/ns;
Imp_cell_ref = Imp_mod_ref/np;
Voc_cell_ref = Voc_mod_ref/ns;
Isc_cell_ref = Isc_mod_ref/np;

%% Physical Constants

q = 1.60217662e-19; % Elementary charge (C)
k = 1.38064852e-23; % Boltzmann constant (J/K)

E_G0 = 1.166;            % Band gap energy at 0K (eV)
k1   = 4.73e-4;          % Coefficient k1 (eV/K)
k2   = 636;              % Coefficient k2 (K)

%% Parameter Estimation via fsolve

x0 = [Isc_cell_ref; 1e-9; 3; 0.01; 2];
opts = optimoptions('fsolve', ...
    'Display', 'off', 'TolFun', 1e-10, 'TolX', 1e-10, 'MaxIter', 1000, 'MaxFunctionEvaluations', 2000);
fun = @(x) residuals_2_20(x, Voc_cell_ref, Isc_cell_ref, Vmp_cell_ref, Imp_cell_ref, q, k, Tref);
[xsol, ~, exitflag] = fsolve(fun, x0, opts);
if exitflag <= 0
    warning('fsolve did not converge (exitflag = %d)', exitflag);
end

Iph_ref = xsol(1);
Is0_ref = xsol(2);
A = xsol(3);
Rs = xsol(4);
Rp = xsol(5);

fprintf('Iph_ref = %.6f A\n', Iph_ref);
fprintf('Is0_ref = %.2e A\n', Is0_ref);
fprintf('A       = %.6f\n', A);
fprintf('Rs      = %.6f Ohms\n', Rs);
fprintf('Rp      = %.6f Ohms\n', Rp);

%% Characteristic Points (used to build the I-V model)

V_cell = [
    0.0;
    0.2 * Vmp_mod_ref / ns;
    0.4 * Vmp_mod_ref / ns;
    0.6 * Vmp_mod_ref / ns;
    0.8 * Vmp_mod_ref / ns;
    0.9 * Vmp_mod_ref / ns;
    Vmp_mod_ref / ns;
    (Vmp_mod_ref / ns + Voc_mod_ref / ns) / 2.0;
    0.9 * Voc_mod_ref / ns;
    0.95 * Voc_mod_ref / ns;
    Voc_mod_ref / ns];


I_cell = arrayfun(@(V) solve_I_V_2_11(V, Iph_ref, Is0_ref, A, Rs, Rp, ...
    q, k, S, Sref, alpha, T, Tref, E_G0, k1, k2), V_cell);
points_V = V_cell * ns;
points_I = I_cell * np;

%% Piecewise PV Model

% Coefficient initialization
coeffs = zeros(length(points_V) - 1, 2);
for i = 1:length(points_V) - 1
    a = (points_I(i + 1) - points_I(i)) / (points_V(i + 1) - points_V(i));
    b = points_I(i) - a * points_V(i);
    coeffs(i, :) = [a, b]; % Store [a, b] in each row
end

% MATLAB function for the PV model (accepts a voltage vector)
modele_pv = @(V) piecewise_pv_model(V, points_V, coeffs);
function I_out = piecewise_pv_model(V_in, points_V, coeffs)
    I_out = zeros(size(V_in)); % Initialize the output current array
    for k = 1:numel(V_in) % Iterate over each voltage input
        V = V_in(k);
        found = false;
        for i = 1:length(points_V) - 1
            % Use a small tolerance for floating-point comparisons
            tol = 1e-9;
            if (points_V(i) - tol <= V && V <= points_V(i + 1) + tol)
                a = coeffs(i, 1);
                b = coeffs(i, 2);
                I_out(k) = a * V + b;
                found = true;
                break;
            end
        end
        if ~found
            I_out(k) = 0; % Default to 0 if V is outside defined segments
        end
    end
end

%% Provided Data: (R, V_test, I_test, V*, I*)
mesures = [
    
    % ---------------------CS6P-250P---------------------
    
    % --- Sref = 1000, Tref = 25, S = 556, T = 33 ---
    % 3.312101911,35.098381,10.830381,15.96633,4.926764;
    % 4.470212766,35.093315,8.052373,21.384974,4.906911;
    % 5.072649573,34.928596,7.081912,24.119793,4.890384;
    % 5.428571429,34.964813,6.637497,25.479982,4.836957;
    % 5.846827133,34.940865,6.137167,27.134899,4.766093;
    % 6.560465116,34.992397,5.522508,28.611193,4.515424;
    % 7.335802469,35.004902,4.953645,30.103748,4.260068;
    % 11.3297491,35.064438,3.307633,31.929382,3.011903;
    % 20.85,35.042049,1.906419,33.620697,1.829092;
    % 62.91512915,35.042969,0.812353,34.34491,0.796171;

    % --- Sref = 1000, Tref = 25, S = 765, T = 44.5 ---
    % 3.424196018,34.855976,10.283596,22.807545,6.728934;
    % 3.684782609,34.869011,9.658541,24.17536,6.696453
    % 4.031847134,34.885811,8.816074,25.775291,6.513734;
    % 4.44278607,34.931274,8.078631,27.240286,6.29992;
    % 5.227356747,34.942135,6.905494,28.715485,5.674942;
    % 5.995125914,34.961349,6.021228,29.939571,5.156351;
    % 7.701643489,34.964321,4.749564,30.825687,4.187371;
    % 13.45162653,34.977104,2.806052,32.122658,2.577053;
    % 68.56557377,35.020119,0.753771,33.615498,0.723538;

    % ----------Vertex N TSM-NEG21C.20 695W------------
    
    % --- Sref = 1000, Tref = 25, S = 200, T = 25 ---
    % 4.873925501,45.740959,9.644203,17.313847,3.650519;
    % 8.299711816,46.032471,5.740947,29.145626,3.634901;
    % 9.694117647,45.859604,4.960142,33.305958,3.602349;
    % 10.51461988,45.897758,4.461115,36.310131,3.529229;
    % 12.40837696,45.938385,3.932939,38.277752,3.277086;
    % 14.13879004,45.959175,3.487038,40.10532,3.042891;
    % 18.57272727,45.955261,2.704127,41.202496,2.424462;
    % 37.6460177,45.97047,1.46015,42.821026,1.360115;
    % 258.6470588,45.98262,0.451527,44.228886,0.434306;

    % ---------------------KC200GT---------------------
    
    % --- Sref = 1000, Tref = 25, S = 511, T = 54,3 ---
    2.969387755,30.078857,10.578335,11.909241,4.188322;
    3.847117794,29.85693,7.953677,15.661912,4.172224;
    4.721649485,30.06168,6.533874,18.635529,4.050412;
    5.596774194,29.935583,5.559875,21.14057,3.926395;
    6.62611276,29.938969,4.761699,22.638481,3.60058;
    7.350346566,29.951893,4.303252,23.602148,3.390971;
    9.27480916,29.965706,3.448748,24.549633,2.825413;
    12.8680203,29.976023,2.564936,25.56798,2.187756;
    17.76870748,29.986725,1.936755,26.319908,1.699926;
    60.97727273,29.984901,0.803749,26.9799,0.723199;

    % ---------------------KC85TS---------------------
    
    % --- Sref = 1000, Tref = 25, S = 1000, T = 25 ---
    % 2.560636183,23.837536,9.623577,13.199326,5.328769;
    % 2.747572816,23.819389,8.729925,14.48859,5.31014;
    % 2.943359375,23.843975,8.196189,15.374142,5.284747;
    % 3.305498982,23.867809,7.438157,16.523491,5.149376;
    % 3.544885177,23.967472,6.982988,17.287716,5.036823;
    % 4.417475728,23.891666,5.552080,18.478067,4.294037;
    % 5.591304348,23.895893,4.412038,19.506063,3.601518;
    % 7.574144487,23.901066,3.292860,20.110668,2.770655;
    % 12.26190476,23.910973,2.129044,20.725580,1.845415;
    % 31.58208955,23.911114,0.999190,21.232008,0.887237;
    % 405.8490566,23.910839,0.326262,21.545561,0.293988;

    % 2.252895753,23.937536,10.608313,12.032783,5.332525;
    % 2.688362919,23.787804,9.02202,14.035207,5.323145;
    % 2.854330709,23.817375,8.457655,14.91899,5.297799;
    % 3.147773279,23.797363,7.802836,15.958516,5.232584;
    % 3.579281184,23.820534,6.915626,17.328434,5.030826;
    % 4.507462687,23.915047,5.544626,18.489063,4.28663;
    % 5.81570997,23.86104,4.327792,19.562981,3.548232;
    % 7.995983936,23.897846,3.208596,20.156059,2.706212;
    % 22.87431694,23.8971677,1.333663,21.079178,1.176395;

    % -------------- KC85TS - charge R+L -------------- 
    
    % --- Sref = 1000, Tref = 25, S = 1000, T = 25 ---
    % 2.291750503,23.776310,10.092688,12.558353,5.330833;
    % 2.628458498,23.813004,8.758986,14.440410,5.311524;
    % 2.784046693,23.824802,8.113836,15.506546,5.280949;
    % 2.940828402,23.840160,7.710728,16.109459,5.210353;
    % 3.298780488,23.937588,6.910994,17.392550,5.021383;
    % 3.82735426,23.850681,5.991889,18.101677,4.547595;
    % 4.816489362,23.860895,4.968657,18.984102,3.953141;
    % 5.238095238,23.863398,4.379676,19.531202,3.584584;
    % 7.25,23.871922,3.343837,20.080969,2.812823;
    % 10.32307692,23.874191,2.385726,20.611444,2.059683;
    % 21.25510204,23.875790, 1.307175,21.090651, 1.154691;
    % 406.0377358,23.883148,0.306018,21.554970,0.276186;

    
    % ---------------Uni-Solar ES-62T-----------------
    
    % --- Sref = 1000, Tref = 25, S = 1000, T = 25 ---
    % 2.68852459,22.999146,8.928772,11.900667,4.620100;
    % 2.795348837,23.152687,8.527935,12.389590,4.563514;
    % 3.073286052,23.148125,7.678766,13.375956,4.437113;
    % 3.569620253,23.013630,6.702841,14.468855,4.214130;
    % 3.780927835,23.023594,6.254211,15.035897,4.084404;
    % 5.068535826,23.038422,4.691563,16.597288,3.379885;
    % 7.243902439,23.059406,3.418242,18.088655,2.681396;
    % 9.658031088,23.066004,2.557652,18.883083,2.093833;
    % 20.5625,23.081512,1.335646,19.927521,1.153136;
    % 92.95454545,23.119362,0.512176,20.578976,0.455897;
    % 92.95454545,23.119362,0.512176,20.578976,0.455897;
    
    ];
if ~isempty(mesures)
    R_data = mesures(:, 1);
    V_violet = mesures(:, 2); % Purple Points (Test Point)
    I_violet = mesures(:, 3); % Purple Points (Test Point)
    V_noir = mesures(:, 4);   % Black Points (Intersection)
    I_noir = mesures(:, 5);   % Black Points (Intersection)
else
    R_data = []; V_violet = []; I_violet = []; V_noir = []; I_noir = [];
end

%% Plotting (Figure Generation)

figure; 
hold on; 

% 1. Plot "Load Lines" FIRST

if isempty(R_data)
    Rs_to_plot = [0.1, 0.5, 1:1:15, 20:5:50, 60:10:200, 300, 400, 500, 1000];
    Rs_to_plot = sort([Rs_to_plot, inf]);
else
    Rs_to_plot = R_data; % Use resistances from your data
    Rs_to_plot = sort(unique(Rs_to_plot));
end
V_line = linspace(0, max(V_violet) * 1.25, 100);
for r_val = Rs_to_plot'
    if isinf(r_val) % Open circuit (I = 0)
        plot(V_line, zeros(size(V_line)), 'k--', 'LineWidth', 0.8, 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
    elseif r_val == 0 % Short circuit (V = 0)
        plot(zeros(size(V_line)), linspace(0, max(I_violet) * 1.1, 100), 'k--', 'LineWidth', 0.8, 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
    else
        I_line = V_line ./ r_val;
        plot(V_line, I_line, 'k--', 'LineWidth', 0.8, 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
    end
end

% 2. Plot the I-V Model (P-V Curve)

V_plot = linspace(0, Voc_mod_ref, 500);
I_plot = modele_pv(V_plot);
h_model = plot(V_plot, I_plot, 'r-', 'LineWidth', 2, 'DisplayName', 'I-V Curve - Single diode model');

% 3. Plot "Test Points" (Purple Points)

h_test_point = []; 
if ~isempty(V_violet)
    h_test_point = scatter(V_violet, I_violet, 40, 'm', 'filled', 'DisplayName', 'Test Point');
end

% 4. Plot "Intersections" (Black Points)

h_intersection = []; 
if ~isempty(V_noir) && ~isempty(I_noir)
    h_intersection = scatter(V_noir, I_noir, 70, 'ks', 'filled', 'DisplayName', 'Intersection');
end

%% Plot Configurations

xlabel('Voltage (V)', 'FontSize', 30); 
ylabel('Current (A)', 'FontSize', 30);

% Axis limits to match the provided figure
xlim([0 1.01*Voc_mod_ref]);
ylim([0 1.1 * max(I_violet)]);

% Font configurations for axis ticks
set(gca, 'FontSize', 30); 
grid off;
% Create the legend and adjust position and font
legend_handles = [h_model];
legend_labels = {'I-V Curve - Single diode model'};
if ~isempty(h_test_point) && ishandle(h_test_point)
    legend_handles = [legend_handles, h_test_point];
    legend_labels = [legend_labels, 'Test Point'];
end
if ~isempty(h_intersection) && ishandle(h_intersection)
    legend_handles = [legend_handles, h_intersection];
    legend_labels = [legend_labels, 'Intersection'];
end
leg = legend(legend_handles, legend_labels, 'Location', 'NorthWest', 'FontSize', 20);
set(leg, 'Box', 'on'); 
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');
set(findall(gca, '-property', 'FontName'), 'FontName', 'Times New Roman');
hold off;

%% --- Auxiliary Functions ---

function F = residuals_2_20(x, Voc, Isc, Vmp, Imp, q, k, Tref)
    Iph = x(1); Is0 = x(2); A = x(3); Rs = x(4); Rp = x(5);

    C   = q/(A*k*Tref);

    cap = 700;
    a_sc = min(C*Rs*Isc, cap);
    a_oc = min(C*Voc, cap);
    a_mp = min(C*(Rs*Imp+Vmp), cap);
    a_eq5 = min(C*Is0, cap); 

    F = zeros(5,1);

    F(1) = Isc - ( Iph - Is0*(exp(a_sc)-1) - (Rs*Isc)/Rp );
    F(2) = 0   - ( Iph - Is0*(exp(a_oc)-1) - (Voc)/Rp );
    F(3) = Imp - ( Iph - Is0*(exp(a_mp)-1) - (Rs*Imp+Vmp)/Rp );
    F(4) = Rs + (q*Is0*Rp*(Rs-Rp)/(A*k*Tref))*exp(a_eq5);
    term_coeff = 1 + q*(Vmp - Rs*Imp)/(A*k*Tref);
    F(5) = Iph - 2*Vmp/Rp + Is0 - Is0 * term_coeff * exp(a_mp);
end

function I = solve_I_V_2_11(V, Iph_ref, Is0_ref, A, Rs, Rp, ...
    q, k, S, Sref, alpha, T, Tref, E_G0, k1, k2) 
    
    % 1. Photocurrent (Iph) 
    Iph = Iph_ref * (S / Sref) * (1 + alpha * (T - Tref));
    
    % 2. Calculation of Saturation Current (Is)
    
    % 2a. Calculation of band gap energy in eV at temperature T: Eg(T)
    Eg_T_eV = E_G0 - (k1 * T^2) / (T + k2);
    
    % 2b. Calculation of the temperature difference term
    temp_diff = (1 / Tref) - (1 / T);
    
    % 2c. Calculation of the exponent 
    exponent_term = (q / (A * k)) * Eg_T_eV * temp_diff; 
    
    % 2d. Final calculation of Is(T)
    Is = Is0_ref * (T / Tref)^3 * exp(exponent_term);
    
    % 3. Iterative Solver (Newton-Raphson)
    I = Iph; % Initial guess
    for it = 1:30
        V_diode = V + I * Rs;
        arg_exp = min(q * V_diode / (A * k * T), 700); 
        expo = exp(arg_exp);
        
        f = Iph - Is * (expo - 1) - V_diode / Rp - I;
        df = -Is * expo * (q * Rs / (A * k * T)) - Rs / Rp - 1;
        
        dI = -f / df;
        I = I + dI;
        
        if abs(dI) < 1e-6, break;end
    end
end
