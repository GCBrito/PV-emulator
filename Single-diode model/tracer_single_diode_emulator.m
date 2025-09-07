clc; clear; close all; format long

%% 1) Module Parameters

ns = 60; % Number of series cells
np = 1; % Number of parallel branches

VMPmod = 30.1; % Voltage at maximum power point of the module (V)
IMPmod = 8.30; % Current at maximum power point of the module (A)
VOCmod = 37.2; % Open-circuit voltage of the module (V)
ISCmod = 8.87; % Short-circuit current of the module (A)

Vmp_cell = VMPmod / ns;
Imp_cell = IMPmod / np;
Voc_cell = VOCmod / ns;
Isc_cell = ISCmod / np;

Tref = 25 + 273.15; % Reference temperature (K)
Sref = 1000; % Reference irradiance (W/m²)

muICC = 0.0058; % Temperature coefficient
%% 2) Operating Conditions

T = 44.5 + 273.15; % Current temperature (K)
S = 765; % Current irradiance (W/m²)

%% 3) Physical Constants

q = 1.60217662e-19; % Elementary charge (C)
k = 1.38064852e-23; % Boltzmann constant (J/K)
Eg = 1.21 * q; % Silicon band gap energy (J)

%% 4) Parameter Estimation via fsolve

x0 = [Isc_cell; 1e-9; 3; 0.01; 2];
opts = optimoptions('fsolve', ...
    'Display', 'off', 'TolFun', 1e-10, 'TolX', 1e-10, 'MaxIter', 1000, 'MaxFunctionEvaluations', 2000);
fun = @(x) residuals_2_20(x, Voc_cell, Isc_cell, Vmp_cell, Imp_cell, q, k, Tref);
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
    0.2 * VMPmod / ns;
    0.4 * VMPmod / ns;
    0.6 * VMPmod / ns;
    0.8 * VMPmod / ns;
    0.9 * VMPmod / ns;
    VMPmod / ns;
    (VMPmod / ns + VOCmod / ns) / 2.0;
    0.9 * VOCmod / ns;
    0.95 * VOCmod / ns;
    VOCmod / ns];

I_cell = arrayfun(@(V) solve_I_V_2_11(V, Iph_ref, Is0_ref, A, Rs, Rp, ...
    q, k, Eg, S, Sref, muICC, T, Tref), V_cell);
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
    % ---------------------KC200GT---------------------
    
    % --- Sref = 1000, Tref = 25, S = 511, T = 54,3 ---
    % 3.409547739,30.042597,9.221842,13.882890,4.261477;
    % 3.837837838,29.833204,7.930077,15.970324,4.245132;
    % 4.282878412,30.036655,7.134727,17.581207,4.176135;
    % 4.674074074,29.913689,6.374939,19.258860,4.104277;
    % 5.430446194,29.912422,5.737882,21.006065,4.029441;
    % 6.206703911,29.915359,4.927178,22.505482,3.706742;
    % 7.051359517,30.027420,4.394181,23.636152,3.458890;
    % 8.111864407,29.946547,3.833536,24.225838,3.101213;
    % 11.20089286,29.959095,2.824838,25.353563,2.390584;
    % 15.69277108,29.962294,2.051082,26.291145,1.799772;
    % 34.67532468,29.970604,1.060586,26.867710,0.950782;
    % 405.9701493,29.976196,0.295321,27.326483,0.269217;

    % --- Sref = 511, Tref = 54,3, S = 511, T = 54,3 ---
    % 3.298200514,30.083759,9.238234,13.129311,4.031798;
    % 4.58778626,29.899143,6.469102,18.305553,3.960665;
    % 5.436997319,29.930714,5.628777,20.560305,3.866576;
    % 5.991689751,29.945515,5.095525,21.877769,3.722718;
    % 6.424068768,29.967121,4.802621,22.681322,3.634977;
    % 7.182098765,29.964804,4.303484,23.515127,3.377195;
    % 8.49825784,29.994026,3.657519,24.605215,3.000399;
    % 9.86328125,29.990040,3.184377,25.463463,2.703740;
    % 12.61927945,30.060379,2.518887,26.104248,2.187386;
    % 16.58125,30.009392,1.938024,26.694448,1.723943;
    % 31.30484988,30.022461, 1.095335,27.238712,0.993773;
    % 65.39379475,30.024586,0.648094,27.517773,0.593983;
    % 40.69117647,30.027824,0.280485,27.751440,0.259222;

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

    % ---------------------CS6P-250P---------------------
    
    % --- Sref = 1000, Tref = 25, S = 556, T = 33 ---
    % 4.966173362,36.916618,7.636598,23.861177,4.935940;
    % 5.235789474,36.925896,7.135579,25.287460,4.886562;
    % 5.670940171,36.925285,6.611194,26.927507,4.821167;
    % 6.189309577,36.964649,6.070982,28.176317,4.627609;
    % 7.035799523,37.076538,5.380336,29.874315,4.335190;
    % 8.924637681,36.992016,4.264987,31.124376,3.588478;
    % 12.23137255,37.003716,3.019732,32.481518,2.650693;
    % 18.71348315,37.020912,2.094456,33.569107,1.899170;
    % 49.46376812,37.018486,0.900303,34.318546,0.834639;
    % 406.4705882,37.029892,0.242107,34.741459,0.227145;

    % --- Sref = 1000, Tref = 25, S = 765, T = 44.5 ---
    3.276853253,36.058861,11.093060,22.148668,6.813763;
    3.473604827,35.893806,10.355574,23.533447,6.789537;
    3.797527048,35.901566,9.535309,25.083090,6.661966;
    4.232854864,36.052902,8.577839,27.035177,6.432309;
    4.633445946,35.951599, 7.810769,27.919050,6.065634;
    5.251838235,35.942787,6.917910,29.022854,5.586030;
    6.098138748,36.042290,5.907085,30.229673,4.954437;
    13.33891213,35.980782, 2.834330,32.149712,2.532543;
    61.90740741,35.987362,0.760594,33.594578,0.710022;
    421.125,35.992996,0.285993,33.809090,0.268640;

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
V_plot = linspace(0, VOCmod, 500);
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
xlim([0 1.2 * VOCmod]);
ylim([0 1.2 * max(I_violet)]);

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
    C = q / (A * k * Tref);
    cap = 700;
    a_sc = min(C * Rs * Isc, cap);
    a_oc = min(C * Voc, cap);
    a_mp = min(C * (Rs * Imp + Vmp), cap);
    F = zeros(5, 1);
    F(1) = Isc - (Iph - Is0 * (exp(a_sc) - 1) - (Rs * Isc) / Rp);
    F(2) = 0 - (Iph - Is0 * (exp(a_oc) - 1) - (Voc) / Rp);
    F(3) = Imp - (Iph - Is0 * (exp(a_mp) - 1) - (Rs * Imp + Vmp) / Rp);
    F(4) = Rs + (q * Is0 * Rp * (Rs - Rp) / (A * k * Tref)) * exp(a_sc);
    term = 1 + q * (Vmp - Rs * Imp) / (A * k * Tref);
    F(5) = Iph - 2 * Vmp / Rp - (Is0 * term * (exp(a_mp) - 1));
end

function I = solve_I_V_2_11(V, Iph_ref, Is0_ref, A, Rs, Rp, ...
    q, k, Eg, S, Sref, muICC, T, Tref)
    Iph = Iph_ref * (S / Sref) + muICC * (T - Tref);
    gamma = (T / Tref)^3 * exp(Eg / (A * k) * (1 / Tref - 1 / T));
    I = Is0_ref * gamma; % Initial estimate
    cap = 700;
    for it = 1:30
        C = q / (A * k * T);
        a = min(C * (Rs * I + V), cap);
        expo = exp(a);
        f = Iph ...
            - (Rs * I + V) / Rp ...
            - Is0_ref * gamma * (expo - 1) ...
            - I;
        df = -Rs / Rp ...
            - Is0_ref * gamma * expo * (q * Rs / (A * k * T)) ...
            - 1;
        dI = -f / df;
        I = I + dI;
        if abs(dI) < 1e-6, break; end
    end
end