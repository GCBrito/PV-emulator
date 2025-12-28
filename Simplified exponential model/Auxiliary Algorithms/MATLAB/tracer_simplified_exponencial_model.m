clc; clear; close all;
%% 1) Panel Parameters

Vmp = 17.4; % Maximum Power Voltage (V)
Imp = 5.02;  % Maximum Power Current (A)
Voc = 21.7; % Open-Circuit Voltage (V)
Isc = 5.34;  % Short-Circuit Current (A)

%% 2) Simplified exponential equation constant

% This constant 'c' is derived from a simplified exponential model,

c = - (Voc - Vmp) / log(1 - Imp / Isc);

%% 3) Characteristic Points for the linear interpolation

points_V = [
    0;
    0.2 * Vmp;
    0.4 * Vmp;
    0.6 * Vmp;
    0.8 * Vmp;       
    0.9 * Vmp;
    0.95 * Vmp;
    Vmp;
    (Vmp + Voc) / 2;          
    0.9 * Voc;
    0.95 * Voc;
    Voc
];
points_I = Isc * (1 - exp((points_V - Voc) / c)); % Current derived from the simplified exp model

%% 4) Linear Interpolation

% Initialization of coefficients for the linear interpolation
coeffs = zeros(length(points_V) - 1, 2);
for i = 1:length(points_V) - 1
    a = (points_I(i + 1) - points_I(i)) / (points_V(i + 1) - points_V(i));
    b = points_I(i) - a * points_V(i);
    coeffs(i,:) = [a, b]; % Store [a, b] in each row
end
pv_model = @(V_in) piecewise_pv_model(V_in, points_V, coeffs);
function I_out = piecewise_pv_model(V_in, points_V, coeffs)
    I_out = zeros(size(V_in)); % Initialize the output current array
    for k = 1:numel(V_in) % Iterate over each voltage input
        V = V_in(k);
        found = false;
        for i = 1:length(points_V) - 1
            % Use a small tolerance for floating-point comparisons
            tol = 1e-9;
            if (V >= (points_V(i) - tol) && V <= (points_V(i + 1) + tol))
                a = coeffs(i, 1);
                b = coeffs(i, 2);
                I_out(k) = a * V + b;
                found = true;
                break;
            end
        end
        if ~found
            % If V is outside all segments, typically 0 for V > Voc and V < 0
            if V < points_V(1) || V > points_V(end)
                 I_out(k) = 0; % Or handle as needed
            else
                 % This case should ideally not be reached if segments cover the range
                 I_out(k) = NaN; 
            end
        end
    end
end

%% 5) Provided Data: (R, V_test, I_test, V*, I*)
measurements = [
    
     % ---------------------KC200GT---------------------
     % 3.699029126,30.056282,7.876191,15.578393,4.082288;
     % 4.289340102,30.050781,7.080119,17.230631,4.059626;
     % 4.693298969,29.872391,6.489758,18.552963,4.030620;
     % 4.913265306,29.879427,6.065750,19.599977,3.978944;
     % 5.397849462,29.889881,5.765186,20.419359,3.938504;
     % 5.911111111,29.917889,5.258528,21.612139,3.798665;
     % 6.507204611,29.941908,4.745583,22.879311,3.626210;
     % 7.993311037,29.952457,3.937661,24.182100,3.179068;
     % 9.636015326,29.963566,3.255278,25.403482,2.759865;
     % 12.20283019,29.978542,2.608978,26.099630,2.271403;
     % 15.66272189,29.987871,2.071892,26.694136,1.844324;
     % 40.19117647,29.994781,0.951167,27.493385,0.871845;
     % 410.2941176,30.012428,0.284112,27.971466, 0.264791;

     % ---------------------KC85TS---------------------
     2.528599606,23.775753,9.540161,13.240793,5.312946;
     2.694552529,23.927464,8.889205,14.250062,5.293989;
     2.948207171,23.923685,8.276764,15.192984,5.256245;
     3.116,23.798820,7.758160,15.953911,5.200803;
     3.481404959,23.811413,6.968660,17.226923,5.041640;
     4.040540541,23.839123,6.009305,18.291187,4.610795;
     4.905370844,23.851524,4.968151,19.491007,4.059877;
     6.158385093,23.908186,3.995351,20.092628,3.357725;
     8.023622047,23.868225,3.134047,20.618151,2.707292;
     15.85606061,23.903660,1.681039,21.106844,1.484351;
     405.6603774,23.868862,0.325685,21.582323,0.294486;

     % ---------------Uni-Solar ES-62T-----------------
     % 2.159645233,22.119232,10.492447,10.140304,4.810140;
     % 2.591224018,21.972612,8.881744,11.602366,4.689895;
     % 2.891203704,22.107285,7.773019,12.875836,4.527201;
     % 3.356968215,21.984802,6.669855,14.143713,4.290987;
     % 3.854497354,22.010578,6.045714,14.959622,4.109006;
     % 5.061919505,22.015034,4.473690,16.689295,3.391443;
     % 6.518518519,22.034122,3.533989,17.921898,2.874441;
     % 9.083373964,22.046028,2.588672,18.896595,2.218861;
     % 18.08256881,22.052059,1.446311,19.905468,1.305524;
     % 58.02857143,22.056940,0.659405,20.491583,0.612608;
];

% Extract columns for easier use
if ~isempty(measurements)
    R_data = measurements(:, 1);
    V_magenta = measurements(:, 2); % V_test -> Magenta Points (Test Point)
    I_magenta = measurements(:, 3); % I_test -> Magenta Points (Test Point)
    V_black = measurements(:, 4);   % V* -> Black Points (Intersection)
    I_black = measurements(:, 5);   % I* -> Black Points (Intersection)
else
    R_data = []; V_magenta = []; I_magenta = []; V_black = []; I_black = [];
end

%% 6) Plotting (Figure Generation)

figure; 
hold on; %

max_V_plot = max([Voc, max(V_magenta), max(V_black)]);
min_I_plot = 0; 
max_I_plot = max([Isc, max(I_magenta), max(I_black)]) * 1.1; 

% 1. Plot "Load Lines" 
Rs_to_plot = unique(R_data); 
if isempty(Rs_to_plot)
    Rs_to_plot = [0.1, 0.5, 1:1:15, 20:5:50, 60:10:200, 300, 400, 500, 1000];
    Rs_to_plot = sort(Rs_to_plot);
end
V_line_range = linspace(0, max_V_plot * 1.1, 100); 
for r_val = Rs_to_plot'
    if isinf(r_val) % Open Circuit (I = 0)
        plot(V_line_range, zeros(size(V_line_range)), '--', 'LineWidth', 0.8, 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
    elseif r_val == 0 % Short Circuit (V = 0)
        plot(zeros(size(V_line_range)), linspace(0, max_I_plot, 100), '--', 'LineWidth', 0.8, 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
    else
        I_line = V_line_range ./ r_val;
        plot(V_line_range, I_line, '--', 'LineWidth', 0.8, 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
    end
end
h_load_line_legend = plot(NaN, NaN, '--', 'Color', [0.5 0.5 0.5], 'DisplayName', 'Load line');

% 2. Plot the I-V Model (PV Curve)
V_plot = linspace(0, Voc, 500);
I_plot = pv_model(V_plot);
h_model = plot(V_plot, I_plot, 'r-', 'LineWidth', 2, 'DisplayName', 'I-V model');

% 3. Plot "Test Points" (Magenta Points)s
h_test_point = []; 
if ~isempty(V_magenta)
    h_test_point = scatter(V_magenta, I_magenta, 40, 'm', 'filled', 'DisplayName', 'Test point');
end

% 4. Plot "Intersections" (Black Points)
h_intersection = []; 
if ~isempty(V_black) && ~isempty(I_black)
    h_intersection = scatter(V_black, I_black, 70, 'ks', 'filled', 'DisplayName', 'Intersection');
end

%% Plot Configurations

xlabel('Voltage (V)', 'FontSize', 30); 
ylabel('Current (A)', 'FontSize', 30); 

% Axis limits
xlim([0 1.2 * Voc]); 
ylim([0 1.2 * max(I_magenta)]); 
set(gca, 'FontSize', 30); 
grid off;

legend_handles = [h_model];
legend_labels = {'I-V model'};

if exist('h_load_line_legend', 'var') && ishandle(h_load_line_legend)
    legend_handles = [legend_handles, h_load_line_legend];
    legend_labels = [legend_labels, 'Load line'];
end
if ~isempty(h_test_point) && ishandle(h_test_point)
    legend_handles = [legend_handles, h_test_point];
    legend_labels = [legend_labels, 'Test point'];
end
if ~isempty(h_intersection) && ishandle(h_intersection)
    legend_handles = [legend_handles, h_intersection];
    legend_labels = [legend_labels, 'Intersection'];
end
leg = legend(legend_handles, legend_labels, 'Location', 'NorthWest', 'FontSize', 20); 

set(leg, 'Box', 'on'); 
set(findall(gcf, '-property','FontName'),'FontName','Times New Roman');
set(findall(gca, '-property','FontName'),'FontName','Times New Roman');
hold off; 