clc; clear; close all; format long

%% Module Parameters

ns = 72;       % number of series cells

Vmp_mod_ref = 17.4;     % voltage at maximum power point (V)
Imp_mod_ref = 5.02;     % current at maximum power point (A)
Voc_mod_ref = 21.7;     % open-circuit voltage (V)
Isc_mod_ref = 5.34;     % short-circuit current (A)

Tref = 25 + 273.15; % Reference temperature (K)
Sref = 1000; % Reference irradiance (W/m²)

alpha   = 0.00040;  % temperature coefficient of Isc (%/K)
beta = -0.0038; % temperature coefficient of Voc (%/K)

%% Operating Conditions

T = 25 + 273.15; % Current temperature (K)
S = 1000; % Current irradiance (W/m²)

%% Physical Constants

q = 1.60217662e-19; % Elementary charge (C)
k = 1.38064852e-23; % Boltzmann constant (J/K)
E_G0 = 1.166;            % Band gap energy at 0K (eV)
k1   = 4.73e-4;          % Coefficient k1 (eV/K)
k2   = 636;              % Coefficient k2 (K)

%% Parameter Estimation via fsolve

Rs_0 = (Voc_mod_ref - Vmp_mod_ref)/Imp_mod_ref;
Rp_0 = (Vmp_mod_ref)/(Isc_mod_ref - Imp_mod_ref);

x0 = [Isc_mod_ref; log10(1e-9); ns; log10(Rs_0); Rp_0];
opts = optimoptions('fsolve', ...
    'Display', 'off', 'TolFun', 1e-10, 'TolX', 1e-10, 'MaxIter', 1000, 'MaxFunctionEvaluations', 2000);
fun = @(x) residuals_2_20(x, Voc_mod_ref, Isc_mod_ref, Vmp_mod_ref, Imp_mod_ref, q, k, Tref);
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

%% Characteristic Points (Adaptive Embedded Strategy - High Precision Knee)

Voc_estimation = Voc_mod_ref * (1 + beta * (T - Tref));

Vmp_estimation = Vmp_mod_ref * (Voc_estimation / Voc_mod_ref);

factors = [
    0.00;  
    0.20;   
    0.40;   
    0.60;   
    0.70;   
    0.80;   
    0.85; 
    0.88; 
    0.90; 
    0.92; 
    0.93; 
    0.94; 
    0.95; 
    0.96; 
    0.97; 
    0.98; 
    0.99; 
    1.00;   
    1.01; 
    1.02; 
    1.03; 
    1.04;
    1.05;   
    1.08;
    1.12;   
];
V_mod = factors * Vmp_estimation;

V_mod = V_mod(V_mod < Voc_estimation); 
V_mod = V_mod(V_mod >= 0);

V_mod = [V_mod; Voc_estimation];

V_mod = unique(V_mod);

I_mod = arrayfun(@(V) solve_I_V_2_11(V, Iph_ref, Is0_ref, A, Rs, Rp, ...
    q, k, S, Sref, alpha, T, Tref, E_G0, k1, k2, ns), V_mod);

points_V = V_mod;
points_I = I_mod;

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
    
    % --- Sref = 1000, Tref = 25, S = 765, T = 44.5 ---
    % 2.890855457	35.784859	12.165794	20.089266	6.829757;
    % 3.030390738	35.786652	11.396839	21.406616	6.817284;
    % 3.191642651	35.784744	10.764368	22.612379	6.802004;
    % 3.379360465	35.78397	10.236521	23.688374	6.776401;
    % 3.54862119	35.786598	9.671537	24.878904	6.723669;
    % 3.739002933	35.785908	9.176556	25.916357	6.645714;
    % 3.90922619	35.783024	8.804281	26.656096	6.558634;
    % 4.1	35.787361	8.393621	27.416744	6.430365;
    % 4.330218069	35.783806	7.966219	28.124405	6.261077;
    % 4.560386473	35.786224	7.636677	28.616514	6.106682;
    % 4.903878583	35.786362	7.083742	29.345022	5.808709;
    % 5.399638336	35.785271	6.442036	30.058086	5.411033;
    % 6.138832998	36.241776	5.862995	30.656363	4.959418;
    % 6.928888889	35.78756	5.102286	31.281769	4.459889;
    % 8.252604167	35.788845	4.317395	31.720699	3.826633;
    % 13.49152542	35.790173	2.6143	32.716629	2.389792;
    % 431.375	35.815018	0.287514	34.182873	0.274412;

    % ---------------- KC85TS - charge R ----------------- 
    
    % --- Sref = 1000, Tref = 25, S = 1000, T = 25 ---
    2.274319066	22.960472	10.105235	12.117974	5.333295;


    % --------------- KC85TS - charge R+L --------------- 
    
    % ---- Sref = 1000, Tref = 25, S = 1000, T = 25 ----
    % 2.111553785	23.090778	11.210976	10.992171	5.336891	2.059658142;
    % 2.520710059	22.955368	9.275552	13.190496	5.32987	2.47482509;
    % 2.990176817	22.999069	7.755183	15.65108	5.277474	2.965638485;
    % 3.31352459	22.993301	7.134839	16.572636	5.142501	3.222680171;
    % 3.611464968	23.009609	6.611756	17.421389	5.005994	3.480105849;
    % 3.87032967	23.01091	5.902045	18.000225	4.61686	3.898802433;
    % 4.431372549	23.01825	5.383958	18.448944	4.3152	4.275339266;
    % 5.054200542	23.015213	4.812998	18.967569	3.966545	4.781886755
    % 5.615835777	23.016911	4.287896	19.471916	3.627488	5.367878819;
    % 6.74137931	23.026741	3.620827	19.864429	3.123571	6.359525364;
    % 8.453389831	23.02433	2.97291	20.225323	2.611501	7.744711949;
    % 11.05978261	23.027668	2.328817	20.59774	2.083076	9.888136583;
    % 38.29090909	23.027573	0.983892	21.222818	0.906781	23.40456847;

    % ---------------------KC200GT---------------------
    
    % --- Sref = 1000, Tref = 25, S = 511, T = 54,3 ---
    % 3.206030151	29.875959	9.626084	13.112031	4.224718;
    % 3.914141414	29.897739	7.932604	15.875197	4.21208;
    % 4.479899497	30.066505	6.832339	18.182568	4.131823;
    % 5.492021277	29.942114	5.744508	21.021599	4.033073;
    % 6.223163842	30.058102	5.098278	22.375362	3.795177;
    % 7.027108434	29.97991	4.516749	23.667114	3.565669;
    % 8.328719723	29.98394	3.9096	24.400532	3.181581;
    % 10.57142857	29.996906	3.105203	25.443676	2.633864;
    % 13.86702128	29.997562	2.43492	26.335112	2.137636;
    % 26.11925709	30.003111	1.506433	26.927042	1.351986;
    % 73.54054054	30.003746	0.800091	27.395672	0.733724;

    % --------------------- ME Solar MESM-50W ---------------------
    
    % --- Sref = 1000, Tref = 25, S = 1000, T = 25 ---
    % 2.564575646	23.063789	9.983421	6.993855	3.027369;
    % 3.848920863	22.905588	6.327226	10.947917	3.02415
    % 5.124087591	22.933201	4.818062	14.322508	3.009032;
    % 5.654676259	22.961166	4.267603	15.97874	2.969837;
    % 6.296992481	22.959677	3.889386	17.005857	2.880805;
    % 6.791505792	22.976776	3.596284	17.847771	2.793502;
    % 7.541666667	22.976551	3.281668	18.351051	2.621022;
    % 9.302439024	22.984993	2.730845	19.306019	2.293747;
    % 11.36571429	22.991455	2.275814	20.117302	1.991315;
    % 16.83606557	23.000027	1.686832	20.741932	1.521222;
    % 26.70886076	23.002243	1.225857	21.237494	1.131808;
    % 135.1875	23.006454	0.591647	21.773888	0.559949;

    % ---------------------Renogy RNG-50DB-H – 50 W---------------------

    % --- Sref = 1000, Tref = 25, S = 1000, T = 25 ---
    % 2.326860841	24.066233	9.584604	7.320326	2.915389;
    % 3.95970696	24.055067	6.341636	11.047176	2.912366;
    % 5.248175182	24.049904	4.764365	14.633799	2.899004;
    % 6.222641509	23.951168	4.069241	16.772703	2.849638;
    % 7.268	23.976217	3.52643	18.451015	2.713782;
    % 8.852534562	23.995972	2.966523	19.458088	2.405523;
    % 11	24.001507	2.436483	20.484665	2.079475;
    % 15.10144928	24.003086	1.895893	21.060669	1.663485;
    % 21.0990099	24.005049	1.487081	21.505325	1.332226;
    % 75.55172414	24.004025	0.724881	22.052029	0.665934;
   
    % ---------------Uni-Solar ES-62T-----------------
    
    % --- Sref = 1000, Tref = 25, S = 1000, T = 25 ---
    % 1.944812362	22.92971	12.203246	9.063915	4.823837;
    % 2.292410714	23.088787	10.189253	10.677806	4.712195;
    % 2.689497717	22.956436	8.67457	12.169374	4.598453;
    % 3.227941176	23.073864	7.521092	13.542478	4.414268;
    % 3.970588235	22.967014	6.090907	15.171309	4.023467;
    % 4.56851312	22.983368	5.228667	16.019594	3.644424;
    % 5.514950166	22.987225	4.418077	16.904531	3.249001;
    % 6.762452107	22.989025	3.577474	17.930962	2.790355;
    % 10.0972973	22.998974	2.547448	18.904875	2.093971;
    % 27.16438356	23.000847	1.247732	20.006903	1.085319;

    % ----------Vertex N TSM-NEG21C.20 695W------------
    
    % --- Sref = 1000, Tref = 25, S = 200, T = 25 ---
    % 4.45505618	46.11351	10.403705	16.184828	3.651472;
    % 7.296829971	45.855804	6.508248	25.680521	3.644799;
    % 9.152974504	45.915001	5.082425	32.682323	3.617673;
    % 10.08235294	45.960556	4.736919	34.661255	3.572358;
    % 10.65671642	45.950321	4.505582	36.09742	3.539472;
    % 11.67823344	45.95752	4.166992	37.399529	3.391034;
    % 13.18983051	46.005436	3.688589	39.285419	3.149796;
    % 14.17437722	46.008614	3.468349	40.211502	3.031335;
    % 18.66666667	46.008034	2.680178	41.241325	2.402495;
    % 39.89719626	46.058121	1.369122	42.957077	1.276941;
    % 199.9090909	46.031281	0.478006	44.203873	0.459029;
    
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
V_plot = linspace(0, max(points_V), 500); 
I_plot = modele_pv(V_plot);
h_model = plot(V_plot, I_plot, 'r-', 'LineWidth', 2, 'DisplayName', 'I-V model');

% 3. Plot "Test Points" (Purple Points)
h_test_point = []; 
if ~isempty(V_violet)
    h_test_point = scatter(V_violet, I_violet, 40, 'm', 'filled', 'DisplayName', 'Test point');
end

% 4. Plot "Intersections" (Black Points)
h_intersection = []; 
if ~isempty(V_noir) && ~isempty(I_noir)
    h_intersection = scatter(V_noir, I_noir, 70, 'ks', 'filled', 'DisplayName', 'Intersection');
end

%% Plot Configurations

xlabel('Voltage (V)', 'FontSize', 30); 
ylabel('Current(A)', 'FontSize', 30);

% Axis limits to match the provided figure
xlim([0 1.2*Voc_estimation]);
ylim([0 1.1 * max(I_violet)]);

% Font configurations for axis ticks
set(gca, 'FontSize', 30); 
grid off;

% Create the legend and adjust position and font
legend_handles = [h_model];
legend_labels = {'I-V model'};
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
    a_eq5 = min(C*Isc, cap); 
    F = zeros(5,1);
    F(1) = Isc - ( Iph - Is0*(exp(a_sc)-1) - (Rs*Isc)/Rp );
    F(2) = 0   - ( Iph - Is0*(exp(a_oc)-1) - (Voc)/Rp );
    F(3) = Imp - ( Iph - Is0*(exp(a_mp)-1) - (Rs*Imp+Vmp)/Rp );
    F(4) = Rs + (q*Is0*Rp*(Rs-Rp)/(A*k*Tref))*exp(a_eq5);
    term_coeff = 1 + q*(Vmp - Rs*Imp)/(A*k*Tref);
    F(5) = Iph - 2*Vmp/Rp + Is0 - Is0 * term_coeff * exp(a_mp);
end

function I = solve_I_V_2_11(V, Iph_ref, Is0_ref, A, Rs, Rp, ...
    q, k, S, Sref, alpha, T, Tref, E_G0, k1, k2, ns) 
    
    % 1. Photocurrent (Iph) 
    Iph = Iph_ref * (S / Sref) * (1 + alpha * (T - Tref));
    
    % 2. Calculation of Saturation Current (Is)
    
    % 2a. Calculation of band gap energy in eV at temperature T: Eg(T)
    Eg_T_eV = E_G0 - (k1 * T^2) / (T + k2);
    
    % 2b. Calculation of the temperature difference term
    temp_diff = (1 / Tref) - (1 / T);
    
    % 2c. Calculation of the exponent 
    exponent_term = (ns*q / (A * k)) * Eg_T_eV * temp_diff; 
    
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
