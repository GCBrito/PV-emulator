clc; clear; close all; format long
%% Module Parameters

ns      = 60;       % number of series cells
np      = 1;        % number of parallel branches

Vmp_mod_ref  = 27.5;     % voltage at maximum power point (V)
Imp_mod_ref  = 6.60;     % current at maximum power point (A)
Voc_mod_ref  = 34.2;     % open-circuit voltage (V)
Isc_mod_ref  = 7.19;     % short-circuit current (A)

Tref = 25 + 273.15; % Reference temperature (K)

%% Operating Conditions

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

