%% Newton-Raphson Power Flow Solution - 3 Bus System
% Three-bus power system with generators at buses 1 and 3
% Bus 1: Slack bus, V1 = 1.05 pu, angle = 0
% Bus 2: Load bus, P = 4.0 pu, Q = 2.5 pu (400 MW, 250 Mvar on 100 MVA base)
% Bus 3: Voltage-regulated bus, |V3| = 1.04 pu, P_gen = 2.0 pu (200 MW)
% Line impedances in pu on 100 MVA base

clc; clear; close all;

fprintf('============================================================\n');
fprintf('   Newton-Raphson Power Flow Solution - 3 Bus System\n');
fprintf('============================================================\n\n');

%% System Data
baseMVA = 100;

% Line impedances (per unit)
% Line 1-2: z12 = 0.02 + j0.04
% Line 1-3: z13 = 0.01 + j0.03
% Line 2-3: z23 = 0.0125 + j0.025
z12 = 0.02 + 1j*0.04;
z13 = 0.01 + 1j*0.03;
z23 = 0.0125 + 1j*0.025;

%% Build Y-Bus Matrix
y12 = 1/z12;
y13 = 1/z13;
y23 = 1/z23;

Ybus = zeros(3,3);
Ybus(1,1) = y12 + y13;
Ybus(2,2) = y12 + y23;
Ybus(3,3) = y13 + y23;
Ybus(1,2) = -y12; Ybus(2,1) = -y12;
Ybus(1,3) = -y13; Ybus(3,1) = -y13;
Ybus(2,3) = -y23; Ybus(3,2) = -y23;

G = real(Ybus);
B = imag(Ybus);

fprintf('Y-Bus Matrix (Real part G):\n');
disp(G);
fprintf('Y-Bus Matrix (Imaginary part B):\n');
disp(B);

%% Bus Data Specification
% Bus types: 1=Slack, 2=Load (PQ), 3=Voltage-regulated (PV)
% Scheduled power in per unit (on 100 MVA base)

% P_scheduled (generation - load) per bus
P_sch = [0; -4.0; 2.0];   % Bus 1: slack; Bus 2: -4 pu load; Bus 3: +2 pu gen
Q_sch = [0; -2.5; 0  ];   % Bus 2: -2.5 pu load; Bus 3: Q free (PV bus)

% Voltage specifications
V    = [1.05; 1.0;  1.04];  % Initial voltages (|V|)
del  = [0;    0;    0   ];  % Initial phase angles (radians)

% Bus 3 voltage magnitude is fixed
V3_fixed = 1.04;

%% Newton-Raphson Iteration
tolerance = 1e-6;
max_iter  = 50;
iter      = 0;
converged = false;

fprintf('Starting Newton-Raphson Iterations...\n\n');
fprintf('%-6s %-16s %-16s %-16s %-16s\n', 'Iter', '|dP2|', '|dQ2|', '|dP3|', '|V2|');
fprintf('%s\n', repmat('-',72,1));

while iter < max_iter
    iter = iter + 1;

    %% Calculate P and Q injections at each bus
    P_calc = zeros(3,1);
    Q_calc = zeros(3,1);
    for i = 1:3
        for k = 1:3
            P_calc(i) = P_calc(i) + V(i)*V(k)*(G(i,k)*cos(del(i)-del(k)) + B(i,k)*sin(del(i)-del(k)));
            Q_calc(i) = Q_calc(i) + V(i)*V(k)*(G(i,k)*sin(del(i)-del(k)) - B(i,k)*cos(del(i)-del(k)));
        end
    end

    %% Mismatch vectors
    % For PQ bus (bus 2): dP2, dQ2
    % For PV bus (bus 3): dP3 only
    dP2 = P_sch(2) - P_calc(2);
    dQ2 = Q_sch(2) - Q_calc(2);
    dP3 = P_sch(3) - P_calc(3);

    mismatch = [dP2; dP3; dQ2];   % Order: [dP2, dP3, dQ2]

    fprintf('%-6d %-16.8f %-16.8f %-16.8f %-16.8f\n', iter, abs(dP2), abs(dQ2), abs(dP3), V(2));

    if max(abs(mismatch)) < tolerance
        converged = true;
        break;
    end

    %% Build Jacobian Matrix [J1 J2; J3 J4]
    % State variables: [del2, del3, |V2|]  (del1=0 fixed, |V3| fixed)
    %
    % J = [ dP2/d_del2  dP2/d_del3  dP2/d_V2 ]
    %     [ dP3/d_del2  dP3/d_del3  dP3/d_V2 ]
    %     [ dQ2/d_del2  dQ2/d_del3  dQ2/d_V2 ]

    J = zeros(3,3);

    % --- dP/d_delta (off-diagonal, i~=k) ---
    % dPi/d_delk = Vi*Vk*(Gik*sin(deli-delk) - Bik*cos(deli-delk))  i~=k
    % dPi/d_deli = -Qi - Bii*Vi^2
    % --- dQ/d_delta ---
    % dQi/d_delk = -Vi*Vk*(Gik*cos(deli-delk) + Bik*sin(deli-delk)) i~=k
    % dQi/d_deli =  Pi - Gii*Vi^2
    % --- dP/d_V, dQ/d_V ---
    % dPi/d_Vk = Vi*(Gik*cos(deli-delk) + Bik*sin(deli-delk))       i~=k
    % dPi/d_Vi = Pi/Vi + Gii*Vi
    % dQi/d_Vk = Vi*(Gik*sin(deli-delk) - Bik*cos(deli-delk))       i~=k
    % dQi/d_Vi = Qi/Vi - Bii*Vi

    % Row 1: dP2/d[del2, del3, V2]
    % dP2/d_del2
    J(1,1) = -Q_calc(2) - B(2,2)*V(2)^2;
    % dP2/d_del3
    J(1,2) = V(2)*V(3)*(G(2,3)*sin(del(2)-del(3)) - B(2,3)*cos(del(2)-del(3)));
    % dP2/d_V2 (multiply by V2 to use standard form d_P/d_|V| * |V|, then divide)
    J(1,3) = P_calc(2)/V(2) + G(2,2)*V(2);

    % Row 2: dP3/d[del2, del3, V2]
    % dP3/d_del2
    J(2,1) = V(3)*V(2)*(G(3,2)*sin(del(3)-del(2)) - B(3,2)*cos(del(3)-del(2)));
    % dP3/d_del3
    J(2,2) = -Q_calc(3) - B(3,3)*V(3)^2;
    % dP3/d_V2
    J(2,3) = V(3)*(G(3,2)*cos(del(3)-del(2)) + B(3,2)*sin(del(3)-del(2)));

    % Row 3: dQ2/d[del2, del3, V2]
    % dQ2/d_del2
    J(3,1) = P_calc(2) - G(2,2)*V(2)^2;
    % dQ2/d_del3
    J(3,2) = -V(2)*V(3)*(G(2,3)*cos(del(2)-del(3)) + B(2,3)*sin(del(2)-del(3)));
    % dQ2/d_V2
    J(3,3) = Q_calc(2)/V(2) - B(2,2)*V(2);

    %% Solve Jacobian equation: J * [d_del2; d_del3; d_V2] = mismatch
    correction = J \ mismatch;

    %% Update state variables
    del(2) = del(2) + correction(1);
    del(3) = del(3) + correction(2);
    V(2)   = V(2)   + correction(3);
    % Bus 3 voltage magnitude stays fixed; Bus 1 is slack
end

fprintf('%s\n\n', repmat('-',72,1));

%% Results
if converged
    fprintf('>>> Converged in %d iterations (tolerance = %g)\n\n', iter, tolerance);
else
    fprintf('>>> WARNING: Did not converge in %d iterations\n\n', max_iter);
end

%% Final Power Flow Calculations
for i = 1:3
    for k = 1:3
        P_calc(i) = P_calc(i) + V(i)*V(k)*(G(i,k)*cos(del(i)-del(k)) + B(i,k)*sin(del(i)-del(k)));
        Q_calc(i) = Q_calc(i) + V(i)*V(k)*(G(i,k)*sin(del(i)-del(k)) - B(i,k)*cos(del(i)-del(k)));
    end
end
% Recalculate cleanly
P_calc = zeros(3,1); Q_calc = zeros(3,1);
for i = 1:3
    for k = 1:3
        P_calc(i) = P_calc(i) + V(i)*V(k)*(G(i,k)*cos(del(i)-del(k)) + B(i,k)*sin(del(i)-del(k)));
        Q_calc(i) = Q_calc(i) + V(i)*V(k)*(G(i,k)*sin(del(i)-del(k)) - B(i,k)*cos(del(i)-del(k)));
    end
end

fprintf('============================================================\n');
fprintf('              BUS VOLTAGES (Final Solution)\n');
fprintf('============================================================\n');
fprintf('%-6s %-12s %-14s %-12s %-12s\n','Bus','|V| (pu)','Angle (deg)','P (pu)','Q (pu)');
fprintf('%s\n', repmat('-',60,1));
for i = 1:3
    fprintf('%-6d %-12.6f %-14.6f %-12.6f %-12.6f\n', i, V(i), rad2deg(del(i)), P_calc(i), Q_calc(i));
end
fprintf('%s\n\n', repmat('-',60,1));

%% Line Flows and Losses
fprintf('============================================================\n');
fprintf('              LINE FLOWS AND LINE LOSSES\n');
fprintf('============================================================\n');
fprintf('%-8s %-14s %-14s %-14s %-14s\n','Line','P_from(pu)','Q_from(pu)','P_to(pu)','Q_to(pu)');
fprintf('%-8s %-14s %-14s\n','','P_loss(pu)','Q_loss(pu)');
fprintf('%s\n', repmat('-',72,1));

lines = [1 2 z12; 1 3 z13; 2 3 z23];
total_Ploss = 0; total_Qloss = 0;

for l = 1:3
    i   = lines(l,1);
    k   = lines(l,2);
    zij = lines(l,3);
    yij = 1/zij;

    % Complex voltages
    Vi = V(i)*exp(1j*del(i));
    Vk = V(k)*exp(1j*del(k));

    % Current from i to k
    Iik = (Vi - Vk) * yij;
    Iki = -Iik;

    % Complex power from i to k
    Sik = Vi * conj(Iik);
    Ski = Vk * conj(Iki);

    Ploss = real(Sik) + real(Ski);
    Qloss = imag(Sik) + imag(Ski);
    total_Ploss = total_Ploss + Ploss;
    total_Qloss = total_Qloss + Qloss;

    fprintf('%-8s %-14.6f %-14.6f %-14.6f %-14.6f\n', ...
        sprintf('%d-%d',i,k), real(Sik), imag(Sik), real(Ski), imag(Ski));
    fprintf('%-8s %-14.6f %-14.6f\n','  Loss:', Ploss, Qloss);
    fprintf('%s\n', repmat('-',72,1));
end

fprintf('\nTotal Active Power Loss  : %.6f pu = %.4f MW\n', total_Ploss, total_Ploss*baseMVA);
fprintf('Total Reactive Power Loss: %.6f pu = %.4f Mvar\n\n', total_Qloss, total_Qloss*baseMVA);

%% Summary in MW and Mvar
fprintf('============================================================\n');
fprintf('          RESULTS IN PHYSICAL UNITS (100 MVA Base)\n');
fprintf('============================================================\n');
fprintf('%-6s %-12s %-14s %-14s %-14s\n','Bus','|V| (pu)','Angle (deg)','P (MW)','Q (Mvar)');
fprintf('%s\n', repmat('-',65,1));
for i = 1:3
    fprintf('%-6d %-12.6f %-14.6f %-14.4f %-14.4f\n', i, V(i), rad2deg(del(i)), P_calc(i)*baseMVA, Q_calc(i)*baseMVA);
end
fprintf('%s\n', repmat('-',65,1));