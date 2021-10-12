addpath("Helpers\")
format long

g = 9.81; %m/s²
%--- konstante Werte des Flugzeuges ---
%Struktur
AC.m = 100000; %kg
AC.S = 268; %m²
AC.I_y = 6.78*(10)^6; %kgm²
AC.l_my = 6.39; %m
%Antrieb
AC.F_TBPmax = 320800; %N
AC.n_rho = 0.7;
AC.n_v = 0.0;
AC.rho_TBP = 1.225;
%Aero - Auftrieb
AC.C_A_Alpha0Eta0 = 0.239;
AC.C_A_Alpha = 4.6;
AC.C_A_Eta = 0.189;
%Aero - Widerstand
AC.C_W0 = 0.015;
AC.k = 0.04;
%Aero - Nickmoment
AC.C_m_Alpha0Eta0 = 0.07;
AC.C_m_Alpha = -0.824;
AC.C_m_Eta = -0.638;
AC.C_m_q = -3.0;
AC.C_m_alphaPunkt = -1.0;

%--- Bezugsflugzustand ---
h = UnitConversion.ft2m(30000); % 30000ft
v_kt = 275; % kt (IAS)
v_ms = UnitConversion.kts2ms(v_kt); % m/s (IAS)

% Berechnung True Airspeed (TAS)
rho_ref = Aerodynamics.rho(h);
V_ref = Aerodynamics.ias2tas(v_ms,rho_ref);

C_A_ref = Aerodynamics.C_A_ref_initial(AC.m,rho_ref,V_ref,AC.S);
q_quer_ref = Aerodynamics.q_quer(V_ref,rho_ref);

for i = 1:10 %Iteriere
    eta_ref = Aerodynamics.eta_ref(AC.C_m_Alpha0Eta0, C_A_ref, AC.C_A_Alpha0Eta0, AC.C_m_Alpha, AC.C_A_Alpha, AC.C_m_Eta, AC.C_A_Eta);
    alpha_ref = Aerodynamics.alpha_ref(C_A_ref,AC.C_A_Alpha0Eta0,AC.C_A_Eta,eta_ref,AC.C_A_Alpha);
    C_W_ref = Aerodynamics.C_W_ref(AC.C_W0,AC.k,C_A_ref);
    W_ref = Aerodynamics.W_ref(C_W_ref,q_quer_ref,AC.S);

    F_ref = Aerodynamics.F_ref(W_ref,alpha_ref);
    delta_ref = Aerodynamics.delta_ref(F_ref,AC.F_TBPmax,rho_ref,AC.rho_TBP,AC.n_rho);
    A_ref = Aerodynamics.A_ref(AC.m,F_ref,alpha_ref);
    C_A_ref = Aerodynamics.C_A_ref(A_ref,q_quer_ref,AC.S);
end

fprintf("\nC_A = %1.15f\neta = %1.15f\nalpha = %1.15f\nC_W = %1.15f\nW = %1.15e\nA = %1.15e\nF_Ref = %1.15e\ndelta = %1.15f\n", ...
    C_A_ref,eta_ref,alpha_ref,C_W_ref,W_ref,A_ref,F_ref,delta_ref);


