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