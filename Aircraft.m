classdef Aircraft
    properties

        %--- konstante Werte des Flugzeuges ---
        %Struktur
        m = 100000; %kg
        S = 268; %m²
        I_y = 6.78*(10)^6; %kgm²
        l_my = 6.39; %m
        %Antrieb
        F_TBPmax = 320800; %N
        n_rho = 0.7;
        n_v = 0.0;
        rho_TBP = 1.225;
        %Aero - Auftrieb
        C_A_Alpha0Eta0 = 0.239;
        C_A_Alpha = 4.6;
        C_A_Eta = 0.189;
        %Aero - Widerstand
        C_W0 = 0.015;
        k = 0.04;
        %Aero - Nickmoment
        C_m_Alpha0Eta0 = 0.07;
        C_m_Alpha = -0.824;
        C_m_Eta = -0.638;
        C_m_q = -3.0;
        C_m_alphaPunkt = -1.0;
    end
end