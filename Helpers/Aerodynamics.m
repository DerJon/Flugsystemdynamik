classdef Aerodynamics
    properties
        g = 9.81; % m/s²
        rho_0 = 1.225; % kg/m³
        R = 287.058; % J/kgK
        T_0 = 288.15; % K
        gamma_H = -0.0065; % K/m
    end
    methods(Static)
        function [rho] = rho(height)
            %berechne Dichte aus Höhe (Standard-Atmosphäre)
            g = 9.81; % m/s²
            rho_0 = 1.225; % kg/m³
            R = 287.058; % J/kgK
            T_0 = 288.15; % K
            gamma_H = -0.0065; % K/m

            rho = rho_0*(1 +(gamma_H/T_0)*height)^(-(g/(R*gamma_H)) - 1);
        end

        function [q_quer] = q_quer(v, rho)
            %berechnet q_quer aus TrueAirspeed(v) und Luftdichte (rho)
            q_quer = rho/2 * v^2;
        end

        function [tas] = ias2tas(ias, rho)
            %Umrechnung des IAS in TAS
            tas = sqrt(1.225/rho) * ias;
        end
        
        function [C_A_ref] = C_A_ref_initial(m,rho_ref,V_ref,S)
            %Berechnung von C_Aref aus dem vertikalen Gleichgewicht
            %Verwendung für die initale Bestimmung
            C_A_ref = (2*m*9.81)/(rho_ref*(V_ref^2)*S);
        end

        function [C_A_ref] = C_A_ref(A_ref, q_quer_ref,S)
            %Berechnung von C_Aref aus Trimmgleichungen
            %Verwendung für Iteration
            C_A_ref = A_ref/(q_quer_ref*S);
        end

        function [eta_ref] = eta_ref(C_M_Alpha0Eta0,C_A_ref,C_A_Alpha0Eta0,C_m_Alpha,C_A_Alpha,C_m_Eta,C_A_Eta)
            %Berechnung des Höhenruderausschlags aus dem vertikalen GGW &
            %dem Bezugsflugzustand
            eta_ref = -((C_M_Alpha0Eta0+(C_A_ref-C_A_Alpha0Eta0)*(C_m_Alpha/C_A_Alpha))/(C_m_Eta-(C_A_Eta*(C_m_Alpha/C_A_Alpha))));
        end

        function [alpha_ref] = alpha_ref(C_A_ref,C_A_Alpha0Eta0,C_A_Eta,eta_ref,C_A_Alpha)
            %Berechnung des Bezugs-Anstellwinkels 
            alpha_ref = (C_A_ref-C_A_Alpha0Eta0-(C_A_Eta*eta_ref))/C_A_Alpha;
        end

        function [A_ref] = A_ref(m,F_ref,alpha_ref)
            %Berechnung des Bezugs-Auftriebs
            A_ref = m*9.81-F_ref*sin(alpha_ref);
        end

        function [C_W_ref] = C_W_ref(C_W0,k,C_A_ref)
            %Berechnung Bezugs-Widerstandsbeiwerts
            C_W_ref = C_W0+(k*(C_A_ref)^2);
        end

        function [W_ref] = W_ref(C_W_ref,q_quer_ref,S)
            %Berechnung Bezugs-Widerstand
            W_ref = C_W_ref*q_quer_ref*S;
        end

        function [F_ref] = F_ref(W_ref,alpha_ref)
            %Berechnung Bezugs-Schub
            F_ref = W_ref / cos(alpha_ref);
        end

        function [delta_ref] = delta_ref(F_ref,F_TBP_max,rho_ref,rho_TBP,n_rho)
            %Berechnung Bezugs-Schubdrosselstellung
                delta_ref = (F_ref/F_TBP_max)/((rho_ref/rho_TBP)^(-n_rho));
        end
    end
end

