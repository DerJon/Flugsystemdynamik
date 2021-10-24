classdef Aero
    methods(Static)
        function [rho] = rho(height,ENV)
            %berechne Dichte aus Höhe (Standard-Atmosphäre)
            rho = ENV.rho_0*(1 +(ENV.gamma_H/ENV.T_0)*height)^(-(ENV.g/(ENV.R*ENV.gamma_H)) - 1);
        end

        function [q_quer] = q_quer(BFZ)
            %berechnet q_quer aus TrueAirspeed(v) und Luftdichte (rho)
            q_quer = BFZ.rho/2 * BFZ.V^2;
        end

        function [tas] = ias2tas(ias, rho)
            %Umrechnung des indicated Airspeed (IAS) in true Airspeed
            %(TAS)
            tas = sqrt(1.225/rho) * ias;
        end
        
        function [C_A_ref] = C_A_ref_initial(AC,BFZ)
            %Berechnung von C_Aref aus dem vertikalen Gleichgewicht
            %Verwendung für die initale Bestimmung
            C_A_ref = (2*AC.m*9.81)/(BFZ.rho*(BFZ.V^2)*AC.S);
        end

        function [C_A_ref] = C_A_ref(AC,BFZ)
            %Berechnung von C_Aref aus Trimmgleichungen
            %Verwendung für Iteration
            C_A_ref = BFZ.A/(BFZ.q_quer*AC.S);
        end

        function [eta_ref] = eta_ref(AC,BFZ)
            %Berechnung des Höhenruderausschlags aus dem vertikalen GGW &
            %dem Bezugsflugzustand
            eta_ref = -((AC.C_m_Alpha0Eta0+(BFZ.C_A-AC.C_A_Alpha0Eta0)*(AC.C_m_Alpha/AC.C_A_Alpha))/(AC.C_m_Eta-(AC.C_A_Eta*(AC.C_m_Alpha/AC.C_A_Alpha))));
        end

        function [alpha_ref] = alpha_ref(AC,BFZ)
            %Berechnung des Bezugs-Anstellwinkels 
            alpha_ref = (BFZ.C_A-AC.C_A_Alpha0Eta0-(AC.C_A_Eta*BFZ.eta))/AC.C_A_Alpha;
        end

        function [A_ref] = A_ref(AC,BFZ)
            %Berechnung des Bezugs-Auftriebs
            A_ref = AC.m*9.81-BFZ.F*sin(BFZ.alpha);
        end

        function [C_W_ref] = C_W_ref(AC,BFZ)
            %Berechnung Bezugs-Widerstandsbeiwerts
            C_W_ref = AC.C_W0+(AC.k*(BFZ.C_A)^2);
        end

        function [W_ref] = W_ref(AC,BFZ)
            %Berechnung Bezugs-Widerstand
            W_ref = BFZ.C_W*BFZ.q_quer*AC.S;
        end

        function [F_ref] = F_ref(BFZ)
            %Berechnung Bezugs-Schub
            F_ref = BFZ.W / cos(BFZ.alpha);
        end

        function [delta_ref] = delta_ref(AC,BFZ)
            %Berechnung Bezugs-Schubdrosselstellung
                delta_ref = (BFZ.F/AC.F_TBPmax)*((BFZ.rho/AC.rho_TBP)^(-AC.n_rho));
        end
        
        %Ersatzgrößen der Widerstandsgleichung
        function [X_V] = X_V(AC,BFZ)
            X_V = -((BFZ.q_quer*AC.S)/(AC.m*BFZ.V))*(2-AC.n_v)*BFZ.C_W;
        end

        function [X_alpha] = X_alpha(AC,BFZ)
            AC.C_W_alpha = 2*AC.k*BFZ.C_A*AC.C_A_Alpha;
            X_alpha = -((BFZ.q_quer*AC.S)/AC.m)*AC.C_W_alpha;
        end

        function [X_eta] = X_eta(AC,BFZ)
            AC.C_W_Eta = 2*AC.k*BFZ.C_A*AC.C_A_Eta;
            X_eta = -((BFZ.q_quer*AC.S)/AC.m)*AC.C_W_Eta;
        end

        function [X_delta] = X_delta(AC,BFZ)
            X_delta = ((BFZ.q_quer*AC.S)/AC.m)*(BFZ.C_W/BFZ.delta);
        end

        %Ersatzgrößen der Auftriebsgleichung
        function [Z_V] = Z_V(AC,BFZ)
             Z_V = -2*((BFZ.q_quer*AC.S)/(AC.m*BFZ.V))*(BFZ.C_A/BFZ.V);
        end

        function [Z_alpha] = Z_alpha(AC,BFZ)
            Z_alpha =  -((BFZ.q_quer*AC.S)/(AC.m*BFZ.V))*AC.C_A_Alpha;
        end

        function [Z_eta] = Z_eta(AC,BFZ)
             Z_eta = -((BFZ.q_quer*AC.S)/(AC.m*BFZ.V))*AC.C_A_Eta;
        end

        function [Z_delta] = Z_delta(BFZ)
            Z_delta = -BFZ.X_delta*(tan(BFZ.alpha)/BFZ.V);
        end

        %Ersatzgrößen der Momentengleichung
        function [M_V] = M_V(AC,BFZ)
            M_V = BFZ.Z_V*((BFZ.q_quer*AC.S*AC.l_my)/AC.I_y)*(AC.l_my/BFZ.V)*AC.C_m_alphaPunkt;
        end

        function [M_q] = M_q(AC,BFZ)
            M_q = ((BFZ.q_quer*AC.S*AC.l_my)/AC.I_y)*(AC.l_my/BFZ.V)*(AC.C_m_q+AC.C_m_alphaPunkt);
        end

        function [M_alpha] = M_alpha(AC,BFZ)
            M_alpha = ((BFZ.q_quer*AC.S*AC.l_my)/AC.I_y)*(AC.C_m_Alpha+BFZ.Z_alpha*(AC.l_my/BFZ.V)*AC.C_m_alphaPunkt);
        end

        function [M_eta] = M_eta(AC,BFZ)
            M_eta =  ((BFZ.q_quer*AC.S*AC.l_my)/AC.I_y)*(AC.C_m_Eta+BFZ.Z_eta*(AC.l_my/BFZ.V)*AC.C_m_alphaPunkt);
        end

        function [M_delta] = M_delta(AC,BFZ)
            M_delta = BFZ.Z_delta*((BFZ.q_quer*AC.S*AC.l_my)/AC.I_y)*(AC.l_my/BFZ.V)*AC.C_m_alphaPunkt;
        end
 
    end
end

