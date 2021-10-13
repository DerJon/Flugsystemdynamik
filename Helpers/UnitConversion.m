classdef UnitConversion
    methods(Static)
        %Geschwindigkeit
        function [v_ms] = kts2ms(v_kts)
            v_ms=(1852/3600)*v_kts; 
        end
        function [v_kts] = ms2kts(v_ms)
            v_kts = (3600/1852) * v_ms;
        end
        
        %Entfernung
        function [m] = ft2m(ft)
            m = ft * 0.30479;
        end
        function [ft] = m2ft(m)
            ft = m / 0.30479;
        end

    end
end