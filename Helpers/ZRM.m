classdef ZRM
    methods(Static)
        function [omega] = omega(omega_0,sigma)
            %Kreisfrequenz
            omega=sqrt(omega_0^2-sigma^2);
        end
        function [D] = D(sigma,omega_0)
            %relative Dämpfung
            D = -sigma/omega_0;
        end
        function [T] = omega2T(omega)
            %Periodendauer aus Kreisfrequenz
            T = 2*pi/omega;
        end 
    end
end
