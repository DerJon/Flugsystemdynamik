addpath("Helpers\")
addpath("Structs\")
format long
sympref('FloatingPointOutput',true);

run("Structs\Aircraft.m");
run("Structs\Environment.m");

%Preallocation of Structs
BFZ(1:4) = struct();
AS(1:4) = struct();
BS(1:4) = struct();
DGS(1:4) = struct();
NickD(1:4) = struct();

input=[30000, 275; 30000, 325; 25000, 250; 25000, 300];

for i = 1:length(input)
    %% --- Bezugsflugzustand ---
    %Laden der vorgegebenen Parameter aus VAR
    BFZ(i).h = UnitConversion.ft2m(input(i,1)); % 30000ft
    BFZ(i).v_kt =input(i,2); % kt (IAS(i))
    BFZ(i).v_ms = UnitConversion.kts2ms(BFZ(i).v_kt); % m/s (IAS(i))
    
    % Berechnung True Airspeed (TAS(i))
    BFZ(i).rho = Aero.rho(BFZ(i).h,ENV);
    BFZ(i).V = Aero.ias2tas(BFZ(i).v_ms,BFZ(i).rho);
    
    BFZ(i).C_A = Aero.C_A_ref_initial(AC,BFZ(i));
    BFZ(i).q_quer = Aero.q_quer(BFZ(i));
    
    for j = 1:10 %Iteriere
        BFZ(i).eta= Aero.eta_ref(AC,BFZ(i));
        BFZ(i).alpha = Aero.alpha_ref(AC,BFZ(i));
        BFZ(i).C_W= Aero.C_W_ref(AC,BFZ(i));
        BFZ(i).W= Aero.W_ref(AC,BFZ(i));
    
        BFZ(i).F= Aero.F_ref(BFZ(i));
        BFZ(i).delta= Aero.delta_ref(AC,BFZ(i));
        BFZ(i).A= Aero.A_ref(AC,BFZ(i));
        BFZ(i).C_A= Aero.C_A_ref(AC,BFZ(i));
    end
    
    %% --- Berechnung der Ersatzgrößen ---
    BFZ(i).X_V = Aero.X_V(AC,BFZ(i));
    BFZ(i).X_alpha = Aero.X_alpha(AC,BFZ(i));
    BFZ(i).X_eta = Aero.X_eta(AC,BFZ(i));
    BFZ(i).X_delta = Aero.X_delta(AC,BFZ(i));
    
    BFZ(i).Z_V = Aero.Z_V(AC,BFZ(i));
    BFZ(i).Z_alpha = Aero.Z_alpha(AC,BFZ(i));
    BFZ(i).Z_eta = Aero.Z_eta(AC,BFZ(i));
    BFZ(i).Z_delta = Aero.Z_delta(BFZ(i));
    
    BFZ(i).M_alpha = Aero.M_alpha(AC,BFZ(i));
    BFZ(i).M_eta = Aero.M_eta(AC,BFZ(i));
    BFZ(i).M_delta = Aero.M_delta(AC,BFZ(i));
    BFZ(i).M_q = Aero.M_q(AC,BFZ(i));
    BFZ(i).M_V = Aero.M_V(AC,BFZ(i));

    %% -- Zustandsraumdarstellung ---
    %-- 2x2 Näherung --
    %Anstellwinkelschwingung
    AS(i).A = [BFZ(i).M_q BFZ(i).M_alpha ; 1 BFZ(i).Z_alpha];
    AS(i).B = [BFZ(i).M_eta BFZ(i).M_delta ; BFZ(i).Z_eta BFZ(i).Z_delta]; 
    AS(i).Sys = ss(AS(i).A,AS(i).B,eye(2,2),0,'StateName',{'Δq','Δα'},'InputName',{'Δη','Δδ'},'OutputName',{'Δq','Δα'});
    AS(i).Eig = eig(AS(i).A);
    AS(i).TF = tf(AS(i).Sys);
    
    AS(i).sigma = 0.5*(BFZ(i).M_q+BFZ(i).Z_alpha);
    AS(i).omega_0 = sqrt(BFZ(i).M_q*BFZ(i).Z_alpha-BFZ(i).M_alpha);
    AS(i).omega = ZRM.omega(AS(i).omega_0,AS(i).sigma);
    AS(i).D = ZRM.D(AS(i).sigma,AS(i).omega_0);
    AS(i).T = ZRM.omega2T(AS(i).omega);
    
    %Bahnschwingung
    BS(i).A = [BFZ(i).X_V -ENV.g ; -BFZ(i).Z_V 0];
    BS(i).B = [BFZ(i).X_eta BFZ(i).X_delta ; -BFZ(i).Z_eta -BFZ(i).Z_delta];
    BS(i).Sys = ss(BS(i).A,BS(i).B,eye(2,2),0,'StateName',{'ΔV','Δγ'},'InputName',{'Δη','Δδ'},'OutputName',{'ΔV','Δγ'});
    BS(i).Eig = eig(BS(i).A);
    BS(i).TF = tf(BS(i).Sys);
    
    BS(i).sigma = BFZ(i).X_V/2;
    BS(i).omega_0 = sqrt(-ENV.g*BFZ(i).Z_V);
    BS(i).omega = ZRM.omega(BS(i).omega_0,BS(i).sigma);
    BS(i).D = ZRM.D(BS(i).sigma,BS(i).omega_0);
    BS(i).T = ZRM.omega2T(BS(i).omega);
    
    %--- 4x4 Lösung / lineares Differntialgleichungssystem --
    DGS(i).A = [BFZ(i).X_V -ENV.g 0 BFZ(i).X_alpha 
            -BFZ(i).Z_V 0 0 -BFZ(i).Z_alpha
            BFZ(i).M_V 0 BFZ(i).M_q BFZ(i).M_alpha
            BFZ(i).Z_V 0 1 BFZ(i).Z_alpha];
    DGS(i).B = [BFZ(i).X_eta BFZ(i).X_delta ; -BFZ(i).Z_eta -BFZ(i).Z_delta ; BFZ(i).M_eta BFZ(i).M_delta ; BFZ(i).Z_eta BFZ(i).Z_delta];
    DGS(i).Sys = ss(DGS(i).A,DGS(i).B,eye(4,4),0,'StateName',{'ΔV','Δγ','Δq','Δα'},'InputName',{'Δη','Δδ'},'OutputName',{'ΔV','Δγ','Δq','Δα'});
    DGS(i).Eig = eig(DGS(i).A);
    DGS(i).TF = tf(DGS(i).Sys);
    
    DGS(i).sigma_AS(i) = real(DGS(i).Eig(1));
    DGS(i).omega_AS(i) = abs(imag(DGS(i).Eig(1)));
    DGS(i).omega_0_AS(i) = sqrt(DGS(i).sigma_AS(i)^2+DGS(i).omega_AS(i)^2);
    DGS(i).D_AS(i) = ZRM.D(DGS(i).sigma_AS(i),DGS(i).omega_0_AS(i));
    DGS(i).T_AS(i) = ZRM.omega2T(DGS(i).omega_AS(i));
    
    DGS(i).sigma_BS(i) = real(DGS(i).Eig(3));
    DGS(i).omega_BS(i) = abs(imag(DGS(i).Eig(3)));
    DGS(i).omega_0_BS(i) = sqrt(DGS(i).sigma_BS(i)^2+DGS(i).omega_BS(i)^2);
    DGS(i).D_BS(i) = ZRM.D(DGS(i).sigma_BS(i),DGS(i).omega_0_BS(i));
    DGS(i).T_BS(i) = ZRM.omega2T(DGS(i).omega_BS(i));
    
    %% --- Proportionalrückführung ---
    [AS(i).TF_Z,AS(i).TF_N] = tfdata(AS(i).Sys);
    NickD(i).N = AS(i).TF_N(1,1);
    NickD(i).Z = AS(i).TF_Z(1,1);
    syms k s;
    
    NickD(i).sigma = -(1/2)*(NickD(i).N{1}(2)+k*NickD(i).Z{1}(2));
    NickD(i).omega_0 = sqrt(NickD(i).N{1}(3)+k*NickD(i).Z{1}(3));
    NickD(i).D = -NickD(i).sigma/NickD(i).omega_0;
    
    NickD(i).k_p = solve(NickD(i).D==1/sqrt(2),k);

    %% --- Lag-Filter Rückführung ---
    [NickD(i).NST,NickD(i).P,NickD(i).k0]=zpkdata(AS(i).TF(1,1),'v');
    NickD(i).Lag.d = -1;

    %Konstruieren von s_soll
    NickD(i).Lag.D = 1/sqrt(2);
    NickD(i).Lag.omega_0 = 3;
    NickD(i).Lag.sigma=-NickD(i).Lag.D*NickD(i).Lag.omega_0;
    NickD(i).Lag.omega=sqrt(NickD(i).Lag.omega_0^2 - NickD(i).Lag.sigma^2);
    NickD(i).Lag.s_soll = complex(NickD(i).Lag.sigma,NickD(i).Lag.omega);
    
    %Anwendung Betrags- und PhAS(i)enbedingung
    syms psi_c k;
    NickD(i).Lag.psi_c = solve(pi+sum(angle(NickD(i).Lag.s_soll-NickD(i).NST))- ...
        sum(angle(NickD(i).Lag.s_soll-NickD(i).P))+pi+psi_c- ...
        angle(NickD(i).Lag.s_soll-NickD(i).Lag.d)== pi,psi_c);
    NickD(i).Lag.c = real(NickD(i).Lag.s_soll) - imag(NickD(i).Lag.s_soll) / tan(NickD(i).Lag.psi_c);

    NickD(i).Lag.k = solve(abs(NickD(i).k0)*prod(abs(k)* ...
        abs(NickD(i).Lag.s_soll-NickD(i).Lag.c)/abs(NickD(i).Lag.s_soll-NickD(i).Lag.d))==1,k);
    if(NickD(i).Lag.k > 0)
        NickD(i).Lag.k=NickD(i).Lag.k*-1;
    end
end
