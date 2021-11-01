addpath("Helpers/")
addpath("Structs/")
format long
sympref('FloatingPointOutput',true);

run("Structs/Aircraft.m");
run("Structs/Environment.m");

input=[30000, 275; 30000, 325; 25000, 250; 25000, 300];

for i = 1:1%length(input)
    %% --- Bezugsflugzustand ---
    %Laden der vorgegebenen Parameter aus VAR
    BFZ.h = UnitConversion.ft2m(input(i,1)); % 30000ft
    BFZ.v_kt =input(i,2); % kt (IAS)
    BFZ.v_ms = UnitConversion.kts2ms(BFZ.v_kt); % m/s (IAS)
    
    % Berechnung True Airspeed (TAS)
    BFZ.rho = Aero.rho(BFZ.h,ENV);
    BFZ.V = Aero.ias2tas(BFZ.v_ms,BFZ.rho);
    
    BFZ.C_A = Aero.C_A_ref_initial(AC,BFZ);
    BFZ.q_quer = Aero.q_quer(BFZ);
    
    for j = 1:10 %Iteriere
        BFZ.eta= Aero.eta_ref(AC,BFZ);
        BFZ.alpha = Aero.alpha_ref(AC,BFZ);
        BFZ.C_W= Aero.C_W_ref(AC,BFZ);
        BFZ.W= Aero.W_ref(AC,BFZ);
    
        BFZ.F= Aero.F_ref(BFZ);
        BFZ.delta= Aero.delta_ref(AC,BFZ);
        BFZ.A= Aero.A_ref(AC,BFZ);
        BFZ.C_A= Aero.C_A_ref(AC,BFZ);
    end
    
    %% --- Berechnung der Ersatzgrößen ---
    BFZ.X_V = Aero.X_V(AC,BFZ);
    BFZ.X_alpha = Aero.X_alpha(AC,BFZ);
    BFZ.X_eta = Aero.X_eta(AC,BFZ);
    BFZ.X_delta = Aero.X_delta(AC,BFZ);
    
    BFZ.Z_V = Aero.Z_V(AC,BFZ);
    BFZ.Z_alpha = Aero.Z_alpha(AC,BFZ);
    BFZ.Z_eta = Aero.Z_eta(AC,BFZ);
    BFZ.Z_delta = Aero.Z_delta(BFZ);
    
    BFZ.M_alpha = Aero.M_alpha(AC,BFZ);
    BFZ.M_eta = Aero.M_eta(AC,BFZ);
    BFZ.M_delta = Aero.M_delta(AC,BFZ);
    BFZ.M_q = Aero.M_q(AC,BFZ);
    BFZ.M_V = Aero.M_V(AC,BFZ);

    %% -- Zustandsraumdarstellung ---
    %-- 2x2 Näherung --
    %Anstellwinkelschwingung
    AS.A = [BFZ.M_q BFZ.M_alpha ; 1 BFZ.Z_alpha];
    AS.B = [BFZ.M_eta BFZ.M_delta ; BFZ.Z_eta BFZ.Z_delta]; 
    AS.Sys = ss(AS.A,AS.B,eye(2,2),0,'StateName',{'Δq','Δα'},'InputName',{'Δη','Δδ'},'OutputName',{'Δq','Δα'});
    AS.Eig = eig(AS.A);
    AS.TF = tf(AS.Sys);
    
    AS.sigma = 0.5*(BFZ.M_q+BFZ.Z_alpha);
    AS.omega_0 = sqrt(BFZ.M_q*BFZ.Z_alpha-BFZ.M_alpha);
    AS.omega = ZRM.omega(AS.omega_0,AS.sigma);
    AS.D = ZRM.D(AS.sigma,AS.omega_0);
    AS.T = ZRM.omega2T(AS.omega);
    
    %Bahnschwingung
    BS.A = [BFZ.X_V -ENV.g ; -BFZ.Z_V 0];
    BS.B = [BFZ.X_eta BFZ.X_delta ; -BFZ.Z_eta -BFZ.Z_delta];
    BS.Sys = ss(BS.A,BS.B,eye(2,2),0,'StateName',{'ΔV','Δγ'},'InputName',{'Δη','Δδ'},'OutputName',{'ΔV','Δγ'});
    BS.Eig = eig(BS.A);
    BS.TF = tf(BS.Sys);
    
    BS.sigma = BFZ.X_V/2;
    BS.omega_0 = sqrt(-ENV.g*BFZ.Z_V);
    BS.omega = ZRM.omega(BS.omega_0,BS.sigma);
    BS.D = ZRM.D(BS.sigma,BS.omega_0);
    BS.T = ZRM.omega2T(BS.omega);
    
    %--- 4x4 Lösung / lineares Differntialgleichungssystem --
    DGS.A = [BFZ.X_V -ENV.g 0 BFZ.X_alpha 
            -BFZ.Z_V 0 0 -BFZ.Z_alpha
            BFZ.M_V 0 BFZ.M_q BFZ.M_alpha
            BFZ.Z_V 0 1 BFZ.Z_alpha];
    DGS.B = [BFZ.X_eta BFZ.X_delta ; -BFZ.Z_eta -BFZ.Z_delta ; BFZ.M_eta BFZ.M_delta ; BFZ.Z_eta BFZ.Z_delta];
    DGS.Sys = ss(DGS.A,DGS.B,eye(4,4),0,'StateName',{'ΔV','Δγ','Δq','Δα'},'InputName',{'Δη','Δδ'},'OutputName',{'ΔV','Δγ','Δq','Δα'});
    DGS.Eig = eig(DGS.A);
    DGS.TF = tf(DGS.Sys);
    
    DGS.sigma_as = real(DGS.Eig(1));
    DGS.omega_as = abs(imag(DGS.Eig(1)));
    DGS.omega_0_as = sqrt(DGS.sigma_as^2+DGS.omega_as^2);
    DGS.D_as = ZRM.D(DGS.sigma_as,DGS.omega_0_as);
    DGS.T_as = ZRM.omega2T(DGS.omega_as);
    
    DGS.sigma_bs = real(DGS.Eig(3));
    DGS.omega_bs = abs(imag(DGS.Eig(3)));
    DGS.omega_0_bs = sqrt(DGS.sigma_bs^2+DGS.omega_bs^2);
    DGS.D_bs = ZRM.D(DGS.sigma_bs,DGS.omega_0_bs);
    DGS.T_bs = ZRM.omega2T(DGS.omega_bs);
    
    %% --- Proportionalrückführung ---
    [AS.TF_Z,AS.TF_N] = tfdata(AS.Sys);
    NickD.N = AS.TF_N(1,1);
    NickD.Z = AS.TF_Z(1,1);
    syms k s;
    
    NickD.sigma = -(1/2)*(NickD.N{1}(2)+k*NickD.Z{1}(2));
    NickD.omega_0 = sqrt(NickD.N{1}(3)+k*NickD.Z{1}(3));
    NickD.D = -NickD.sigma/NickD.omega_0;
    
    NickD.k_p = solve(NickD.D==1/sqrt(2),k);

    %% --- Lag-Filter Rückführung ---
    [NickD.NST,NickD.P,NickD.k0]=zpkdata(AS.TF(1,1),'v');
    NickD.Lag.d = -1;

    %Konstruieren von s_soll
    NickD.Lag.D = 1/sqrt(2);
    NickD.Lag.omega_0 = 3;
    NickD.Lag.sigma=-NickD.Lag.D*NickD.Lag.omega_0;
    NickD.Lag.omega=sqrt(NickD.Lag.omega_0^2 - NickD.Lag.sigma^2);
    NickD.Lag.s_soll = complex(NickD.Lag.sigma,NickD.Lag.omega);
    
    %Anwendung Betrags- und Phasenbedingung
    syms psi_c k;
    NickD.Lag.psi_c = solve(pi+sum(angle(NickD.Lag.s_soll-NickD.NST))- ...
        sum(angle(NickD.Lag.s_soll-NickD.P))+pi+psi_c- ...
        angle(NickD.Lag.s_soll-NickD.Lag.d)== pi,psi_c);
    NickD.Lag.c = real(NickD.Lag.s_soll) - imag(NickD.Lag.s_soll) / tan(NickD.Lag.psi_c);

    NickD.Lag.k = solve(abs(NickD.k0)*prod(abs(NickD.Lag.s_soll-NickD.NST))/ ...
        prod(abs(NickD.Lag.s_soll-NickD.P))*abs(k)* ...
        abs(NickD.Lag.s_soll-NickD.Lag.c)/abs(NickD.Lag.s_soll-NickD.Lag.d)==1,k);
    if(NickD.Lag.k > 0)
        NickD.Lag.k=NickD.Lag.k*-1;
    end
end
