addpath("Helpers\")
addpath("Structs\")
format long

run("Structs\Aircraft.m");
run("Structs\Environment.m");

%--- Bezugsflugzustand ---
h = UnitConversion.ft2m(30000); % 30000ft
v_kt = 275; % kt (IAS)
v_ms = UnitConversion.kts2ms(v_kt); % m/s (IAS)

% Berechnung True Airspeed (TAS)
BFZ.rho = Aero.rho(h,ENV);
BFZ.V = Aero.ias2tas(v_ms,BFZ.rho);

BFZ.C_A = Aero.C_A_ref_initial(AC,BFZ);
BFZ.q_quer = Aero.q_quer(BFZ);

for i = 1:10 %Iteriere
    BFZ.eta= Aero.eta_ref(AC,BFZ);
    BFZ.alpha = Aero.alpha_ref(AC,BFZ);
    BFZ.C_W= Aero.C_W_ref(AC,BFZ);
    BFZ.W= Aero.W_ref(AC,BFZ);

    BFZ.F= Aero.F_ref(BFZ);
    BFZ.delta= Aero.delta_ref(AC,BFZ);
    BFZ.A= Aero.A_ref(AC,BFZ);
    BFZ.C_A= Aero.C_A_ref(AC,BFZ);
end

%--- Berechnung der Ersatzgrößen ---
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


%--- Zustandsraumdarstellung ---
%-- 2x2 Näherung --
%Anstellwinkelschwingung
AS.A = [BFZ.M_q BFZ.M_alpha ; 1 BFZ.Z_alpha];
AS.B = [BFZ.M_eta BFZ.M_delta ; BFZ.Z_delta BFZ.Z_delta]; 
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
BS.B = [BFZ.X_eta BFZ.X_delta ; -BFZ.Z_delta -BFZ.Z_delta];
BS.Sys = ss(BS.A,BS.B,eye(2,2),0,'StateName',{'ΔV','Δγ'},'InputName',{'Δη','Δδ'},'OutputName',{'ΔV','Δγ'});
BS.Eig = eig(BS.A);
BS.TF = tf(BS.Sys);

BS.sigma = BFZ.X_V/2;
BS.omega_0 = sqrt(-ENV.g*BFZ.Z_V);
BS.omega = ZRM.omega(BS.omega_0,BS.sigma);
BS.D = ZRM.D(BS.sigma,BS.omega_0);
BS.T = ZRM.omega2T(BS.omega);

%-- 4x4 Lösung / lineares Differntialgleichungssystem --
DGS.A = [BFZ.X_V -ENV.g 0 BFZ.X_alpha 
        -BFZ.Z_V 0 0 -BFZ.alpha
        BFZ.M_V 0 BFZ.M_q BFZ.M_alpha
        BFZ.Z_V 0 1 BFZ.Z_alpha];
DGS.B = [BFZ.X_eta BFZ.X_delta ; -BFZ.Z_eta -BFZ.Z_delta ; BFZ.M_eta BFZ.M_delta ; BFZ.Z_eta BFZ.Z_delta];
DGS.Sys = ss(DGS.A,DGS.B,eye(4,4),0,'StateName',{'ΔV','Δγ','Δq','Δα'},'InputName',{'Δη','Δδ'},'OutputName',{'ΔV','Δγ','Δq','Δα'});
DGS.Eig = eig(DGS.A);
DGS.TF = tf(DGS.Sys);

