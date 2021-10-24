addpath("Helpers\")
addpath("Structs\")
format long

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


