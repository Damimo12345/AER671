
Hauteur = 25000 ; %On peux le modifier de 25000 ft à 29000 ft
rho = density(Hauteur) ;

v_initial = 265*1.687811; %en ft/sec pour le calcul de W/s cruise pour avoir le vrai V
M = v_initial / sqrt(1.4 * 1716 * tempatmstd(Hauteur)*1.8); % Calculate Mach number


C_d0 = 0.020 ; %Le range change BCP quand je change le C_d0
e = 0.78 ;
A = 6;
k = 1/(pi*A*e) ;


reserve_fuel = 0.05 ;
trapped_fuel= 0.01;
W_payload = 500 + 2*180 ; %en lbs

ws_cruising = ws_cruise(C_d0,A,Hauteur,M,'prop');
V = sqrt((2/rho)*ws_cruising)*((k/C_d0)^0.25); %ft/s
V_kts = V/1.687811; % Vitesse en kts

L_D_max = sqrt(1/(4*C_d0*k));
C_BHP=0.45;
rendement = 0.8 ;
C = consommation_spec(C_BHP,V,rendement);


% À CHANGER À CHANQUE ITÉRATION (faire une boucle)
Range_start= 0; %nm commencer à 600nm
Range_iteration = 600 ;

while (abs(Range_start-Range_iteration) > 0.1) 
    
    Range_start = Range_iteration ;

    [W_TO,W_fuel,W_empty] = itertow('aviation-gen(1mot)',M,Hauteur,A,C,0,reserve_fuel,trapped_fuel,W_payload,Range_start); %Aviation-gen(1mot) Je suis pas sûre que ça soit ça (on a prop 1 mot, mais il donne slm prop 2moteurs)
    
    %CALCULS POUR ITÉRATION SUR R (on recalcule R avec W_to et W_fuel)
    
    Wla_Wto = 1- ((W_fuel/W_TO) / (1 + reserve_fuel + trapped_fuel));
    W1_Wto = 0.975 ;
    W2_W1 = 0.975;
    Wla_W_3 = 0.995;
    W2_W3 = (W1_Wto*W2_W1*Wla_W_3)/(Wla_Wto) ;
    
    %Range à remettre en tant que "RANGE_START"

    ws_cruising = ws_cruise(C_d0,A,Hauteur,M,'prop');
    V = sqrt((2/rho)*ws_cruising)*((k/C_d0)^0.25); %ft/s
    V_kts = V/1.687811; % Vitesse en kts
    M = v_initial / sqrt(1.4 * 1716 * tempatmstd(Hauteur)*1.8);
    
    Range_iteration = (V_kts/C)*L_D_max*log(W2_W3); %Ln c'est log sur MatLab
end

%%%%%%   TROUVER S %%%%%%%

S = W_TO/ws_cruising; %en FT^2 , À MODIFIER

% Vérification pour le Distance de décollage et d'atterissage;
%Est-que on peux juste prendre le W_to final avec S qu' en prenant une
%valeur historique pour T/W

rhosl=density(0);

sigma = rho/rhosl;

%%%%%%%%%%%%%%% TAKE OFF DISTANCE %%%%%%%%%%%%%%%%

% s_TO_start = 2500 ; % ft, on choisi au pif !!!!!

C_l_TOmax = 1.7 ; % AVEC RAYMER
T_W = 0.5 ; %Donnée historique (0,25-0,57)
W_S_TO = ws_cruising;

TOP = W_S_TO* (1/C_l_TOmax)*(1/T_W)*(1/sigma) ; 

takeoffDistance = 20.9*TOP+87*sqrt(TOP*T_W); 

display(takeoffDistance) 

%%%%%%%%%%%%%%% LANDING DISTANCE %%%%%%%%%%%%%%%%

C_l_LAmax = 2.7 ; % AVEC RAYMER
% S_L_start = 2800; % ft On choisi au pif !!!

W_S_LA = (W_TO-W_fuel+((reserve_fuel+trapped_fuel)*W_fuel))/S;
LP = W_S_LA*(1/(sigma*C_l_LAmax)); 

landingDistance = 118*LP+400; %Pour vérifier

display(landingDistance)




% Juste pour voir les valeurs

display(rho)
display(V)
display(W2_W3)
display(V_kts)
display(M)
display(C)
display(L_D_max)
display(S)
display(S2)
display(Range_iteration)
display(W_TO)
display(ws_cruising)
