
mylist = [] ;


for Hauteur = 10000:5000:35000 

    rho = density(Hauteur) ;
    
    v_initial = 245*1.687811; %en ft/sec pour le calcul de W/s cruise pour avoir le vrai V
    M = v_initial / sqrt(1.4 * 1716 * tempatmstd(Hauteur)*1.8); % Calculate Mach number
    
    
    C_d0 = 0.0207 ; %Le range change BCP quand je change le C_d0
    e = 0.8 ;
    A = 4;
    k = 1/(pi*A*e) ;
    
    reserve_fuel = 0.05 ;
    trapped_fuel= 0.01;
    W_payload = 500 + 2*180 ; %en lbs
    
    ws_cruising = (rho*(v_initial^2)/2)*sqrt(C_d0/k);
    V = sqrt((2/rho)*ws_cruising)*((k/C_d0)^0.25); %ft/s
    V_kts = V/1.687811; % Vitesse en kts
    
    
    v_initial = 0 ;
    V_kts = 245 ; % vitesse en KTS
    
    while (abs(v_initial-V_kts) > 0.1)
    
        v_initial = V_kts ;
    
        v_initial_1 = v_initial*1.687811; %en ft/sec pour le calcul de W/s cruise pour avoir le vrai V
        M = v_initial_1 / sqrt(1.4 * 1716 * tempatmstd(Hauteur)*1.8); % Calculate Mach number
    
        ws_cruising = (rho*(v_initial_1^2)/2)*sqrt(C_d0/k);
    
        V = sqrt((2/rho)*ws_cruising)*((k/C_d0)^0.25); %ft/s
        V_kts = V/1.687811; % Vitesse en kts
    
        display(V_kts)
    
    end
    
    L_D_max = sqrt(1/(4*C_d0*k));
    C_BHP=0.5;
    rendement = 0.7 ;
    C = consommation_spec(C_BHP,V,rendement);
    
    
    % À CHANGER À CHANQUE ITÉRATION (faire une boucle)
    Range_start= 0; %nm commencer à 600nm
    Range_iteration = 600 ;
    
    while (abs(Range_start-Range_iteration) > 0.1)
    
        Range_start = Range_iteration ;
    
        [W_TO,W_fuel,W_empty] = itertow('Aviation-gen(1mot)',M,Hauteur,A,C,0,reserve_fuel,trapped_fuel,W_payload,Range_start); %Aviation-gen(1mot) Je suis pas sûre que ça soit ça (on a prop 1 mot, mais il donne slm prop 2moteurs)
        
        %CALCULS POUR ITÉRATION SUR R (on recalcule R avec W_to et W_fuel)
        
        Wla_Wto = 1- ((W_fuel/W_TO) / (1 + reserve_fuel + trapped_fuel));
        W1_Wto = 0.975 ;
        W2_W1 = 0.975;
        Wla_W_3 = 0.995;
        W2_W3 = (W1_Wto*W2_W1*Wla_W_3)/(Wla_Wto) ;
    
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
    
    % s_TO_start = 2800 ; % ft, on choisi au pif !!!!!
    
    C_l_TOmax = 1.7 ; % AVEC RAYMER
    T_W = 0.57 ; %Donnée historique (0,25-0,57)
    W_S_TO = ws_cruising;
    
    TOP = W_S_TO* (1/C_l_TOmax)*(1/T_W)*(1/sigma) ; 
    
    takeoffDistance = 20.9*TOP+87*sqrt(TOP*T_W); 
    
    
    
    %%%%%%%%%%%%%%% LANDING DISTANCE %%%%%%%%%%%%%%%%
    
    C_l_LAmax = 3.5 ; % AVEC RAYMER
    % S_L_start = 3000; % ft On choisi au pif !!!
    
    W_S_LA = (W_TO-W_fuel+((reserve_fuel+trapped_fuel)*W_fuel))/S;
    LP = W_S_LA*(1/(sigma*C_l_LAmax)); 
    
    landingDistance = 118*LP+400; %Pour vérifier
    
    V_s = sqrt(ws_cruising*2/(rho*C_l_LAmax));
    Puissance_disponible = T_W * V_s*1.1*0.7 ;
    P_hp = Puissance_disponible/(550*rendement*sigma);

    mylist(end+1) = W_TO;
end

Hauteur = 10000:5000:35000 ; %On peux le modifier de 25000 ft à 29000 ft

figure
plot(Hauteur, mylist, '-o','LineWidth',1.5)
grid on
xlabel('Altitude (ft)')
ylabel('Takeoff Weight W_{TO} (lbf)')
title('W_{TO} vs Altitude')


%%%%%% GRAPHIQUE %%%%%%
