%--------CONSTANTES-----------

M = 0.41 ;
A = 4 ; % (ft)
H = 29000 ; % (ft)
S = 76.5 ;  % (ft^2)
b = 17.5 ; % ft 
c_tip = 3.32 ; % (ft)   SELON LE DESSIN (DONNÉES HISTORIQUES)
c_root = 5.38 ; % (ft)  SELON LE DESSIN (DONNÉES HISTORIQUES)
Q_aile = 1 ; % p.39 Chapitre 4 (aile basse avec congé de raccordement)
Q_empennage = 1.05 ; % Chapitre 6 (tableau)
effilement = c_tip/c_root ;
c_moy = 2*c_root/3 ; % Corde moyenne
fleche_LE = 30 ; % degrée (SELON LES DONNÉES HISTORIQUE ET MON DESSIN) FLÈCHE !
Alpha_HL_volet = 0.75 ; % pourcentage SELON RAYMER (aucune information dans le cours)   
Alpha_HL_bec = 0.20; % poucentage SELON RAYMER (aucune information dans le cours)   


% Calcul pour les flèches H.L (avec la formule du cours 4)

fleche_HL_volet = atand(tand(fleche_LE)-(Alpha_HL_volet*2*c_root*(1-effilement)/b)) ;
fleche_HL_bec = atand(tand(fleche_LE)-(Alpha_HL_bec*2*c_root*(1-effilement)/b));


W_cruise = 3170 ; %(lbs), Même qu'au début (ce ne devrait pas avoir changé de bcp !) HYPOTHÈSE!!!! (demander au prof ....)

fleche_ht = 12.15 ; % En degrée (SELON LES DONNÉES HISTORIQUE ET MON DESSIN) 
fleche_vt = 23.79 ; % En degrée (SELON LES DONNÉES HISTORIQUE ET MON DESSIN) 

c_t_HT = 2.14; % SELON DONNÉES HISTORIQUES ET MON DESSIN
c_r_HT = 4 ;
c_t_VT = 2.23;
c_r_VT = 5.67;

c_ht = 2*c_r_HT/3 ;
c_vt = 2*c_r_VT/3 ;

eff_ht = c_t_HT/c_r_HT; % effilement H tail
eff_vt = c_t_VT/c_r_VT; % effilement V tail

A_ht = 3.73368 ;
A_vt = 2.59823/2;

D_max = 5 ;
L_tot = 29 ;

Lsection_const = 0 ;
k = 0.09947 ;

l_ht = 12.5 ;
l_vt = 12.5 ;

% Profil (À CHOISIR SE LE PROFIL)

t_c_rel = 0.12 ; % On trouve ça sur airfoil tools
t_c_max = 0.30 ; % Position de l'épaisseur maximal (Pour beaucoup de profils NACA 4- et 5-digits, l'épaisseur max se trouve à 30 % de la corde)
c_prime_c_volet = 1.6 ; % Valeur arbitraire donnée par chapGPT par rapport au gain de surface grâce au volets
c_prime_c_bec = 1.05; % Valeur arbitraire donnée par chapGPT par rapport au gain de surface grâce au slats

Cla = 0.105 ; % voir portance aile
a0 = -2.0;
Clmax = 1.50 ;
dy = 3.4 ;

% ----- C_d0 TOTAL AVION (SHOULD BE 0,0207)-------


% C_do aile

[d_aile,d0_aile,C_d0_aile]= trainee_aile(M, H, fleche_LE, S, A, t_c_rel, t_c_max, effilement, Q_aile, W_cruise);

% C_d0 fuselage ***IL N'Y A PAS DE TRAINÉE D'ONDE, ON EST BCP TROP LENT
% POUR ÇA

[d_fuselage, S_wet_fuselage, C_d0_fuselage ]= trainee_fuselage(M, H, S, D_max, L_tot, Lsection_const, k);

% C_d0 HT et VT (empennage)

[d_HT, S_ht, C_d0_HT] = trainee_emp_horiz(M, H, S, A, effilement, c_ht, l_ht, fleche_ht, t_c_rel, t_c_max, eff_ht, A_ht, Q_empennage) ;
[d_VT, S_vt, C_d0_VT] = trainee_emp_vert(M, H, S, A, effilement, c_vt, l_vt, fleche_vt, t_c_rel, t_c_max, eff_vt, A_vt, Q_empennage);

% --------- C_D0 TOTAL -----------------

C_d0_total = C_d0_VT + C_d0_HT + C_d0_fuselage + C_d0_aile ;


% ----- C_L_max lisse-------


[ CL_a, a_trim, C_l_max_lisse, a_dec ] = portance_aile(M, H, fleche_LE, S, A, t_c_rel, t_c_max, effilement, W_cruise, Cla, a0, Clmax, dy);

% ON DOIT AVOIR C_l,LA = 3,5 et C_l_TO = 1.7


% ------- C_l max VOLET, BEC --------

petit_delta_cl_volet = 1.9 ; %p.15 CHAP 9 TABLEAU (tripplet slotted)
petit_delta_cl_bec = 0.4 ; %p.16 CHAP 9 TABLEAU (Slats)

fraction_s_volet_s_aile = 0.8 ; %entre 60% et 80%

delta_cl_volet = 0.9* petit_delta_cl_volet*c_prime_c_volet* fraction_s_volet_s_aile * cosd(fleche_HL_volet);
delta_cl_bec = 0.9* petit_delta_cl_bec*c_prime_c_bec* fraction_s_volet_s_aile * cosd(fleche_HL_bec);


% C_L_MAX TOTAL

pourcentage = 0.6 ; %entre 0.6 et 0.8 

C_l_max_LA = C_l_max_lisse + delta_cl_volet + delta_cl_bec;
C_l_max_TO = C_l_max_lisse + pourcentage*(delta_cl_volet + delta_cl_bec); 


%----------PARTIE EXTRA (CALCUL INVERSE POUR TROUVER LE BON PROFIL)--------
C_l_max_LA_visee = 3.5;

C_l_lisse_2 = C_l_max_LA_visee - delta_cl_bec - delta_cl_volet;

fleche_quart_corde= atand(tand(fleche_LE)-((1/4)*2*c_root*(1-effilement)/b)) ;

C_l_max_2D = C_l_lisse_2/(0.9*cosd(fleche_quart_corde));


%-----------------DISPLAY DES RÉSULTATS--------------------------

display(delta_cl_bec)
display(delta_cl_volet)
display(C_l_max_LA)
display(C_l_max_TO)
display(C_l_max_lisse)
display(C_d0_total)
display (C_d0_aile)
display(C_d0_fuselage)
display(C_d0_HT)
display(C_d0_VT)
display (fleche_HL_volet)
display(C_l_max_2D)
display(S_ht)
display(S_vt)