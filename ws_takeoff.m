function ws = ws_takeoff(H_takeoff, CL, TW, Sto);
% evaluation de la charge alaire necessaire au decollage
% d'apres Corke, Design of Aircraft, chapitre 3, 2003
%
% Copyright 2008: François Morency
%
% ---------------------------------------------------------------------
% Ces valeurs d'entree de la fonction sont definies par l'usager
% ---------------------------------------------------------------------
H_takeoff=25000; % altitude de decollage en ft
CL=1.7; % coefficient de portance de l'aile au decollage
TW=0.57; % ratio poussee poids au decollage lbf/lbf
Sto=2800; % distance de decollage en ft
% ---------------------------------------------------------------------
% Ces valeurs sont retournees par la fonction
% ---------------------------------------------------------------------
% ws charge alaire decollage: lbf / ft^2

% ---------------------------------------------------------------------
% sigma=densite decollage/densite SL
% ---------------------------------------------------------------------
sigma = density(H_takeoff)/density(0);

% ---------------------------------------------------------------------
% evaluation de TOP en lb/ft^2
% ---------------------------------------------------------------------
f = @(x) Sto - 20.9*x - 87*sqrt(x*TW) ;
TOP=fzero(f,[ 0 1e6 ]);

% ---------------------------------------------------------------------
% evaluation de la charge alaire
% ---------------------------------------------------------------------
f = @(x) TOP - x*1/(CL*TW*sigma);
ws = fzero(f, [ 0 1e3 ]);

sprintf('%s %0.1f', 'Charge alaire necessaire au decollage',ws)
