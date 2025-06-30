function [TSFC] = consommation_spec(CBHP,Vit,eff)
%UNTITLED Convert CBHP to TSFC
%   Output variable
%TSFC = lb/hr/lb
%   Input Variables
% CBHP = lb/hr/puissance au frein
% Vit = Vitesse de vol (ft/s)
% eff = rendement de l'hélice
%CBHP=0.45;
%Vit = 100 ft/s;
%eff = 0.75;
%
% conversion
TSFC=CBHP*Vit/(550*eff);

end

