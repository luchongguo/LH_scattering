%% calculate drive efficiency
%% ne(10^19m^-3),P(MW),I(kA)
function [effcd]=calc_cd(ne,P,I)
R=1.91;
effcd=ne*I*R/P*1e-3;
