clear all;
clc;
close all;

global print;
global plota;
print = 1;
plota = 1;

%Dê load no arquivo que tem os resultados
load('2020_2_16-20h57min.mat')
%load('2019_12_27-18h23min.mat')
%load('2020_1_13-23h6min.mat')
%load('AAAA_M_DD-HHhMMmin.mat')

result = Reflectance_adapted(vars); %registra os resultados dos objetivos
original = [0.0118374112 3.6904333333337e+03 0.1201348351 1.34463927383e+09 9.546694853150e-08 0.89104208416836 0.7139]; 
melhora = [result(1)/original(1) result(2)/original(2) result(3)/original(3) original(4)/result(4) original(5)/result(5) original(6)/result(6) original(7)/result(7)] - 1; %registra a melhora em relação aos resultados base (enviado pelo professor cotta)
melhora = melhora*100;


fprintf('\n\nRESULTADOS BASE\t\t\t\t\t\tNOVOS RESULTADOS\t\t\t\t\tPERCENTUAL DE MELHORA\n');
fprintf('beta = 1.183741e-02\t\t\t\t\tbeta = %d \t\t\t\t(%f%%)\n',result(1),melhora(1));
fprintf('Q = 3.690433e+03\t\t\t\t\tQ = %d \t\t\t\t\t(%f%%)\n',result(2),melhora(2));
fprintf('Fp = 1.201348e-01\t\t\t\t\tFp = %d \t\t\t\t\t(%f%%)\n',result(3),melhora(3));
fprintf('DN = 1.344639e+09\t\t\t\t\tDN = %d \t\t\t\t\t(%f%%)\n',result(4),melhora(4));
fprintf('CavityLoss = 9.546695e-08\t\t\tCavityLoss = %d \t\t\t(%f%%)\n',result(5),melhora(5));
fprintf('Distance = 8.910421e-01\t\t\t\tDistance = %d \t\t\t(%f%%)\n',result(6),melhora(6));
fprintf('Reflec_min= 7.139e-01\t\t\t\tReflec_min = %d \t\t\t(%f%%)\n',result(7),melhora(7)); 
