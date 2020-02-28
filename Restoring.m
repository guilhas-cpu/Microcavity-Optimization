clear all;
clc;
close all;
ds = datastore("./Dados",'FileExtensions',{'.jpg','.mat'},'Type','tabulartext');

tamanho = size(ds.Files);
bons=0;
qtdd = 0;

for i=1:tamanho(1)
   dado = string(ds.Files(i)); 
   Teste_Resultado_final_adapted(dado);   %ultima linha desse codigo não comentada
end


