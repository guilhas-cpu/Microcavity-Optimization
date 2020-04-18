clear all;
clc;
close all;
ds = datastore("./Dados",'FileExtensions',{'.jpg','.mat'},'Type','tabulartext');

tamanho = size(ds.Files);

for i=1:tamanho(1)
   dado = string(ds.Files(i)); 
   i = i
   Teste_Resultado_final_adapted_v2(dado,'restore.txt'); 
end


