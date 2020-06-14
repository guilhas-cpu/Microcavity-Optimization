clear all;
clc;
close all;
ds = datastore("./Dados",'FileExtensions',{'.jpg','.mat'},'Type','tabulartext');

tamanho = size(ds.Files);

for i=(tamanho(1)-19):tamanho(1)
   dado = string(ds.Files(i)) 
   i = i
   load(dado)
   s=2;
   save(dado,'vars','z','norm_z','PESO','countSuccesses','countErrors','METODO','maior','range','repeticoes' ,'parametrosOtimizacao','delta','tempoExec','numIter','iter','numRep','dataInicio','dispersion','s')
end


