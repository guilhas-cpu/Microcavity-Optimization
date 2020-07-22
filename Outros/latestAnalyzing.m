clear all;
clc;
close all;
ds = datastore("./Dados",'FileExtensions',{'.mat'},'Type','tabulartext');

tamanho = size(ds.Files);
ultimos=0;

for i=1:tamanho(1)
   dado = string(ds.Files(i)); 
   if(i<tamanho(1)-10)
       bad = 1;
   else
       ultimos = ultimos+1;
       latest(ultimos) = dado;
   end
end


fileID = fopen('latest.txt','w');
fprintf(fileID,'ULTIMOS RESULTADOS\n\n');
fclose(fileID);

for i=1:ultimos
   Teste_Resultado_final_adapted_v2(latest(i),'latest.txt');    
end