clear all;
clc;
close all;
ds = datastore("./Dados",'FileExtensions',{'.mat'},'Type','tabulartext');

tamanho = size(ds.Files);
bons=0;
qtdd = 0;

for i=1:tamanho(1)
   dado = string(ds.Files(i)); 
   load(dado);
   %fator de qualidade
   if(z(2)<10000)
       bad = 1;
   elseif(size(z)<7)
       bad = 1;
   elseif(z(7)>0.7139)
       bad = 1;
   elseif(z(6)>0.1)
       bad = 1;
   else
       bad = 0;
       bons = bons+1;
       vars = vars
       rank(bons) = dado;
   end
end


fileID = fopen('rank.txt','w');
fprintf(fileID,'RANKING MELHORES RESULTADOS\n\n');
fclose(fileID);

original_vars = [0 0.9 0 20 20];
original =  Reflectance_adapted(original_vars);

for i=1:bons
   Teste_Resultado_final_adapted_v2(rank(i),'rank.txt',original);    
end