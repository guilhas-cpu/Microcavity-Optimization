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
   if(exist('parametrosOtimizacao','var')==0)
       bad = 1;
   else
       if(repeticoes==0)
           hold off;
           figure();
           hold on;
       end
       bad = 0;
       bons = bons+1;
       rank(bons) = dado;
       plot(1:parametrosOtimizacao.rodadas, parametrosOtimizacao.maiorFITGer(1:parametrosOtimizacao.rodadas));
       title('Evolucao do melhor fitness ao longo das geracoes');
       xlabel('Geracao');
       ylabel('Fitness');
       %legend(sprintf('%d',repeticoes));
       if(repeticoes==14)
            legend('0','1','2','3','4','5','6','7','8','9','10','11','12','13','14');
       end
   end
end
hold off;

fileID = fopen('rank_v2.txt','w');
fprintf(fileID,'RANKING MELHORES RESULTADOS\n\n');
fclose(fileID);

for i=1:bons
   Teste_Resultado_final_adapted(rank(i),'rank_v2.txt');    
end