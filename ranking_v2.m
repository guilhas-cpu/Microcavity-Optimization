clear all;
clc;
close all;
ds = datastore("./Dados",'FileExtensions',{'.mat'},'Type','tabulartext');

tamanho = size(ds.Files);
bons=0;
qtdd = 0;
k=1;
for i=1:tamanho(1)
   dado = string(ds.Files(i)); 
   load(dado);
   %fator de qualidade
   if(exist('parametrosOtimizacao','var')==0)
       bad = 1;
   else
        if(k==1)
            fitness = zeros(15,parametrosOtimizacao.rodadas); 
        end
       if(size(fitness,2)==parametrosOtimizacao.rodadas)
            fitness(k,:) = parametrosOtimizacao.maiorFITGer(1:parametrosOtimizacao.rodadas);           
       end
       k = k+1;
       if(repeticoes==0)
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
           hold off;
           figure()
           title('Fitness médio');
           hold on;
           a = mean(fitness(1:15,:));
           a = transpose(a);
           plot(1:parametrosOtimizacao.rodadas, a);
           hold off;
           k = 1;
       end
   end
end
hold off;

%fileID = fopen('rank_v2.txt','w');
%fprintf(fileID,'RANKING MELHORES RESULTADOS\n\n');
%fclose(fileID);

%for i=1:bons
   %Teste_Resultado_final_adapted_v2(rank(i),'rank_v2.txt');    
%end