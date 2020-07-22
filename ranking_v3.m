clear all;
clc;
close all;
analisarRank = 1;
analisarFitness = 0;
plotarTodosFitness = 0; %Ou só a media?

ds = datastore("./Dados",'FileExtensions',{'.mat'},'Type','tabulartext');
tamanho = size(ds.Files);
global bons;
bons=0;
qtdd = 0;
k=1;
original_vars = [0 0.9 0 20 20];
global range;
range = 20;
original =  Reflectance_adapted(original_vars);
fileID = fopen('rank.txt','w');
fprintf(fileID,'RANKING MELHORES RESULTADOS\n\n');
fclose(fileID);
warning off;
if(analisarRank)
    for i=1:tamanho(1)
       dado = string(ds.Files(i)); 
       load(dado);
       %fator de qualidade
       if(z(2)<12000)
           bad = 1;
       elseif(size(z)<7)
           bad = 1;
       elseif(z(7)>0.5)
           bad = 1;
       elseif(z(6)>0.1)
           bad = 1;
       else
           bad = 0;
           bons = bons+1;
           rank(bons) = dado;
           
           Teste_Resultado_final_adapted_v3(rank(bons),'rank.txt',original,1);  
           
       end
    end
end
warning on;
if(analisarFitness)
   for i=1:tamanho(1)
       dado = string(ds.Files(i)); 
       load(dado);
       if(exist('parametrosOtimizacao','var')==0)
           bad = 1;
       elseif(i>(tamanho(1)-100))
           if(k==1)
                fitness = zeros(15,parametrosOtimizacao.rodadas); 
           end
           if(size(fitness,2)==parametrosOtimizacao.rodadas)
                fitness(k,:) = parametrosOtimizacao.maiorFITGer(1:parametrosOtimizacao.rodadas);           
           end
           k = k+1;
           
           if (plotarTodosFitness==1)
               if(repeticoes==0)
                   figure();
                   hold on;
               end
               bad = 0;
               plot(1:parametrosOtimizacao.rodadas, parametrosOtimizacao.maiorFITGer(1:parametrosOtimizacao.rodadas));
               title('Evolucao do melhor fitness ao longo das geracoes por repeticao');
               xlabel('Geracao');
               ylabel('Fitness');
               %legend(sprintf('%d',repeticoes));
           end
           if(repeticoes==14)
               legend('0','1','2','3','4','5','6','7','8','9','10','11','12','13','14');
               a = mean(fitness(1:15,:));
               a = transpose(a);
               hold off;
               if(size(a,1) == parametrosOtimizacao.rodadas)
                   figure()
                   title('Media de Evolucao do Melhor Fitness');
                   xlabel('Geracao');
                   ylabel('Fitness');
                   hold on;
                   plot(1:parametrosOtimizacao.rodadas, a);
                   hold off;
               end
               k = 1;
           end
       end
   end
end


