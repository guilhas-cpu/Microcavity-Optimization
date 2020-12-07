clear all;
clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%ANALISES%%%%%%%%%%%%%%%%%%%%%%%%
analisarRank = 1;
analisarFitness = 1;
plotarTodosFitness = 1; %Ou só a media?
numRepMedia = 5; %Numero de Repeticoes com os mesmos parametros para fazer media
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%VALORES DE FILTRO%%%%%%%%%%%%%%%%%%%
Quality_Factor = 12000; %12000
Distance_Peaks = 0.1; %0.1
Minimum_Reflectance = 0.5; %0.5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Data_Import = fileDatastore("./Dados",'ReadFcn',@load,'FileExtensions','.mat');
Data_Read = readall(Data_Import);
Data_Size = size(Data_Read,1);
global bons range;
bons=0;
range = 20;
qtdd = 0;
k=1;
original_vars = [0 0.9 0 20 20];
original =  Reflectance_adapted(original_vars);

fileID = fopen('rank.txt','w');
fprintf(fileID,'RANKING MELHORES RESULTADOS\n\n');
fclose(fileID);

warning off;
if(analisarRank)
    for i=1:Data_Size(1)
       Data_String = string(Data_Import.Files(i)); 
       load(Data_String);
       %fator de qualidade
       if(z(2) < Quality_Factor)
           bad = 1;           
       elseif(size(z) < 7)
           bad = 1;
       elseif(z(6) > Distance_Peaks)
           bad = 1;
       elseif(z(7) > Minimum_Reflectance)
           bad = 1;
       else
           bad = 0;
           bons = bons+1;
           rank(bons) = Data_String;
           Teste_Resultado_final_adapted_v3(rank(bons),'rank.txt',original,1);  
       end
    end
end

warning on;
if(analisarFitness)
   for i=1:Data_Size(1)
       Data_String = string(Data_Import.Files(i)); 
       load(Data_String);
       if(exist('parametrosOtimizacao','var')==0)
           bad = 1;
       elseif(i>(Data_Size(1)-100))
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
           if(repeticoes==numRepMedia-1)
               legend('0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29');
               size_a = size(fitness,1);
               a = mean(fitness(1:size_a,:));
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


