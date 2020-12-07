clear all;
clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%ANALISES%%%%%%%%%%%%%%%%%%%%%%%%
analisarRank = 1;
gerarTabela = 1; %tabela com melhores resultados listados, com melhor, pior e fitness mediano de cada melhor resultado
analisarFitness = 0;
plotarTodosFitness = 0; %Ou só a media?
numRepMedia = 5; %Numero de Repeticoes com os mesmos parametros para fazer media
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%VALORES DE FILTRO%%%%%%%%%%%%%%%%%%%
Quality_Factor = 12500; %12000
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
           if(gerarTabela) %adds data from each best result
                PESO_base = [5 21 1 5 6 15 25]; 
                tablename = 'rankingTable.xlsx';
                [~,tasks] = system('tasklist/fi "imagename eq Excel.exe"');
                if contains(tasks, 'nenhuma tarefa')==0; system('taskkill /F /IM EXCEL.EXE'); end
                parametros = cell2mat(struct2cell(parametrosOtimizacao.populacao));
                [valor,HE_index] = max(parametros(1,1,:)); %maior fitness
                [valor,LE_index] = min(parametros(1,1,:)); %menor fitness
                [valor,AE_index] = min(abs(parametros(1,1,:)-median(parametros(1,1,:)))); %fitness mediano
                
                if(size(parametros,1)==9) %tem 6 vars
                    vars_HE = parametros(4:9,1,HE_index);
                    vars_LE = parametros(4:9,1,LE_index);
                    vars_AE = parametros(4:9,1,AE_index);
                else
                    vars_HE = parametros(4:8,1,HE_index);
                    vars_LE = parametros(4:8,1,LE_index);
                    vars_AE = parametros(4:8,1,AE_index);
                end
                z_HE = Reflectance_final_results(vars_HE,1);
                z_LE = Reflectance_final_results(vars_LE,1);
                z_AE = Reflectance_final_results(vars_AE,1);
                
                %normalizing
                result_HE = [z_HE(1)/(3*0.0118374112) z_HE(2)/(3*5.281276e+03) z_HE(3)/(3*1.718943e-01) ((3*(z_HE(4)+1e10))/(1e10+9.397523e+08))^-1 ((3*(z_HE(5)+2e-07))/(2e-07+9.546695e-08))^-1 ((3*(z_HE(6)+1))/(1+9.381136e-01))^-1 ((3*(z_HE(7)+1))/(1+6.076570e-01))^-1];
                result_LE = [z_LE(1)/(3*0.0118374112) z_LE(2)/(3*5.281276e+03) z_LE(3)/(3*1.718943e-01) ((3*(z_LE(4)+1e10))/(1e10+9.397523e+08))^-1 ((3*(z_LE(5)+2e-07))/(2e-07+9.546695e-08))^-1 ((3*(z_LE(6)+1))/(1+9.381136e-01))^-1 ((3*(z_LE(7)+1))/(1+6.076570e-01))^-1];
                result_AE = [z_AE(1)/(3*0.0118374112) z_AE(2)/(3*5.281276e+03) z_AE(3)/(3*1.718943e-01) ((3*(z_AE(4)+1e10))/(1e10+9.397523e+08))^-1 ((3*(z_AE(5)+2e-07))/(2e-07+9.546695e-08))^-1 ((3*(z_AE(6)+1))/(1+9.381136e-01))^-1 ((3*(z_AE(7)+1))/(1+6.076570e-01))^-1];                          

                Fi_AE(bons,1) = sum((PESO_base).*result_AE);
                Fi_HE(bons,1) = sum((PESO_base).*result_HE);
                Fi_LE(bons,1) = sum((PESO_base).*result_LE);
                Q_AE(bons,1) = z_AE(2);
                Q_HE(bons,1) = z_HE(2);
                Q_LE(bons,1) = z_LE(2);
                dist_AE(bons,1) = z_AE(6);
                dist_HE(bons,1) = z_HE(6);
                dist_LE(bons,1) = z_LE(6);
                Best_Result(bons,1)=bons;

           end
       end
    end
   
    if(gerarTabela) %Generate table on Excel
    	T = table(Best_Result,Q_HE,dist_HE,Fi_HE,Q_AE,dist_AE,Fi_AE,Q_LE,dist_LE,Fi_LE);
        if exist(tablename, 'file') ; delete(tablename); end
    	writetable(T,tablename,'Sheet',1,'Range','D1');
    	writetable(T,tablename,'Sheet','MyNewSheet','WriteVariableNames',false);
    end
    
end


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

warning on;


