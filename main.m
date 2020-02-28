clear all;
clc;
close all;
iter2 = 0;
while(iter2<15)
    iter = 0;
    while(iter<1)
        clc;
        close all;
        global range;
        range = 20+12*iter;
        global METODO;
        global maior;
        maior = [0 0 0 0 0 0 0];
        METODO = 1; %metodo 1: soma todos os objetivos, m�todo 2: pega o pior objetivo e otimiza ele
        %metodo 3: otimiza s� o beta
        global print;
        global plota;
        global erro;
        global countErrors;
        global countSuccesses;
        global PESO;
        global inicio;
        inicio = 0;
        print = 0; %#ok<NASGU>
        plota = 0; %#ok<NASGU>
        erro = 0;
        countErrors = 0;
        countSuccesses = 0;
        PESO = [5 20+iter2 1 5 6 6 21]; %favor escolher valores maiores ou iguais a 1. 
        %Quanto maior o peso, maior a relev�ncia do objetivo caso metodo 1 e o contrario caso metodo 2

        inicializacoes();
        global F;
        global n_individuos;
        global chance_mutacao;
        global erro_parada;
        global geracoes_parada;
        global max_geracoes;


        a = EvolucaoDiferencial(F, chance_mutacao, n_individuos, erro_parada, ...
            geracoes_parada, max_geracoes);

        %%%%%%%%%%VALORES DE RESTRICAO DE GERACAO%%%%%%%%%%%%%
        a.adicionarParametro('c01', 'real', [0 1]);
        a.adicionarParametro('c02', 'real', [0 1]);
        a.adicionarParametro('c03', 'real', [0 1]);
        a.adicionarParametro('ncs', 'inteiro', [10 25]);
        a.adicionarParametro('nci', 'inteiro', [10 25]);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        a.adicionarRestricao(@funcao_restricao_area, 'c01', 'c02', 'c03', 'ncs', 'nci');
        a.adicionarFuncaoFitness(@funcao_fitness);
        a.adicionarFuncaoObjetivo(@funcao_objetivo,'y');

        a.gerarPopulacaoInicial();

        tic()
        % Gera mutacoes NUM_RODADAS vezes
        while ~a.finalizou
            a.gerarNovaPopulacao();
        end

        [FO, individuo] = a.mostrarResultados();

        print = 1;
        plota = 1;
        funcao_objetivo(a.populacao(individuo).parametros);

        fprintf("\nTempo de execucao: %fs \n", toc());

        figure()
        plot(1:a.rodadas, a.maiorFITGer(1:a.rodadas));
        title('Evolucao do fitness ao longo das gera��es');
        xlabel('Gera��o');
        ylabel('Fitness')


        iter = iter + 1;
    end
    iter2 = iter2 +1;
end

 function y = funcao_fitness(vars) 
        y = vars(1);
    end


    function y = funcao_restricao_area(vars)
        global countErrors;
        if((vars(1) >= 0 && vars(1) <= 1) || (vars(2) >= 0 && vars(2) <= 1) || (vars(3) >= 0 && vars(3) <= 1) || (vars(4) >= 1 && vars(1) <= 25) || (vars(5) >= 1 && vars(5) <= 25))
            y = true;
        else
            y = false;
        end
        if(y==false)
            countErrors = countErrors + 1;
            if(mod(countErrors,10000)==0)
                countErrors = countErrors
            end
        end
    end

    function y = funcao_objetivo(vars)
        global erro;
        global countErrors;
        global countSuccesses;
        global PESO;
        global print;
        global METODO;
        global maior;
        global inicio;
        global iter;
        global range;
        inicio = inicio + 1;
        norm_z = zeros(1,7);

        z = Reflectance_adapted(vars);

        %%%%%%%%%%NORMALIZACAO%%%%%%%%%%
        %Deve ficar de 0 a 1, quanto maior, melhor
        %Valores base de normaliza��o de 3 vezes o valor encontrado no codigo enviado pelo
        %cotta
        norm_z(1) = z(1)/(3*0.0118374112);                                   %beta - MAXIMIZE
        norm_z(2) = z(2)/(3*5.281276e+03);                            %Q - MAXIMIZE
        norm_z(3) = z(3)/(3*1.718943e-01);                                   %Fp - MAXIMIZE
        norm_z(4) = ((3*(z(4)+1e10))/(1e10+9.397523e+08))^-1;           %DN -  MINIMIZE
        norm_z(5) = ((3*(z(5)+2e-07))/(2e-07+9.546695e-08))^-1;        %Cav_loss - MINIMIZE
        norm_z(6) = ((3*(z(6)+1))/(1+9.381136e-01))^-1;                  %Dist - MINIMIZE
        norm_z(7) = ((3*(z(7)+1))/(1+6.076570e-01))^-1;                            %Reflec - MINIMIZE
        norm_z_com_PESO = norm_z.*(PESO);

        for i=1:7
            if(norm_z(i)>maior(i))
                maior(i)=norm_z(i) %armazena o melhor resultado obtido pra cada objetivo at� ent�o
                %atencao, n�o se trata do melhor obtido simultaneamente entre
                %todos os objetivos, mas individualmente. Mostra "o melhor" de
                %cada um
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%PEGAR PIOR CASO%%%%%%%%
        if(METODO==1)
            valor = sum(norm_z_com_PESO);
        elseif(METODO==2)    
            [valor indice] = min(norm_z_com_PESO);
        elseif(METODO==3)
            valor = norm_z_com_PESO(1);
        end

        y = [valor erro];


        if(erro)
            %CASO INDIVIDUO GERADO SEJA INVALIDO, erro vira 1
            if(inicio~=1)    
                erro = 0;
                valor = 0;
                y = [valor erro]; %a partir da segunda geracao, d� valor ruim pro fitness mas n�o informa mais como erro e n�o descarta
            end
            countErrors = countErrors+1 %sem ponto e virgula vai imprimir na tela toda hora
            %apesar do individuo n�o ser descartado a partir da segunda rodada, ainda vai
            %informar que deu erro pra gente ter uma no��o boa
            date = clock;
            errorName = sprintf('./Erros/%d_%d_%d-%dh%dmin.mat',date(1),date(2),date(3),date(4),date(5));
            save(errorName,'vars')
        else
            norm_z = norm_z %informa normas obtidas para individuos validos
            countSuccesses = countSuccesses+1 %sem ponto e virgula vai imprimir na tela toda hora
        end

        if(print&&countSuccesses>50) %Se tiver menos que 50 sucessos, resultado � ruim
            date = clock;
            matFileName = sprintf('./Dados/%d_%d_%d-%dh%dmin.mat',date(1),date(2),date(3),date(4),date(5));
            save(matFileName,'vars','z','norm_z','PESO','countSuccesses','countErrors','METODO','maior','range','iter')
        end
    end