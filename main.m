clear all;
clc;
close all;

global range delta maior METODO print plota saveErrors countGen repeticoes;
global erro countErrors countSuccesses PESO geracao numIndiv;
global tempoExec F n_individuos chance_mutacao erro_parada geracoes_parada max_geracoes parametrosOtimizacao;
global iter numIter numRep dataInicio dispersion isLambdaR1;
isLambdaR1=0;
iter = 0;
numIter = 2;
numRep = 30;
inicio = clock;
dataInicio = sprintf('%.2d/%.2d/%.2d-%.2dh%.2dmin',inicio(3),inicio(2),inicio(1),inicio(4),inicio(5));

while(iter<numIter)
    repeticoes = 0;
    while(repeticoes<numRep)
        clc;
        close all;
        delta = 0.25; %25%
        range = 20;
        maior = [0 0 0 0 0 0 0];
        METODO = 1; %metodo 1: soma todos os objetivos, método 2: pega o pior objetivo e otimiza ele
        %metodo 3: otimiza só o beta
        inicio = 0;
        print = 0; %#ok<NASGU>
        plota = 0; %#ok<NASGU> 
        erro = 0;
        geracao = 1;
        numIndiv = 1;
        countErrors = 0;
        countSuccesses = 0;
        countGen = 0;
        saveErrors = 0;
        PESO = [5 21 1 5 6 15 25]; %favor escolher valores maiores ou iguais a 1.

        
        %Quanto maior o peso, maior a relevância do objetivo caso metodo 1 e o contrario caso metodo 2

        inicializacoes();
        a = EvolucaoDiferencial(F, chance_mutacao, n_individuos, erro_parada, ...
            geracoes_parada, max_geracoes);

        %%%%%%%%%%VALORES DE RESTRICAO DE GERACAO%%%%%%%%%%%%%
        a.adicionarParametro('c01', 'real', [0 1]);
        a.adicionarParametro('c02', 'real', [0 1]);
        a.adicionarParametro('c03', 'real', [0 1]);
        switch iter
            case 0
                a.adicionarParametro('ncs', 'inteiro', [10 25]);
                a.adicionarParametro('nci', 'inteiro', [10 25]);
            case 1
                a.adicionarParametro('ncs', 'inteiro', [10 20]);
                a.adicionarParametro('nci', 'inteiro', [10 20]);
        end
        a.adicionarParametro('shift', 'real', [0 20]);
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
        parametrosOtimizacao = a;
        funcao_objetivo(a.populacao(individuo).parametros);

        fprintf("\nTempo de execucao: %fs \n", toc());
        tempoExec = toc();


        figure()
        plot(1:a.rodadas, a.maiorFITGer(1:a.rodadas));
        title('Evolucao do melhor fitness ao longo das geracoes');
        xlabel('Geracao');
        ylabel('Fitness')


        repeticoes = repeticoes +1;
    end
    iter = iter+1
end

 function y = funcao_fitness(vars) 
        y = vars(1);
    end


    function y = funcao_restricao_area(vars)
        global countGen;
        countGen = countGen + 1;
        global delta; %delta da diferença entre c01 e c02, tem que ser maior que ele a diferença
        if((vars(1) >= 0 && vars(1) <= 1) && (vars(2) >= 0 && vars(2) <= 1) && (abs(vars(2)-vars(1)) >= delta) && (vars(3) >= 0 && vars(3) <= 1) && (vars(4) >= 1 && vars(1) <= 25) && (vars(5) >= 1 && vars(5) <= 25))
            y = true;
        else
            y = false;
        end
    end

    function y = funcao_objetivo(vars)
        global erro repeticoes saveErrors countErrors countSuccesses PESO print METODO maior geracao numIndiv;
        global range countGen delta parametrosOtimizacao tempoExec isLambdaR1; 
        global iter numIter numRep n_individuos max_geracoes dataInicio deltaTempo tempoPrevisto;
        dispersion =1;
        s = 1; %default
        norm_z = zeros(1,7);
        z = Reflectance_adapted(vars);
        %%%%%%%%ANALISE DE TEMPO%%%%%%%%%%%%
        if(numIndiv==1&&geracao>1)
            deltaTempo = 0;
            tic();
        elseif(numIndiv==2&&geracao>1)
            deltaTempo = toc();
            temposIterCompletas = (deltaTempo/60)*n_individuos*max_geracoes*numRep*(numIter-(iter+1));
            temposRepCompletas = (deltaTempo/60)*n_individuos*max_geracoes*(numRep-(repeticoes+1))*1;
            temposGeracoesCompletas = (deltaTempo/60)*n_individuos*(max_geracoes-(geracao))*1*1;
            temposIndividuosCompletos = (deltaTempo/60)*(n_individuos-(numIndiv))*1*1*1;
            aux = temposIterCompletas + temposRepCompletas + temposGeracoesCompletas + temposIndividuosCompletos;
            dias = floor(aux/(60*24));
            horas = floor(aux/60-24*dias);
            minutos = ceil(aux-24*60*dias-60*horas);
            tempoPrevisto = sprintf("%.2d dias, %.2d horas e %.0f minutos",dias,horas,minutos);
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%LIMITA�OES REALISTICAS%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%0.2 nm j� � limite real�stico%%%%%%%%%%%%%%%%
        if(z(6)<0.1)
            z(6) = 0.1;
        end        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%NORMALIZACAO%%%%%%%%%%
        %Deve ficar de 0 a 1, quanto maior, melhor
        %Valores base de normalização de 3 vezes o valor encontrado no codigo enviado pelo
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
                maior(i)=norm_z(i); %armazena o melhor resultado obtido pra cada objetivo até então
                %atencao, não se trata do melhor obtido simultaneamente entre
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


        if(erro)%||z(2)<4000||z(7)>0.7139||z(6)>0.1)
            %CASO INDIVIDUO GERADO SEJA INVALIDO, erro vira 1
            erro = 0;
            if(geracao>1)    
                valor = 0;
                y = [valor erro]; %a partir da segunda geracao, dá valor ruim pro fitness mas não informa mais como erro e não descarta
            end
            countErrors = countErrors+1; %sem ponto e virgula vai imprimir na tela toda hora
            %apesar do individuo não ser descartado a partir da segunda
            %rodada, ainda vai salvar como erro ocorrido
            if(saveErrors)
                %informar que deu erro pra gente ter uma noção boa
                date = clock;
                errorName = sprintf('./Erros/%d_%d_%d-%dh%dmin%ds.mat',date(1),date(2),date(3),date(4),date(5),ceil(date(6)));
                save(errorName,'vars')
            end
        else
            norm_z = norm_z; %informa normas obtidas para individuos validos
            countSuccesses = countSuccesses+1; %sem ponto e virgula vai imprimir na tela toda hora
        end

        if(print&&countSuccesses>1) %Se tiver menos que 50 sucessos, resultado eh ruim
            date = clock;
            matFileName = sprintf('./Dados/%.2d_%.2d_%.2d-%.2dh%.2dmin.mat',date(1),date(2),date(3),date(4),date(5));
            save(matFileName,'vars','z','norm_z','PESO','countSuccesses','countErrors','METODO','maior','range','repeticoes' ,'parametrosOtimizacao','delta','tempoExec','numIter','iter','numRep','dataInicio','dispersion','s','isLambdaR1')
        end
        
        
        
        %update de informacao dos valores de interesse
        fprintf('\n....................................................................................................');
        fprintf('\nRESULTADOS PARCIAIS:\nIter:%d/%d\t\t\tRepeticao %d/%d\t\t\tGeracao %d/%d\t\t\tIndividuo %d/%d',iter+1,numIter,repeticoes+1,numRep,geracao,max_geracoes,numIndiv,n_individuos);
        fprintf('\nAcertos: %d\t\t\tErros: %d\t\t\t', countSuccesses, countErrors);
        fprintf('\nCome�ou: %s\t\t\tTempo RESTANTE Previsto: %s', dataInicio, tempoPrevisto);
        %fprintf('\t\t\tGerados: %d', countGen);
        %fprintf(' (delta = %.1f%%).', delta*100);
        norm_z = norm_z
        maior = maior
    end