% Gera a populacao inicial
% Aqui pegamos um ponto de referencia (candidato_inicial) e adicionamos um
% valor aleatorio dentro de limites para gerar populacoes iniciais
function gerarPopulacaoInicial(obj)
    global numIndiv;
    global erro;

    % Inicializa a matriz com informacao de candidatos a solucao
    obj.populacao(obj.num_individuos).FO = {};
    obj.populacao(obj.num_individuos).FITNESS = {};
    obj.populacao(obj.num_individuos).erro = {};
    obj.populacao(obj.num_individuos).parametros= zeros(obj.num_parametros,1);
    
    % Nao houveram rodadas de execucao do algoritmo
    obj.rodadas = 0;
    obj.finalizou = false;
    
    i = 1;
    while i~=obj.num_individuos+1
        numIndiv = i;
        for j=1:obj.num_parametros
            % Limites inferior e superior do parametro
            limi = obj.parametros(j).liminferior;
            lims = obj.parametros(j).limsuperior;
            
            % tipo inteiro
            if obj.parametros(j).tipo == 0
                obj.populacao(i).parametros(j) = randi([limi lims]);
                
            % tipo real
            else
                % Aqui geramos um novo ponto, com um numero aleatorio
                % acrescentado rand() que e uma funcao que varia de 0 a 1,
                % e da forma que fizemos, varia de lims a limi
                obj.populacao(i).parametros(j) = ...
                ((lims-limi)*rand()+limi);
            end
            
        end
        
        % Verifica se o individuo esta dentro do limite de restricao, e se
        % nao tiver, refaz
        if(obj.verificarRestricoes(obj.populacao(i).parametros))
            
            % Armazena o resultado da funcao objetivo
            FOeErro = obj.calcularFuncaoObjetivo(obj.populacao(i).parametros(:));
            
            % Se nao teve erro, armazena e continua
            if(FOeErro(obj.funcao_objetivo.num_saidas+1) == 0)
                FO = FOeErro(1:obj.funcao_objetivo.num_saidas);
                obj.populacao(i).FO = FO;
                obj.populacao(i).FITNESS = obj.calcularFitness(FO);
                % Indica que o erro em relacao ao resultado da geracao
                % anterior e infinito
                obj.populacao(i).erro = inf;
                i = i+1;
            end
        end
        
    end
    
    % Obtem o maior FITNESS
    vetorFIT = zeros(obj.num_individuos, 1);
    for i=1:obj.num_individuos
        vetorFIT(i) = obj.populacao(i).FITNESS;
    end

    obj.maiorFIT = max(vetorFIT);
    obj.maiorFITGer(1) = obj.maiorFIT;
    
    obj.rodadas = 1;
    
    %DEBUG
    %obj.mostrarResultados();
end