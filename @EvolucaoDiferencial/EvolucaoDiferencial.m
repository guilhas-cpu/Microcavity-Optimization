classdef EvolucaoDiferencial < handle
    
    %% Variaveis da classe
    properties(Access = public)
        
        % Armazena a funcao objetivo
        funcao_objetivo = struct(),
        
        % Armazena a funcao fitness
        funcao_fitness,
        
        % Armazena as funcoes de restricao
        funcao_restricao = struct(), num_restricoes = 0,
        
        % Lista contendo informacoes dos nomes dos parametros
        parametros = struct(), num_parametros = 0,
        
        % Matriz contendo os individuos de nossa populacao
        populacao = struct(), num_individuos,
        
        % Maior FITNESS encontrado
        maiorFIT, maiorFITGer,
        
        % Peso para calculo do vetor ruido
        F,
        
        % Chance de mutacao
        chance_mutacao,
        
        % Quantifica o erro desejado para realizar a parada do algoritmo
        erro_parada, quant_erro,
        
        % Indica por quantas geracos o erro se deve perpetuar para que pare
        % o algoritmo
        quant_geracos_erro,
        
        % Indica quantas rodadas o algoritmo executou
        rodadas,
        
        % Numero maximo de geracoes do algoritmo
        max_geracoes,
        
        % Indica que o algoritmo finalizou
        finalizou = false
        
    end
    
    %% Metodos publicos da classe (podem ser acessados por qualquer um)
    methods(Access = public)
        
        %% Metodos da classe
        adicionarParametro(obj, nome_parametro, tipo, limites),
        adicionarRestricao(obj, funcao_restricao, varargin),
        adicionarFuncaoObjetivo(obj, funcao_objetivo, varargin),
        adicionarFuncaoFitness(obj, funcao_fitness),
        gerarPopulacaoInicial(obj),
        gerarNovaPopulacao(obj),
        [FO, individuo] = mostrarResultados(obj)
        
        % Construtor da classe
        function obj = EvolucaoDiferencial(F, ...
                chance_mutacao, num_individuos, ...
                erro_parada, quant_geracos_erro, max_geracoes)
            
            obj.F = F;
            obj.chance_mutacao = chance_mutacao;
            obj.num_individuos = num_individuos;
            obj.erro_parada = erro_parada;
            obj.quant_geracos_erro = quant_geracos_erro;
            obj.max_geracoes = max_geracoes;
            obj.maiorFITGer = zeros(max_geracoes, 1);
        end
        
    end
    
    %% Metodos privados da classe (podem ser acessados somente pela propria classe)
    methods(Access = private)
        
        vetor_ruido = gerarMutacao(obj, vetor_base, indice_vetor_base),
        vetor = gerarCruzamento(obj, parametros_individuo, vetor_ruido),
        y = verificarRestricoes(obj, individuo),
        FO = calcularFuncaoObjetivo(obj, parametros),
        FITNESS = calcularFitness(obj, FO)
        
    end
end