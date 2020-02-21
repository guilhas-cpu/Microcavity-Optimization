function gerarNovaPopulacao(obj)

    i = 1;
    while i~=obj.num_individuos+1
        
        individuo = obj.populacao(i).parametros;
        vetor_ruido = gerarMutacao(obj, individuo, i);
        cruzamento = obj.gerarCruzamento(individuo, vetor_ruido);
        
        % Verifica se foi atendido as restricos, e se foi possivel fazer a selecao,
        % Se nao, refaz os calculos
        if(obj.verificarRestricoes(cruzamento) && ...
            obj.fazerSelecao(i, cruzamento))
                
                i = i+1;
        end
        
    end
    
    % Atualiza o maior fitness
    vetorFIT = zeros(obj.num_individuos, 1);
    vetorERRO = zeros(obj.num_individuos, 1);
    for i = 1:obj.num_individuos
        vetorFIT(i) = obj.populacao(i).FITNESS;
        vetorERRO(i) = abs(obj.populacao(i).erro);
    end

    obj.maiorFIT = max(vetorFIT);
    
    % Verifica se todos erros sao menores do que o erro de parada, e se
    % forem, aumenta o numero de rodadas que o erro foi menor
    if(max(vetorERRO) <= obj.erro_parada)
        obj.quant_erro = obj.quant_erro + 1;
        
        % Finaliza o algoritmo se o erro se perpetuou durante geracoes
        % suficientes
        if(obj.quant_erro == obj.quant_geracos_erro)
            obj.finalizou = true;
        end
    else
        % Reseta o contador de geracos de erro
        obj.quant_erro = 0;
    end
    
    % Indica que houve uma nova rodada
    obj.rodadas = obj.rodadas + 1;
    
    obj.maiorFITGer(obj.rodadas) = obj.maiorFIT;
    
    % Finaliza por ter chegado ao numero limite de geracoes
    if(obj.rodadas >= obj.max_geracoes)
        obj.finalizou = true;
    end
    
    % DEBUG
    %obj.mostrarResultados();
    fprintf("Fim da Geracao %d\n-----------------------\n", obj.rodadas) % Numero de geracoes
end