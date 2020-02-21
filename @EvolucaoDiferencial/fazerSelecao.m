% Faz a selecao entre os individuos
% Retorna se foi possivel fazer a selecao ou nao
function selecionou = fazerSelecao(obj, individuo, mutacao)
    
    % A funcao retorna os valores da FO e tambem como ultimo elemento se a
    % funcao teve erro ou nao
    FOeErro = calcularFuncaoObjetivo(obj, mutacao);
    
    % Verifica se teve erro
    if(FOeErro(obj.funcao_objetivo.num_saidas+1) == 1)
        selecionou = false;
        
    % Nao teve erro
    else
    
        % Obtem somente os parametros da FO, sem o erro, e o FITNESS
        FO = FOeErro(1:obj.funcao_objetivo.num_saidas);
        FITNESS = calcularFitness(obj, FO);
        FITNESS_ANTIGA = obj.populacao(individuo).FITNESS;
        
        % Se o FITNESS for melhor do que o individuo atual, faz a troca
        if(FITNESS >= FITNESS_ANTIGA)
            
            % Indica o quanto o individuo aprimorou do individuo passado
            obj.populacao(individuo).erro = ...
                (FITNESS-FITNESS_ANTIGA)/FITNESS_ANTIGA;
            obj.populacao(individuo).FO = FO;
            obj.populacao(individuo).FITNESS = FITNESS;
            obj.populacao(individuo).parametros = mutacao;
        end
        
        selecionou = true;
    end
end