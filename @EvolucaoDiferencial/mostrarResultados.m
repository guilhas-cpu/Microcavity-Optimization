function [FO, individuo] = mostrarResultados(obj)

    FOs = zeros(obj.num_individuos, obj.funcao_objetivo.num_saidas);
    FIs = zeros(obj.num_individuos, 1);
    ERROS = zeros(obj.num_individuos, 1);
    
    fprintf("Individuos da ultima geracao:\n");
    for i = 1:obj.num_individuos
        for j = 1:obj.num_parametros
            fprintf("%s: %f ", obj.parametros(j).nome, ...
                 obj.populacao(i).parametros(j))
        end
        FOs(i, :) = obj.populacao(i).FO(:);
        FIs(i) = obj.populacao(i).FITNESS;
        ERROS(i) = obj.populacao(i).erro;        
        fprintf(" FI: %f ERRO: %f%%\n", FIs(i), ERROS(i));
    end
    
    [maiorFI, i] = max(FIs);
    
    fprintf("\nMelhor Individuo:\n");
    for j = 1:obj.num_parametros
        fprintf("%s: %f ", obj.parametros(j).nome, ...
             obj.populacao(i(1)).parametros(j))
    end
    fprintf(" FI: %f ERRO: %f%%\n", maiorFI, ERROS(i(1)));
    fprintf("Numero de geracoes: %d\n", obj.rodadas);
    
    FO = FOs(i(1), :);
    individuo = i(1);
end