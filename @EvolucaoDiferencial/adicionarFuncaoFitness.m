function adicionarFuncaoFitness(obj, funcao_fitness)

    % Verifica erros
    assert(isa(funcao_fitness, 'function_handle'), ...
        "[ERRO] funcao_fitness deve ser do tipo 'function_handle'");

    obj.funcao_fitness = funcao_fitness;
end

