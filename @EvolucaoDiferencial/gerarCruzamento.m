function vetor = gerarCruzamento(obj, parametros_individuo, vetor_ruido)
    
    vetor = zeros(obj.num_parametros, 1);
    
    % Seleciona pelo menos 1 posicao para realizar mutacao, de forma
    % garantida
    i = randi([1 obj.num_parametros]);
    vetor(i) = vetor_ruido(i);

    % Faz a mutacao
    for i = 1:obj.num_parametros
        
        % Se o numero aleatorio gerado for menor ou igual a mutacao,
        % pega do vetor_ruido, se nao, pega do individuo
        if(rand() <= obj.chance_mutacao)
            vetor(i) = vetor_ruido(i);
        else
            vetor(i) = parametros_individuo(i);
        end
        
    end
end