% Gera a mutacao de determinado individuo
% Retorna o vetor ruido criado
function vetor_ruido = gerarMutacao(obj, vetor_base, indice_vetor_base)

    % Gera 2 numeros aleatorios para escolher 2 individuos da nossa
    % populacao. Os individuos devem ser diferentes entre si e tambem do
    % vetor base
    while 1
        indices = randi([1 obj.num_parametros], 2);

        if(indice_vetor_base ~= indices(1) && ...
           indice_vetor_base ~= indices(2) && ...
                indices(1) ~= indices(2))
            break
        end
    end

    % Gera um vetor ruido a partir do vetor base, e a diferença vezes um
    % peso entre os dois individuos aleatorios obtidos da populacao
    vetor_ruido = vetor_base + ...
        obj.F*(obj.populacao(indices(1)).parametros-...
        obj.populacao(indices(2)).parametros);
    
    % Verifica se existe posicoes com tipo inteiro, e se existir,
    % o peso deve ser inteiro e a diferença deve ser obrigatoriamente
    % positiva
    for i = 1:obj.num_parametros
        
        % Tipo inteiro
        if(obj.parametros(i).tipo == 0)
            % Gera o vetor ruido inteiro
            vetor_ruido(i) = vetor_base(i) + ...
            abs(randi([1 2])*(obj.populacao(indices(1)).parametros(i)-...
            obj.populacao(indices(2)).parametros(i)));
        
            limi = obj.parametros(i).liminferior;
            lims = obj.parametros(i).limsuperior;
            
            % Se exceder o limite inferior, gera um valor dentro do limite
            % entre o limite inferior e metade dos limites
            if(limi > vetor_ruido(i))
                vetor_ruido(i) = randi([limi (floor((lims-limi)/2)+limi)]);
            
            % Mesma logica para o superior
            elseif(lims < vetor_ruido(i))
                vetor_ruido(i) = randi([(lims-floor((lims-limi)/2)) lims]);
            end
        end
    end
    
end