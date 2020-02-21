% Esta funcao verifica se os parametros dados da populacao respeitam as restricoes
% Retorno: restricoes respeitadas ou nao
% Esta funcao verifica se os parametros dados da populacao respeitam as restricoes
% Retorno: restricoes respeitadas ou nao
function y = verificarRestricoes(obj, parametros)

    % Verifica se os parametros estao nos limites
    for i=1:obj.num_parametros
        if(parametros(i) > obj.parametros(i).limsuperior || ...
            parametros(i) < obj.parametros(i).liminferior)
            y = false;
            return;
        end
    end

    % Verifica se as funcoes de restricao foram atendidas
    for i=1:size(obj.funcao_restricao, 1)
        
        % Armazena os parametros necessarios para dada restricao
        posparametros = zeros(size(obj.funcao_restricao(i).restricoes, 1), 1);
        posparametros(:) = obj.funcao_restricao(i).restricoes(:);
        valor_parametros = zeros(size(posparametros, 1), 1);
        valor_parametros(:) = parametros(posparametros(:));
        
        % Executa a funcao de restricao com o valor de parametros lido
        if(~obj.funcao_restricao(i).funcao(valor_parametros))
            y = false;
            return;
        end
    end
    
    y = true;
    
end