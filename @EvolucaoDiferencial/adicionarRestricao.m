function adicionarRestricao(obj, funcao_restricao, varargin)

    % Verifica erros
    assert(isa(funcao_restricao, 'function_handle'), ...
        "[ERRO] funcao_restricao deve ser do tipo 'function_handle'");
    
    obj.num_restricoes = obj.num_restricoes + 1;
    
    % Armazena a funcao restricao
    obj.funcao_restricao(obj.num_restricoes).funcao = funcao_restricao;
    obj.funcao_restricao(obj.num_restricoes).restricoes = zeros(nargin-2, 1);
    
    % Busca os parametros e armazena seus identificadores
    for i = 1:nargin-2
        
        encontrou = false;
        for j = 1:obj.num_parametros
            if varargin{i} == obj.parametros(j).nome
                encontrou = true;
                break;
            end
        end
        
        assert(encontrou, ...
        sprintf("[ERRO] O parametro %s nao foi encontrado", varargin{i}))
    
        % Armazena onde esta o argumento
        obj.funcao_restricao(obj.num_restricoes).restricoes(i) = j;
    end
end