function adicionarFuncaoObjetivo(obj, funcao_objetivo, varargin)
    
    % Verifica erros
    assert(isa(funcao_objetivo, 'function_handle'), ...
        "[ERRO] funcao_objetivo deve ser do tipo 'function_handle'");
    
    % Armazena a funcao objetivo
    obj.funcao_objetivo.funcao = funcao_objetivo;
    obj.funcao_objetivo.nome_saidas = varargin;
    obj.funcao_objetivo.num_saidas = nargin-2;
end