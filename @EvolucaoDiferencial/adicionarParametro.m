function adicionarParametro(obj, nome_parametro, tipo, limites)
    
    assert(strcmp(tipo, 'real') || strcmp(tipo, 'inteiro'),...
    "[ERRO] O tipo do parametro deve ser real ou inteiro")

    for i = 1:obj.num_parametros
        er = sprintf("[ERRO] O parametro %s ja existe", nome_parametro);
        assert(~strcmp(obj.parametros(i).nome, nome_parametro), er)
    end
    
    obj.num_parametros = obj.num_parametros + 1;
    
    if(strcmp(tipo, 'inteiro'))
        t = 0;
    else
        t = 1;
    end
    
    % Armazena as informacoes do parametro
    obj.parametros(obj.num_parametros).nome = nome_parametro;
    obj.parametros(obj.num_parametros).tipo = t;
    obj.parametros(obj.num_parametros).liminferior = limites(1,1);
    obj.parametros(obj.num_parametros).limsuperior = limites(1,2);
end