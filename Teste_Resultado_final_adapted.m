function fim = Teste_Resultado_final_adapted(dado,name)

    global print;
    global plota;
    global range;
    global iter;
    print = 0;
    plota = 0;
    range = 30; %valor default
    iter = 0;
    maior = zeros(1,7);
    METODO = 1; %valor default
    load(dado);

    result = Reflectance_adapted(vars); %registra os resultados dos objetivos
    original = [1.183741e-02 5.281276e+03 1.718943e-01 9.397523e+08 9.546695e-08 9.381136e-01 6.076570e-01]; 
    melhora = [result(1)/original(1) result(2)/original(2) result(3)/original(3) original(4)/result(4) original(5)/result(5) original(6)/result(6) original(7)/result(7)] - 1; %registra a melhora em relação aos resultados base (enviado pelo professor cotta)
    melhora = melhora*100;

    fileID = fopen(name,'a');
    
    fprintf(fileID,'\n\n\ndado: %s',dado)
    fprintf(fileID,'\nRESULTADOS BASE\t\t\t\t\t\tNOVOS RESULTADOS\t\t\t\t\tPERCENTUAL DE MELHORA\n');
    fprintf(fileID,'beta = 1.183741e-02\t\t\t\t\tbeta = %d \t\t\t\t\t(%f%%)\n',result(1),melhora(1));
    fprintf(fileID,'Q = 5.281276e+03\t\t\t\t\tQ = %d \t\t\t\t\t(%f%%)\n',result(2),melhora(2));
    fprintf(fileID,'Fp = 1.718943e-01\t\t\t\t\tFp = %d \t\t\t\t\t(%f%%)\n',result(3),melhora(3));
    fprintf(fileID,'DN = 9.397523e+08\t\t\t\t\tDN = %d \t\t\t\t\t(%f%%)\n',result(4),melhora(4));
    fprintf(fileID,'CavityLoss = 9.546695e-08\t\t\t\tCavityLoss = %d \t\t\t\t(%f%%)\n',result(5),melhora(5));
    fprintf(fileID,'Distance = 9.381136e-01\t\t\t\t\tDistance = %d \t\t\t\t(%f%%)\n',result(6),melhora(6));
    fprintf(fileID,'Reflec_min= 6.076570e-01\t\t\t\tReflec_min = %d \t\t\t\t(%f%%)\n',result(7),melhora(7)); 
    
    fclose(fileID);
    
    %save(dado,'vars','z','norm_z','PESO','countSuccesses','countErrors','METODO','maior','range','iter')
end
