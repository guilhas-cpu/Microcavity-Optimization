function fim = Teste_Resultado_final_adapted_v2(dado,name,original)

    global print;
    global plota;
    global range;
    global iter;
    print = 0;
    plota = 0;
    range = 50; %valor default
    iter = 0;
    load(dado);

    result = Reflectance_adapted(vars); %registra os resultados dos objetivos
    melhora = [result(1)/original(1) result(2)/original(2) result(3)/original(3) original(4)/result(4) original(5)/result(5) original(6)/result(6) original(7)/result(7)] - 1; %registra a melhora em relação aos resultados base (enviado pelo professor cotta)
    melhora = melhora*100;

    fileID = fopen(name,'a');
    z = result;
    
    fprintf(fileID,'\n\n\ndado: %s',dado)
    fprintf(fileID,'\nRESULTADOS BASE\t\t\t\t\t\tNOVOS RESULTADOS\t\t\t\t\tPERCENTUAL DE MELHORA\n');
    fprintf(fileID,'beta = %d\t\t\t\t\tbeta = %d \t\t\t\t\t(%.2f%%)\n',original(1),result(1),melhora(1));
    fprintf(fileID,'Q = %d\t\t\t\t\tQ = %d \t\t\t\t\t(%.2f%%)\n',original(2),result(2),melhora(2));
    fprintf(fileID,'Fp = %d\t\t\t\t\tFp = %d \t\t\t\t\t(%.2f%%)\n',original(3),result(3),melhora(3));
    fprintf(fileID,'DN = %d\t\t\t\t\tDN = %d \t\t\t\t\t(%.2f%%)\n',original(4),result(4),melhora(4));
    fprintf(fileID,'CavityLoss = %d\t\t\t\tCavityLoss = %d \t\t\t\t(%.2f%%)\n',original(5),result(5),melhora(5));
    fprintf(fileID,'Distance = %d\t\t\t\t\tDistance = %d \t\t\t\t(%.2f%%)\n',original(6),result(6),melhora(6));
    fprintf(fileID,'Reflec_min= %d\t\t\t\tReflec_min = %d \t\t\t\t(%.2f%%)\n',original(7),result(7),melhora(7)); 
    
    fclose(fileID);
    
    if(strcmp(name,'restore.txt'))
        save(dado,'vars','z','norm_z','PESO','countSuccesses','countErrors','METODO','maior','range','iter')
    end
end
