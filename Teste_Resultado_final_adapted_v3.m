function fim = Teste_Resultado_final_adapted_v3(dado,name,original,plotar)
    global print;
    global plota;
    global range;
    global iter;
    global bons;
    warning off;
    print = 0;
    plota = plotar;
    range = 50; %valor default
    iter = 0;
    load(dado);

    result = Reflectance_final_results(vars); %registra os resultados dos objetivos
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
    
    
    str1 = sprintf('dado: %s',dado);
    str2 = sprintf('\n\n                           RESULTADOS BASE                                        NOVOS RESULTADOS                          PERCENTUAL DE MELHORA\n');
    str3 = sprintf('beta = %.2d                                                 beta = %.2d                                        (%.2f%%)\n',original(1),result(1),melhora(1));
    str4 = sprintf('Q = %.2d                                                    Q = %.2d                                            (%.2f%%)\n',original(2),result(2),melhora(2));
    str5 = sprintf('Fp = %.2d                                                    Fp = %.2d                                           (%.2f%%)\n',original(3),result(3),melhora(3));
    str6 = sprintf('DN = %.2d                                                  DN = %.2d                                         (%.2f%%)\n',original(4),result(4),melhora(4));
    str7 = sprintf('CavityLoss = %.2d                                      CavityLoss = %.2d                              (%.2f%%)\n',original(5),result(5),melhora(5));
    str8 = sprintf('Distance = %.2d                                         Distance = %.2d                                 (%.2f%%)\n',original(6),result(6),melhora(6));
    str9 = sprintf('Reflec_min= %.2d                                      Reflec_min = %.2d                              (%.2f%%)\n',original(7),result(7),melhora(7)); 
    str = [str1,str2,str3,str4,str5,str6,str7,str8,str9];
    

    figure(bons)
    dim = [.1 .05 .8 .30];
    a = annotation('textbox',dim,'String',str);
    a.HorizontalAlignment = 'center';
    a.VerticalAlignment = 'middle';
    warning on;
end
