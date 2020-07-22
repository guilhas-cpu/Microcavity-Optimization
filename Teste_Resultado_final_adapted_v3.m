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
    isLambdaR1 = 0; %valor default
    dispersion = 2; %valor default
    iter = 0;
    s = 1;
    load(dado);

    %result = Reflectance_adapted(vars); %registra os resultados dos objetivos
    result = Reflectance_final_results(vars,s); %registra os resultados dos objetivos
    melhora = [result(1)/original(1) result(2)/original(2) result(3)/original(3) original(4)/result(4) original(5)/result(5) original(6)/result(6) original(7)/result(7)] - 1; %registra a melhora em relação aos resultados base (enviado pelo professor cotta)
    melhora = melhora*100;
    if(size(vars,1)==6)
        shift = vars(6); %Cavity shift [nm]
        %lambdaR1 = vars(6);
    else
        shift = 0;
        %lambdaR1 = 898;
    end
    fileID = fopen(name,'a');
    z = result;
    
    fprintf(fileID,'\n\n\ndado: %s',dado);
    fprintf(fileID,'\nVARIAVEIS DE ENTRADA:');
    if(isLambdaR1)
        fprintf(fileID,'\nc1 = %.2f      c2 = %.2f      c3= %.2f      ncs= %d      nci= %d      shift= 0      lambdaR1= %d      range=  %d      s= %d',vars(1),vars(2),vars(3),vars(4),vars(5),shift+898,range,s);
    else
        fprintf(fileID,'\nc1 = %.2f      c2 = %.2f      c3= %.2f      ncs= %d      nci= %d      shift= %d      lambdaR1= 898      range=  %d      s= %d',vars(1),vars(2),vars(3),vars(4),vars(5),shift,range,s);
    end
    fprintf(fileID,'\nRESULTADOS BASE\t\t\t\t\t\tNOVOS RESULTADOS\t\t\t\t\tPERCENTUAL DE MELHORA\n');
    fprintf(fileID,'beta = %d\t\t\t\t\tbeta = %d \t\t\t\t\t(%.2f%%)\n',original(1),result(1),melhora(1));
    fprintf(fileID,'Q = %d\t\t\t\t\tQ = %d \t\t\t\t\t(%.2f%%)\n',original(2),result(2),melhora(2));
    fprintf(fileID,'Fp = %d\t\t\t\t\tFp = %d \t\t\t\t\t(%.2f%%)\n',original(3),result(3),melhora(3));
    fprintf(fileID,'DN = %d\t\t\t\t\tDN = %d \t\t\t\t\t(%.2f%%)\n',original(4),result(4),melhora(4));
    fprintf(fileID,'CavityLoss = %d\t\t\t\tCavityLoss = %d \t\t\t\t(%.2f%%)\n',original(5),result(5),melhora(5));
    fprintf(fileID,'Distance = %d\t\t\t\t\tDistance = %d \t\t\t\t(%.2f%%)\n',original(6),result(6),melhora(6));
    fprintf(fileID,'Reflec_min= %d\t\t\t\tReflec_min = %d \t\t\t\t(%.2f%%)\n',original(7),result(7),melhora(7)); 
    fprintf(fileID,'PESOS: %s',evalc('disp(PESO)'));
    fprintf(fileID,'isLambdaR1= %d, dispesion= %d, METODO= %d, range= %d, s= %d\n',isLambdaR1,dispersion,METODO,range,s);
    
    fclose(fileID);
    

    
    
    str1 = sprintf('\ndado: %s',dado);
    str2 = sprintf('\nVARIAVEIS DE ENTRADA:');
    if(isLambdaR1)
        str3 = sprintf('\nc1 = %.2f      c2 = %.2f      c3= %.2f      ncs= %d      nci= %d      shift= 0      lambdaR1= %d      range=  %d      s= %d',vars(1),vars(2),vars(3),vars(4),vars(5),shift+898,range,s);
    else
        str3 = sprintf('\nc1 = %.2f      c2 = %.2f      c3= %.2f      ncs= %d      nci= %d      shift= %d      lambdaR1= 898      range=  %d      s= %d',vars(1),vars(2),vars(3),vars(4),vars(5),shift,range,s);
    end
    str4 = sprintf('\n\n                           RESULTADOS BASE                                        NOVOS RESULTADOS                          PERCENTUAL DE MELHORA\n');
    str5 = sprintf('beta = %.2d                                                 beta = %.2d                                        (%.2f%%)\n',original(1),result(1),melhora(1));
    str6 = sprintf('Q = %.2d                                                    Q = %.2d                                            (%.2f%%)\n',original(2),result(2),melhora(2));
    str7 = sprintf('Fp = %.2d                                                    Fp = %.2d                                           (%.2f%%)\n',original(3),result(3),melhora(3));
    str8 = sprintf('DN = %.2d                                                  DN = %.2d                                         (%.2f%%)\n',original(4),result(4),melhora(4));
    str9 = sprintf('CavityLoss = %.2d                                      CavityLoss = %.2d                              (%.2f%%)\n',original(5),result(5),melhora(5));
    str10 = sprintf('Distance = %.2d                                         Distance = %.2d                                 (%.2f%%)\n',original(6),result(6),melhora(6));
    str11 = sprintf('Reflec_min= %.2d                                      Reflec_min = %.2d                              (%.2f%%)\n',original(7),result(7),melhora(7)); 
    str = [str1,str2,str3,str4,str5,str6,str7,str8,str9,str10,str11];
    

    figure(bons)
    dim = [.1 .02 .8 .33];
    a = annotation('textbox',dim,'String',str);
    a.HorizontalAlignment = 'center';
    a.VerticalAlignment = 'middle';
end
