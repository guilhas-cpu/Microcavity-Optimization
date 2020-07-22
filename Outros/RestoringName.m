clear all;
clc;
close all;
ds = datastore("./Dados",'FileExtensions',{'.jpg','.mat'},'Type','tabulartext');

tamanho = size(ds.Files);

for i=1:tamanho(1)
   dado = char(string(ds.Files(i))); 
   i = i
   load(dado);
   sizeDado = size(dado);
   
   if (isletter(dado(sizeDado(2)-8)))
       novoDado = insert(dado,sizeDado(2)-8,'0');
       novoDado = string(novoDado);
       save(novoDado,'vars','z','norm_z','PESO','countSuccesses','countErrors','METODO','maior','range','iter')
       delete(dado);
       dado = char(novoDado);
   end
   sizeDado = size(dado);
   if ((dado(sizeDado(2)-11))=='-')
       novoDado = insert(dado,sizeDado(2)-11,'0');
       novoDado = string(novoDado);
       save(novoDado,'vars','z','norm_z','PESO','countSuccesses','countErrors','METODO','maior','range','iter')
       delete(dado);
       dado = char(novoDado);
   end
   sizeDado = size(dado);
   if ((dado(sizeDado(2)-14))=='_')
       novoDado = insert(dado,sizeDado(2)-14,'0');
       novoDado = string(novoDado);
       save(novoDado,'vars','z','norm_z','PESO','countSuccesses','countErrors','METODO','maior','range','iter')
       delete(dado);
       dado = char(novoDado);
   end
   sizeDado = size(dado);
   if ((dado(sizeDado(2)-16))=='_')
       novoDado = insert(dado,sizeDado(2)-16,'0');
       novoDado = string(novoDado);
       save(novoDado,'vars','z','norm_z','PESO','countSuccesses','countErrors','METODO','maior','range','iter')
       delete(dado);
       dado = char(novoDado);
   end
    sizeDado = size(dado);
   if ((dado(sizeDado(2)-17))=='_')
       novoDado = insert(dado,sizeDado(2)-17,'0');
       novoDado = string(novoDado);
       save(novoDado,'vars','z','norm_z','PESO','countSuccesses','countErrors','METODO','maior','range','iter')
       delete(dado);
       dado = char(novoDado);
   end
end


        

