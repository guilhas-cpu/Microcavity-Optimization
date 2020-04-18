close all;
clc;
clear all;
global print;
global plota;
global range;
global erro;
print = 0;
plota = 0;
range = 20;
erro = 0;

ds = datastore("./Erros",'FileExtensions',{'.mat'},'Type','tabulartext','TextscanFormats',{'%f','%f'});
tamanho = size(ds.Files);


for i=1:tamanho(1)
   dado = string(ds.Files(i)); 
   load(dado)
   dado = dado;
   vars = vars
end
%z = Reflectance_adapted_v2(vars)