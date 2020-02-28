close all;
clc;
clear all;
global print;
global plota;
global range;
global erro;
print = 1;
plota = 1;
range = 20;
erro = 0;
load('./Erros/erro5.mat')


z = Reflectance_adapted(vars)