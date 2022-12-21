
clear all
close all
clc 

%Entrada de dos do IEEE - 123 - Regulador 2

%Coordenadas/Configura��o (1-N�, 2-Linha, 3-Transformador, 4-Regulador e 9-N� auxiliar)
M = [1,0,0;2,0,0;1,2,1;2,0,0;1,0,0];

%Quantidade de configura��es
NC = 1;

%Configura��es Y e Z
MC_Y(1:3,1:3) = (10^-6)*[4.5193 0 0; 0 0 0; 0 0 0];    %Configura��o 9
MC_Z(1:3,4:6) = [1.3292+1i*1.3475 0 0; 0 0 0; 0 0 0];

%Dist�ncias e configura��o entre os n�s

Dist(1,1)=425*0.000189394;           %9-14 - 1-2 - Regulador 2 (9,23-9,25) - 11
Dist(1,2)=1;            %Configura��o 9

Dist(2,1)=250*0.000189394;             %14-10 - 2-3 (9,21-9,23) - 10
Dist(2,2)=1;             %Configura��o 9

Dist(3,1)=250*0.000189394;           %14-11 - 2-4 (9,23-11,23) - 15
Dist(3,2)=1;            %Configura��o 9

%Pot�ncia nos n�s em estrela e delta (TC = Tipo de carga: 1-S, 2-I ou 3-Z)
%PPA = Percentual de Pot�ncia Ativa
PPA = 0             %N� 9 - 1
SN(1:3,1)=[((40+1i*20)*1e3),((40+1i*20)*1e3),0];                                                                                           
TC_Y(1) = 1;
SN(4:6,1)=[0,0,0];

PPA = 0             %N� 14 - 2
SN(1:3,2)=[0,((40+1i*20)*1e3),0]; 
TC_Y(2) = 1;
SN(4:6,2)=[0,0,0];

PPA = 0             %N� 10 - 3
SN(1:3,3)=[((20+1i*10)*1e3),0,0];                                                 
TC_Y(3) = 2;
SN(4:6,3)=[0,0,0];

PPA = 0             %N� 11 - 4
SN(1:3,4)=[((40+1i*20)*1e3),0,0];                                                                                                         
TC_Y(4) = 3;
SN(4:6,4)=[0,0,0];

%Tens�o de refer�ncia
Vref = (4160)/sqrt(3);

%Erro
Erro = 0.001;

