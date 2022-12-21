
clear all
close all
clc 

%Entrada de dos do IEEE - 123 - Regulador 2

%Coordenadas/Configura��o (1-N�, 2-Linha, 3-Transformador, 4-Regulador e 9-N� auxiliar)
M = [1,0,0,0,0;2,0,0,0,0;1,2,1,2,1;2,0,0,0,0;1,0,0,0,0;2,0,0,0,0;1,0,0,0,0];

%Quantidade de configura��es
NC = 3;

%Configura��es Y e Z
MC_Y(1:3,19:21) = (10^-6)*[5.1154 0 -1.0549; 0 0 0; 0 0 5.1704];                  %Configura��o 7
MC_Z(1:3,19:21) = [0.4676+j*1.0780 0 0.1535+j*0.3849; 0 0 0; 0 0 0.4515+j*1.0651];

MC_Y(1:3,25:27) = (10^-6)*[4.5193 0 0; 0 0 0; 0 0 0];                             %Configura��o 9
MC_Z(1:3,25:27) = [1.3292+j*1.3475 0 0; 0 0 0; 0 0 0];

MC_Y(1:3,31:33) = (10^-6)*[0 0 0; 0 0 0; 0 0 4.5193];                             %Configura��o 11
MC_Z(1:3,31:33) = [0 0 0; 0 0 0; 0 0 1.3292+j*1.3475];

%Dist�ncias e configura��o entre os n�s
Dist(56,1)=350*0.000189394;           %26-25 - Regulador 3 - 1
Dist(56,2)=1;             %Configura��o 1-7

Dist(57,1)=275*0.000189394;           %26-27 - 2
Dist(57,2)=1;             %Configura��o 1-7

Dist(62,1)=500*0.000189394;             %27-33 - 3
Dist(62,2)=2;             %Configura��o 2-9

Dist(55,1)=225*0.000189394;             %31-26 - 4
Dist(55,2)=3;             %Configura��o 3-11

Dist(54,1)=300*0.000189394;            %32-31 - 5
Dist(54,2)=3;             %Configura��o 3-11

%Pot�ncia nos n�s em estrela e delta (TC = Tipo de carga: 1-S, 2-I ou 3-Z)
%PPA = Percentual de Pot�ncia Ativa
PPA = 0             %N� 25 - 1
SN(1:3,56)=[0,0,0];                                                                                                                   
SN(4:6,56)=[0,0,0];

PPA = 0             %N� 26 - 2
SN(1:3,55)=[0,0,0];
SN(4:6,55)=[0,0,0];

PPA = 0             %N� 27 - 3
SN(1:3,58)=[0,0,0];                                                                                                                      
SN(4:6,58)=[0,0,0];

PPA = 0             %N� 33 - 4 
SN(1:3,63)=[((40+j*20)*1e3),0,0];                                                                                                                    
TC_Y(63) = 2;
SN(4:6,63)=[0,0,0];

PPA = 0             %N� 31 - 5
SN(1:3,54)=[0,0,((20+j*10)*1e3)];                                                                                                          
TC_Y(54) = 1;
SN(4:6,54)=[0,0,0];

PPA = 0             %N� 32 - 6
SN(1:3,53)=[0,0,((20+j*10)*1e3)];                                                                                                         
TC_Y(53) = 1;
SN(4:6,53)=[0,0,0];

%Tens�o de refer�ncia
Vref = (4160)/sqrt(3);

%Erro
Erro = 0.001;

