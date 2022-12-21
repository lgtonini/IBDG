    
clear all
close all
clc 

%Entrada de dos do IEEE - 123 - Regulador 1

%Coordenadas/Configura��o (1-N�, 2-Linha, 3-Transformador, 4-Regulador e 9-N� auxiliar)
M = [0,0,1,0,0,0,0;0,0,2,0,0,0,0;0,0,1,2,1,2,1;0,0,2,0,0,0,0;1,2,1,2,1,3,1;0,0,2,0,0,0,0;0,0,1,0,0,0,0;0,0,2,0,0,0,0;0,0,1,0,0,0,0;0,0,2,0,0,0,0;0,0,1,0,0,0,0;0,0,2,0,0,0,0;0,0,1,0,0,0,0;0,0,2,0,0,0,0;0,0,1,0,0,0,0];

%Quantidade de configura��es
NC = 4; 

%Configura��es Y e Z
MC_Y(1:3,7:9) = (10^-6)*[5.3971 -0.6982 -1.1645; 0 5.6765 -1.8319; 0 0 5.9809];   %Configura��o 3
MC_Z(1:3,7:9) = [0.4615+j*1.0651 0.1535+j*0.3849 0.1580+j*0.4236; 0 0.4576+j*1.0780 0.1560+j*0.5017; 0 0 0.4666+j*1.0482];

MC_Y(1:3,13:15) = (10^-6)*[5.9809 -1.8319 -1.1645; 0 5.6765 -0.6982; 0 0 5.3971]; %Configura��o 5
MC_Z(1:3,13:15) = [0.4666+j*1.0482 0.1560+j*0.5017 0.1580+j*0.4236; 0 0.4576+j*1.0780 0.1535+j*0.3849; 0 0 0.4615+j*1.0651];

MC_Y(1:3,28:30) = (10^-6)*[0 0 0; 0 4.5193 0; 0 0 0];                             %Configura��o 10
MC_Z(1:3,28:30) = [0 0 0; 0 1.3292+j*1.3475 0; 0 0 0];

MC_Y(1:3,34:36) = (10^-6)*[67.2242 0 0; 0 67.2242 0; 0 0 67.2242];                %Configura��o 12
MC_Z(1:3,34:36) = [1.5209+j*0.7521 0.5198+j*0.2775 0.4924+j*0.2157; 0 1.5329+j*0.7162 0.5198+j*0.2775; 0 0 1.5209+j*0.7521];

%Dist�ncias e configura��o entre os n�s
Dist(64,1)=350*0.000189394;            %54-57 - 1
Dist(64,2)=1;             %Configura��o 1-3

Dist(66,1)=250*0.000189394;            %57-58 - 2
Dist(66,2)=3;            %Configura��o 3-10

Dist(65,1)=250*0.000189394;            %58-59 - 3
Dist(65,2)=3;             %Configura��o 3-10

Dist(68,1)=750*0.000189394;            %57-60 - 4
Dist(68,2)=1;            %Configura��o 1-3

Dist(77,1)=0;                           %60-160 - 5
Dist(77,2)=1;             %Configura��o 1-3

Dist(70,1)=550*0.000189394;            %60-61 - 6
Dist(70,2)=2;             %Configura��o 2-5 

Dist(69,1)=0;            %61-610 - 7
Dist(69,2)=2;             %Configura��o 2-5 (Transformador)     

Dist(71,1)=250*0.000189394;             %60-62 - 8
Dist(71,2)=4;             %Configura��o 4-12

Dist(72,1)=175*0.000189394;            %62-63 - 9
Dist(72,2)=4;             %Configura��o 4-12 

Dist(73,1)=350*0.000189394;             %63-64 - 10
Dist(73,2)=4;             %Configura��o 4-12 

Dist(74,1)=425*0.000189394;            %64-65 - 11
Dist(74,2)=4;             %Configura��o 4-12

Dist(75,1)=325*0.000189394;             %65-66 - 13
Dist(75,2)=4;             %Configura��o 4-12

%Pot�ncia nos n�s em estrela e delta (TC = Tipo de carga: 1-S, 2-I ou 3-Z)
%PPA = Percentual de Pot�ncia Ativa
PPA = 0             %N� 54 - 1
SN(1:3,23)=[0,0,0];                                                                                                        
SN(4:6,23)=[0,0,0];

PPA = 0             %N� 57 - 2
SN(1:3,67)=[0,0,0];                                                                                                      
SN(4:6,67)=[0,0,0];

PPA = 0             %N� 58 - 3
SN(1:3,66)=[0,(20+j*10)*1e3,0];                                                                                   
TC_Y(66) = 2;
SN(4:6,66)=[0,0,0];

PPA = 0             %N� 59 - 4
SN(1:3,65)=[0,((20+j*10)*1e3),0];                                                                         
TC_Y(65) = 1;
SN(4:6,65)=[0,0,0];

PPA = 0             %N� 160 - 5
SN(1:3,78)=[0,0,0];
SN(4:6,78)=[0,0,0];

PPA = 0             %N� 60 - 6
SN(1:3,71)=[(20+j*10)*1e3,0,0];                                                                                         
TC_Y(71) = 1;
SN(4:6,71)=[0,0,0];

PPA = 0             %N� 61 - 7
SN(1:3,70)=[0,0,0];                                                                                         
SN(4:6,70)=[0,0,0];

PPA = 0             %N� 610 - 8 
SN(1:3,69)=[0,0,0];                                                                                                                          
SN(4:6,69)=[0,0,0];

PPA = 0             %N� 62 - 9
SN(1:3,72)=[0,0,((40+j*20)*1e3)];                                                                                         
TC_Y(72) = 3;
SN(4:6,72)=[0,0,0];

PPA = 0             %N� 63 - 10
SN(1:3,73)=[((40+j*20)*1e3),0,0];
TC_Y(73) = 1;
SN(4:6,73)=[0,0,0];

PPA = 0             %N� 64 - 11 
SN(1:3,74)=[0,((75+j*35)*1e3),0];                                                                                                
TC_Y(74) = 1;
SN(4:6,74)=[0,0,0];

PPA = 0             %N� 65 - 12 
SN(1:3,75)=[0,0,0];                                                                                
SN(4:6,75)=[((35+j*25)*1e3),((35+j*25)*1e3),((70+j*50)*1e3)];
TC_D(75) = 3;

PPA = 0             %N� 66 - 13
SN(1:3,76)=[0,0,((75+j*35)*1e3)];                                                                                                                 
TC_Y(76) = 1;
SN(4:6,76)=[0,0,0];

%Transformadores:
    %Configura��es 
1;
    %Tipo de transformador(1-DYg e 2-YgYg)
2;
    %Pot�ncia(VA)
150*10^3;
    %Tens�o de entrada(V)
4.16*10^3;
    %Tens�o de sa�da(V)
0.48*10^3;
    %Imped�cia
0.0127+1i*0.0272;

%Tens�o de refer�ncia
Vref = (4.16*10^3)/sqrt(3);

%Erro
Erro = 0.001;

