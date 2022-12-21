    
clear all
close all
clc 

%Entrada de dos do IEEE - 123 - Regulador 1

%Coordenadas/Configura��o (1-N�, 2-Linha, 3-Transformador, 4-Regulador e 9-N� auxiliar)
M = [1,0,0,0,0,0,0,0,0;2,0,0,0,0,0,0,0,0;1,2,1,2,1,2,1,2,1;2,0,0,0,0,0,0,0,0;1,2,1,2,1,2,1,2,1;2,0,0,0,0,0,0,0,0;1,0,0,0,0,0,0,0,0;2,0,0,0,0,0,0,0,0;1,2,1,2,1,2,1,0,0;2,0,0,0,0,0,0,0,0;1,2,1,2,1,0,0,0,0;2,0,0,0,0,0,0,0,0;1,2,1,0,0,0,0,0,0;2,0,0,0,0,0,0,0,0;1,0,0,0,0,0,0,0,0;2,0,0,0,0,0,0,0,0;1,2,1,0,0,0,0,0,0;2,0,0,0,0,0,0,0,0;1,0,0,0,0,0,0,0,0;2,0,0,0,0,0,0,0,0;1,0,0,0,0,0,0,0,0;2,0,0,0,0,0,0,0,0;1,0,0,0,0,0,0,0,0];

%Quantidade de configura��es
NC = 5; 

%Configura��es Y e Z
MC_Y(1:3,7:9) = (10^-6)*[5.3971 -0.6982 -1.1645; 0 5.6765 -1.8319; 0 0 5.9809];   %Configura��o 3
MC_Z(1:3,7:9) = [0.4615+j*1.0651 0.1535+j*0.3849 0.1580+j*0.4236; 0 0.4576+j*1.0780 0.1560+j*0.5017; 0 0 0.4666+j*1.0482];

MC_Y(1:3,16:18) = (10^-6)*[5.6765 -0.6982 -1.8319; 0 5.3971 -1.1645; 0 0 5.9809]; %Configura��o 6
MC_Z(1:3,16:18) = [0.4576+j*1.0780 0.1535+j*0.3849 0.1560+j*0.5017; 0 0.4615+j*1.0651 0.1580+j*0.4236; 0 0 0.4666+j*1.0482];

MC_Y(1:3,25:27) = (10^-6)*[4.5193 0 0; 0 0 0; 0 0 0];                             %Configura��o 9
MC_Z(1:3,25:27) = [1.3292+j*1.3475 0 0; 0 0 0; 0 0 0];

MC_Y(1:3,28:30) = (10^-6)*[0 0 0; 0 4.5193 0; 0 0 0];                             %Configura��o 10
MC_Z(1:3,28:30) = [0 0 0; 0 1.3292+j*1.3475 0; 0 0 0];

MC_Y(1:3,31:33) = (10^-6)*[0 0 0; 0 0 0; 0 0 4.5193];                             %Configura��o 11
MC_Z(1:3,31:33) = [0 0 0; 0 0 0; 0 0 1.3292+j*1.3475];

%Dist�ncias e configura��o entre os n�s
Dist(78,1)=350*0.000189394;             %160-67 - Regulador 4 - 1
Dist(78,2)=2;             %Configura��o 2-6

Dist(88,1)=200*0.000189394;             %67-68 - 2
Dist(88,2)=3;             %Configura��o 3-9

Dist(89,1)=275*0.000189394;             %68-69 -3
Dist(89,2)=3;             %Configura��o 3-9

Dist(90,1)=325*0.000189394;             %69-70 - 4
Dist(90,2)=3;             %Configura��o 3-9

Dist(91,1)=275*0.000189394;             %70-71 - 5
Dist(91,2)=3;             %Configura��o 3-9

Dist(98,1)=250*0.000189394;             %67-97 - 6
Dist(98,2)=1;             %Configura��o 1-3

Dist(99,1)=275*0.000189394;             %97-98 - 7
Dist(99,2)=1;             %Configura��o 1-3

Dist(100,1)=550*0.000189394;             %98-99 - 8
Dist(100,2)=1;             %Configura��o 1-3

Dist(101,1)=300*0.000189394;             %99-100 - 9
Dist(101,2)=1;             %Configura��o 1-3

Dist(102,1)=800*0.000189394;             %100-450 - 10
Dist(102,2)=1;             %Configura��o 1-3

Dist(105,1)=250*0.000189394;             %97-197 - 11
Dist(105,2)=1;             %Configura��o 1-3

Dist(109,1)=250*0.000189394;             %197-101 - 12
Dist(109,2)=1;             %Configura��o 1-3

Dist(110,1)=225*0.000189394;             %101-102 - 13
Dist(110,2)=5;             %Configura��o 5-11

Dist(111,1)=325*0.000189394;             %102-103 - 14
Dist(111,2)=5;             %Configura��o 5-11

Dist(112,1)=700*0.000189394;             %103-104 - 15
Dist(112,2)=5;             %Configura��o 5-11

Dist(114,1)=275*0.000189394;             %101-105 - 16
Dist(114,2)=1;             %Configura��o 1-3

Dist(117,1)=225*0.000189394;             %105-106 - 17
Dist(117,2)=4;             %Configura��o 4-10

Dist(118,1)=575*0.000189394;             %106-107 - 18
Dist(118,2)=4;             %Configura��o 4-10

Dist(120,1)=325*0.000189394;             %105-108 - 19
Dist(120,2)=1;             %Configura��o 1-3

Dist(121,1)=1000*0.000189394;             %300-108 - 20
Dist(121,2)=1;             %Configura��o 1-3

Dist(123,1)=450*0.000189394;             %108-109 - 21
Dist(123,2)=3;             %Configura��o 3-9

Dist(124,1)=300*0.000189394;             %109-110 - 22
Dist(124,2)=3;             %Configura��o 3-9

Dist(125,1)=575*0.000189394;             %111-110 - 23
Dist(125,2)=3;             %Configura��o 3-9

Dist(126,1)=125*0.000189394;             %110-112 - 24
Dist(126,2)=3;             %Configura��o 3-9

Dist(127,1)=525*0.000189394;             %112-113 - 25
Dist(127,2)=3;             %Configura��o 3-9

Dist(128,1)=325*0.000189394;             %113-114 - 26
Dist(128,2)=3;             %Configura��o 3-9

%Pot�ncia nos n�s em estrela e delta (TC = Tipo de carga: 1-S, 2-I ou 3-Z)
%PPA = Percentual de Pot�ncia Ativa
PPA = 0             %N� 160 - 1
SN(1:3,78)=[0,0,0];
SN(4:6,78)=[0,0,0];

PPA = 0             %N� 67 - 2
SN(1:3,88)=[0,0,0];                                                                                                         
SN(4:6,88)=[0,0,0];

PPA = 0             %N� 68 - 3
SN(1:3,89)=[((20+j*10)*1e3),0,0];                                                                                                         
TC_Y(89) = 1;
SN(4:6,89)=[0,0,0];

PPA = 0             %N� 69 - 4
SN(1:3,90)=[((40+j*20)*1e3),0,0];                                                                                                          
TC_Y(90) = 1;
SN(4:6,90)=[0,0,0];

PPA = 0             %N� 70 - 5
SN(1:3,91)=[((20+j*10)*1e3),0,0]; 
TC_Y(91) = 1;
SN(4:6,91)=[0,0,0];

PPA = 0             %N� 71 - 6
SN(1:3,92)=[((40+j*20)*1e3),0,0];                                                                                                                  
TC_Y(2) = 1;
SN(4:6,92)=[0,0,0];

PPA = 0             %N� 97 - 7  
SN(1:3,99)=[0,0,0];                                                                                                                    
SN(4:6,99)=[0,0,0];

PPA = 0             %N� 98 - 8
SN(1:3,100)=[((40+1i*20)*1e3),0,0];                                                                     
TC_Y(100) = 1;
SN(4:6,100)=[0,0,0];

PPA = 0             %N� 99 - 9 
SN(1:3,101)=[0,((40+1i*20)*1e3),0];                                                                         
TC_Y(101) = 1;
SN(4:6,101)=[0,0,0];

PPA = 0             %N� 100 - 10
SN(1:3,102)=[0,0,(40+1i*20)*1e3];                                                                                   
TC_Y(102) = 3;
SN(4:6,102)=[0,0,0];

PPA = 0             %N� 450 - 11
SN(1:3,103)=[0,0,0];                                                                                                     
SN(4:6,103)=[0,0,0];

PPA = 0             %N� 197 - 12
SN(1:3,107)=[0,0,0];                                                                                         
SN(4:6,107)=[0,0,0]; 

PPA = 0             %N� 101 - 13
SN(1:3,110)=[0,0,0];                                                                                                
SN(4:6,110)=[0,0,0];

PPA = 0             %N� 102 - 14
SN(1:3,111)=[0,0,((20+1i*10)*1e3)];                                                                                  
TC_Y(111) = 1;
SN(4:6,111)=[0,0,0];

PPA = 0             %N� 103 - 15
SN(1:3,112)=[0,0,((40+1i*20)*1e3)];                                                                                                                
TC_Y(112) = 1;
SN(4:6,112)=[0,0,0];

PPA = 0             %104 - N� 16
SN(1:3,113)=[0,0,((40+1i*20)*1e3)];                                                                              
TC_Y(113) = 1;
SN(4:6,113)=[0,0,0];

PPA = 0            %N� 105 - 17
SN(1:3,117)=[0,0,0];                                                                   
SN(4:6,117)=[0,0,0]; 

PPA = 0             %N� 106 - 18
SN(1:3,118)=[0,((40+1i*20)*1e3),0];                                                                                                  
TC_Y(118) = 1;
SN(4:6,118)=[0,0,0];

PPA = 0             %N� 107 - 19 
SN(1:3,119)=[0,((40+1i*20)*1e3),0];                                              
TC_Y(119) = 1;
SN(4:6,119)=[0,0,0];

PPA = 0             %N� 108 - 20
SN(1:3,122)=[0,0,0];                                                                                                      
SN(4:6,122)=[0,0,0];

PPA = 0             %N� 300 - 21
SN(1:3,121)=[0,0,0];                                                                                          
SN(4:6,121)=[0,0,0];

PPA = 0             %N� 109 - 22
SN(1:3,124)=[((40+j*20)*1e3),0,0];                                                                                                          
TC_Y(124) = 1;
SN(4:6,124)=[0,0,0];

PPA = 0              %N� 110 - 23 
SN(1:3,126)=[0,0,0];                                                                                                         
SN(4:6,126)=[0,0,0];

PPA = 0             %N� 111 - 24
SN(1:3,125)=[((20+j*10)*1e3),0,0];                                                                                                         
TC_Y(125) = 1;
SN(4:6,125)=[0,0,0];

PPA = 0             %N� 112 - 25
SN(1:3,127)=[((20+j*10)*1e3),0,0]; 
TC_Y(127) = 2;
SN(4:6,127)=[0,0,0];

PPA = 0             %N� 113 - 26
SN(1:3,128)=[((40+j*20)*1e3),0,0];                                                                                                                     
TC_Y(128) = 3;
SN(4:6,128)=[0,0,0]; 

PPA = 0             %N� 114 - 27
SN(1:3,129)=[((20+j*10)*1e3),0,0]; 
TC_Y(129) = 1;
SN(4:6,129)=[0,0,0];

%Tens�o de refer�ncia
Vref = (4.16*10^3)/sqrt(3);

%Erro
Erro = 0.001;
