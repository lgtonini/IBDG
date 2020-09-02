

prompt = 'Número de configurações de linha: ';
NC = input(prompt);

for c = 1:NC
    prompt = sprintf('%dº Configuração Y (em S/km): ',c);
    MC_Y(1:3,3*c-2:3*c) = input(prompt);
    
    prompt = sprintf('%dº Configuração Z (em Ohms/km): ',c);
    MC_Z(1:3,3*c-2:3*c) = input(prompt);
end



%%
nc = 2;
MY1=zeros(T_l*3,T_c*3);
MZ1_CC=zeros(T_l*3,T_c*3);
cont=1;
for c = 1:NL
    prompt1 = sprintf('Distância entre %d,%d e %d,%d: ',amr(nc,1),amr(nc,2),amr(nc+1,1),amr(nc+1,2));
    Dist(c,1) = input(prompt1);
    prompt2 = sprintf('Qual a configuração da linha %d - %d,%d e %d,%d: ',c, amr(nc,1),amr(nc,2),amr(nc+1,1),amr(nc+1,2));
    a = input(prompt2);
    if a <= NC
        Dist(c,2) = a;
    else error('Número inválido')
    end
    nc = nc + 3;
    
    % Configuração x distância
   MZ(1:3,3*c-2:3*c) = MC_Z(1:3,3*a-2:3*a)*Dist(c,1);
   MY(1:3,3*c-2:3*c) = MC_Y(1:3,3*a-2:3*a)*Dist(c,1)/2;
   
   MZ_CC(1:3,3*c-2:3*c) = MC_Z(1:3,3*a-2:3*a)*Dist(c,1)/2;
   
   eixoX1=amr(cont+1,1);
   eixoY1=amr(cont+1,2);
   
   eixoX2=amr(cont+2,1);
   eixoY2=amr(cont+2,2);
   
   %Calculo de Y para cada nó
   
   for a2=1:3 
      for a1=1:3 
          
          if a1==1
              a11=3*eixoX1-2;
              a33=3*eixoX2-2;
              b11=1;
          end
          if a1==2
              a11=3*eixoX1-1;
              a33=3*eixoX2-1;
              b11=2;
          end
          if a1==3
              a11=3*eixoX1-0;
              a33=3*eixoX2-0;
              b11=3;
          end
          
          if a2==1
              a22=3*eixoY1-2;
              a44=3*eixoY2-2;
              b22=3*c-2;
          end
          if a2==2
              a22=3*eixoY1-1;
              a44=3*eixoY2-1;
              b22=3*c-1;
          end
          if a2==3
              a22=3*eixoY1-0;
              a44=3*eixoY2-0;
              b22=3*c-0;
          end

          MY1(a11,a22)=MY1(a11,a22)+MY(b11,b22);
          MY1(a33,a44)=MY1(a33,a44)+MY(b11,b22);
          
          MZ1_CC(a11,a22)=MZ1_CC(a11,a22)+MZ(b11,b22);
          MZ1_CC(a33,a44)=MZ1_CC(a33,a44)+MZ(b11,b22);
      end
   end
        
   cont=cont+3;
end

%%
cont=1;
for a1=1:T_l
    for a2=1:T_c
        if M(a1,a2)==1
          MY2(1:3,3*cont-2:3*cont)=MY1(3*a1-2:3*a1,3*a2-2:3*a2); %ordenar Y
          MZ2(1:3,3*cont-2:3*cont)=MZ1_CC(3*a1-2:3*a1,3*a2-2:3*a2); %ordenar Z
          %matriz posicoes
          posicoes(1,cont)= a1;
          posicoes(2,cont)= a2;
%             prompt = sprintf('Qual a porcentagem de potência ATIVA na GD, entre 0 e 1, neste nó: ');
%             PPA = input(prompt);
%             prompt = sprintf('Qual a porcentagem de potência REATIVA na GD, entre 0 e 0.44, neste nó: ');
%             PPR = input(prompt);
%             prompt = sprintf('Qual a potência no nó %d (%d,%d),cargas em estrela: ',cont, a1, a2);
%             PN_Y = input(prompt);
%             SN(1:3,cont) = (1-PPA).*real(PN_Y)+((1-PPR).*imag(PN_Y))*j;
%             SN_S_Y(1:3,cont) = abs(PPA.*real(PN_Y)+PPR.*imag(PN_Y)*j);
%             SN_P_Y(1:3,cont) = PPA.*real(PN_Y);          
%             prompt = sprintf('Qual a potência no nó %d (%d,%d),cargas em delta: ',cont, a1, a2);
%             PN_D = input(prompt);
%             SN(4:6,cont) = (1-PPA).*real(PN_D)+((1-PPR).*imag(PN_D))*j;
%             SN_S_D(1:3,cont) = abs(PPA.*real(PN_D)+PPR.*imag(PN_D)*j);
%             SN_P_D(1:3,cont) = PPA.*real(PN_D);
            
            prompt = sprintf('Qual a porcentagem de potência ATIVA na GD, entre 0 e 1, neste nó: ');
            PPA(1:3,cont) = input(prompt);
            prompt = sprintf('Qual a potência no nó %d (%d,%d),cargas em estrela: ',cont, a1, a2);
            SN_Y = input(prompt);
            SN_P_Y(1:3,cont) = real(SN_Y);          
            SN_Q_Y(1:3,cont) = imag(SN_Y);
            if SN_P_Y(1,cont)+SN_P_Y(2,cont)+SN_P_Y(3,cont)+SN_Q_Y(1,cont)+SN_Q_Y(2,cont)+SN_Q_Y(3,cont) ~= 0
                prompt = sprintf('Qual o tipo de carga no nó %d (%d,%d), em estrela: ',cont, a1, a2);
                TC_Y(cont) = input(prompt);
            else
                TC_Y(cont) = 1;
            end
            
            prompt = sprintf('Qual a potência no nó %d (%d,%d),cargas em delta: ',cont, a1, a2);
            SN_D = input(prompt);
            SN_P_D(1:3,cont) = real(SN_D);
            SN_Q_D(1:3,cont) = imag(SN_D);
            if SN_P_D(1,cont)+SN_P_D(2,cont)+SN_P_D(3,cont)+SN_Q_D(1,cont)+SN_Q_D(2,cont)+SN_Q_D(3,cont) ~= 0
                prompt = sprintf('Qual o tipo de carga no nó %d (%d,%d), em delta: ',cont, a1, a2);
                TC_D(cont) = input(prompt);
            else
                TC_D(cont) = 1;
            end

            SN_P_T(1:3,cont) = SN_P_Y(1:3,cont)+SN_P_D(1:3,cont);          
            SN_Q_T(1:3,cont) = SN_Q_Y(1:3,cont)+SN_Q_D(1:3,cont);   
            
            SN_S_T(1:3,cont) = SN_P_T(1:3,cont)+SN_Q_T(1:3,cont)*j;
            
            FP = 1;
            PPR(1:3,cont) = abs((PPA(1:3,cont).*SN_P_T(1:3,cont).*(1-FP)./(FP.*SN_Q_T(1:3,cont)*j)));
            
            PPR(isnan(PPR))=0;
                     
            if PPR(1,cont)<=0 && PPR(1,cont)>=0.44 && PPR(2,cont)<=0 && PPR(2,cont)>=0.44 && PPR(3,cont)<=0 && PPR(3,cont)>=0.44 
               disp('Saiu da faixa');
            end
            
            SN(1:3,cont) = (1-PPA(1:3,cont)).*real(SN_Y)'+(1-PPR(1:3,cont)).*imag(SN_Y)'*j;
            SN(4:6,cont) = (1-PPA(1:3,cont)).*real(SN_D)'+(1-PPR(1:3,cont)).*imag(SN_D)'*j;
            
            
            if abs(SN(1:3,cont)+SN(4:6,cont)) == [0;0;0]
                SN_FP(1:3,cont) = 0;
            else
                SN_FP(1:3,cont) = ((PPA(1:3,cont)).*SN_P_T(1:3,cont))./(abs((PPA(1:3,cont)).*SN_P_T(1:3,cont)+(PPR(1:3,cont)).*SN_Q_T(1:3,cont)));
            
            end
            
            if abs(PPA(1:3,cont)+PPR(1:3,cont)) == [0;0;0]
                SN_FP(1:3,cont) = 0;
            end
            
            SN_FP(isnan(SN_FP))=0;
            
            cont=cont+1;  
        end
               
    end
end
%%
%Entrada de dados do transformador
cont=1;
if T_Traf >= 1
    prompt = 'Quantas configurações de transformador? ';
    NCT = input(prompt);
end
%CT (configuração do transformador)
% 1- Tipo (T_T)
% 2- Potência (PPT)
% 3- Tensão de entrada (VTH)
% 4- Tensão de saída (VTL)
% 5- Impedâcia (ZTT)
% 6- Impedância de base (ZBASE)
% 7- Relação de transformação (NT)

%CTM (Configuração do transformador em matriz)
% 1:3- Impedância real (ZTT)
% 4:6- aTT
% 7:9- bTT
% 10:12- dTT
% 13:15- A_TT
% 16:18- B_TT
% 19:21- AV_TT

for a1=1:T_l
    for a2=1:T_c
        if M(a1,a2)==3
            prompt = sprintf('Qual o tipo do transformador %d (1-DYg e 2-YgYg)? ', cont);
            CT(1,cont) = input(prompt);
            prompt = sprintf('Qual a potência do transformador %d [VA]? ', cont);
            CT(2,cont) = input(prompt);
            prompt = sprintf('Qual a tensão de entrada do transformador %d [V]? ', cont);
            CT(3,cont) = input(prompt);
            prompt = sprintf('Qual a tensão de saída do transformador %d [V]? ', cont);
            CT(4,cont) = input(prompt);
            prompt = sprintf('Qual a impedancia do transformador %d [R+j*X]? ', cont);
            CT(5,cont) = input(prompt);
            
            CT(6,cont)=(CT(4,cont)^2)/CT(2,cont);
            
            CTM(1:3,3*cont-2:3*cont)=CT(5,cont)*CT(6,cont)*[1 0 0; 0 1 0; 0 0 1]; 
            
            if CT(1,cont)==1 % se for DYg
                CT(7,cont) = CT(3,cont)/(sqrt(3)*CT(4,cont));
                
                CTM(4:6,3*cont-2:3*cont) = (-CT(7,cont)/3)*[0 2 1; 1 0 2; 2 1 0];
                CTM(7:9,3*cont-2:3*cont) = (-CT(7,cont)/3)*CTM(1:3,1:3)*[0 2 1; 1 0 2; 2 1 0];
                CTM(10:12,3*cont-2:3*cont) = (1/CT(7,cont))*[1 -1 0; 0 1 -1; -1 0 1];
                CTM(13:15,3*cont-2:3*cont) = (1/CT(7,cont))*[1 0 -1; -1 1 0; 0 -1 1];
                CTM(16:18,3*cont-2:3*cont) = CTM(1:3,1:3);
                CTM(19:21,3*cont-2:3*cont) = [0 -1/CT(7,cont) 0; 0 0 -1/CT(7,cont); -1/CT(7,cont) 0 0];
%                 CTM(19:21,3*cont-2:3*cont) = [0 -CT(7,cont) 0; 0 0 -CT(7,cont); -CT(7,cont) 0 0];
            end
            
            if CT(1,cont)==2 % se for YgYg
                CT(7,cont) = CT(3,cont)/CT(4,cont);
                
                CTM(4:6,3*cont-2:3*cont) = (CT(7,cont))*[1 0 0; 0 1 0; 0 0 1];
                CTM(7:9,3*cont-2:3*cont) = CTM(4:6,1:3)*CTM(1:3,1:3);
                CTM(10:12,3*cont-2:3*cont) = (1/CT(7,cont))*[1 0 0; 0 1 0; 0 0 1];
                CTM(13:15,3*cont-2:3*cont) = [1/CT(7,cont) 0 0; 0 1/CT(7,cont) 0; 0 0 1/CT(7,cont)];
                CTM(16:18,3*cont-2:3*cont) = CTM(1:3,1:3);
                CTM(19:21,3*cont-2:3*cont) = CTM(13:15,1:3);
            end
            
            cont=cont+1;
        end
    end
end

%Entrada de dados do regulador

% CR (Configuração do regulador)
% 1- Relaçao de trasnformação do regulador (NTR)
% 2- Impedância do regulador (ZTR)
cont=1;

for a1=1:T_l
    for a2=1:T_c
        if M(a1,a2)==4
            prompt = sprintf('Qual a relação de transformação do regulador %d? ', cont);
            relacao = input(prompt);
            CR(1,cont) = 1+(1/((relacao)));
            prompt = sprintf('Qual a impedância do regulador %d? ', cont);
            CR(2,cont) = input(prompt);
            cont=cont+1;
        end
    end
end
%%
prompt = sprintf('Qual a tensão de referência: ');
Vref = input(prompt);

%%
for f = 1:NN
    VN(3*f-2:3*f,1) = [Vref*1;Vref*(-0.5-j*0.8660254);Vref*(-0.5+j*0.8660254)];
end

prompt = 'Qual o erro: ';
Erro = input(prompt);

prompt = 'Qual a configuração de Falta: ';
ZFalta = input(prompt);