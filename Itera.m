tic


M3=M; % M deve conter a matriz com 0,1,2,3 e 4

size_M=size(M);
size_M1=size_M(1,1);
size_M2=size_M(1,2);


for x=1:size_M1
    for y=1:size_M2
        
        if M(x,y) ~=1 && M(x,y) ~=0
            M(x,y)=2;
        end
    end
end

conti1=1;
conti2=1;
conti3=1;

for x=1:size_M1
    for y=1:size_M2
        
        if M(x,y) ==2
            if M3(x,y) ==3
                
                MZ(1:3,3*conti1-2:3*conti1)=MZ(1:3,3*conti1-2:3*conti1)+CTM(1:3,3*conti2-2:3*conti2);
                
                conti2=conti2+1;
            end
            
            if M3(x,y) ==4
                
                MZ(1:3,3*conti1-2:3*conti1)=MZ(1:3,3*conti1-2:3*conti1)+CR(2,conti3)*[1 0 0; 0 1 0; 0 0 1];
                
                conti3=conti3+1;
            end
            
            
            conti1=conti1+1;
        end
        
    end

end

Modifica = 1;

while Modifica == 1;
    
    MOD = 0;

%Vref=4160/sqrt(3)*[1.069 1.05555 1.08265];
%Erro=1e-3;
%Cálculo da tensão nodal a partir da tensão de referência
for f = 1:NN
    
    VN(3*f-2:3*f,1) = [Vref*1;Vref*(-0.5-j*0.8660254);Vref*(-0.5+j*0.8660254)];
   % VN(3*f-2:3*f,1) = [Vref(1,1)*1;Vref(1,2)*(-0.5-j*0.8660254);Vref(1,3)*(-0.5+j*0.8660254)];
end

erro=10; % Atribuição de um valor inicial para o erro que seja maior que o critério de parada
V_nos2(:,1)=VN;
ite=2; %Variável que determina o número de iterações

while erro>Erro

    %Cálculo da corrente nodal
    for f = 1:NN
        
        IN(:,f) = -[conj((1)*j*MY2(1:3,3*f-2:3*f))*VN(3*f-2:3*f,1)];      
        
        if TC_D(f) == 1
            IN(:,f) = IN(:,f)+[(conj(SN(4,f)/(VN(3*f-2,1)-VN(3*f-1,1)))-conj(SN(6,f)/(VN(3*f,1)-VN(3*f-2,1))));(conj(SN(5,f)/(VN(3*f-1,1)-VN(3*f,1)))-conj(SN(4,f)/(VN(3*f-2,1)-VN(3*f-1,1))));(conj(SN(6,f)/(VN(3*f,1)-VN(3*f-2,1)))-conj(SN(5,f)/(VN(3*f-1,1)-VN(3*f,1))))];
        end
        
        if TC_D(f) == 2
            IN(:,f) = IN(:,f)+[conj(SN(4,f)/(Vref*(cosd(30)+j*sind(30))))-conj(SN(6,f)/(Vref*(cosd(150)+j*sind(150))));conj(SN(5,f)/(Vref*(cosd(-90)+j*sind(-90))))-conj(SN(4,f)/(Vref*(cosd(30)+j*sind(30))));conj(SN(6,f)/(Vref*(cosd(150)+j*sind(150))))-conj(SN(4,f)/(Vref*(cosd(-90)+j*sind(-90))))];
        end
        
        if TC_D(f) == 3
            IN(:,f) = IN(:,f)+[(conj(SN(4,f)*conj(VN(3*f-2,1)-VN(3*f-1,1)))-conj(SN(6,f)*conj(VN(3*f,1)-VN(3*f-2,1))))/(Vref*conj(Vref));(conj(SN(5,f)*conj(VN(3*f-1,1)-VN(3*f,1)))-conj(SN(4,f)*conj(VN(3*f-2,1)-VN(3*f-1,1))))/(Vref*conj(Vref));(conj(SN(6,f)*conj(VN(3*f,1)-VN(3*f-2,1)))-conj(SN(5,f)*conj(VN(3*f-1,1)-VN(3*f,1))))/(Vref*conj(Vref))];
        end
        
        if TC_Y(f) == 1 
            IN(:,f) = IN(:,f)+[conj(SN(1,f)/VN(3*f-2,1));conj(SN(2,f)/VN(3*f-1,1));conj(SN(3,f)/VN(3*f,1))];
        end
        
        if TC_Y(f) == 2
            IN(:,f) = IN(:,f)-[conj(SN(1,f)/(Vref*(cosd(0)+j*sind(0))));conj(SN(2,f)/(Vref*(cosd(-120)+j*sind(-120))));conj(SN(3,f)/(Vref*(cosd(120)+j*sind(120))))];
        end
        
        if TC_Y(f) == 3
            IN(:,f) = IN(:,f)+[conj(SN(1,f)*VN(3*f-2,1)/(Vref*conj(Vref)));conj(SN(2,f)*VN(3*f-1,1)/(Vref*conj(Vref)));conj(SN(3,f)*VN(3*f,1)/(Vref*conj(Vref)))];
        end
        
         %IN(:,f)=[(conj(SN(4,f)      /(VN(3*f-2,1)-VN(3*f-1,1)))            -conj(SN(6,f)      /(VN(3*f,1)         -VN(3*f-2,1))))      +conj(SN(1,f)        /VN(3*f-2,1))      ;(conj(SN(5,f)      /(VN(3*f-1,1)       -VN(3*f,1)))      -conj(SN(4,f)      /(VN(3*f-2,1)      -VN(3*f-1,1))))      +conj(SN(2,f)        /VN(3*f-1,1))      ;(conj(SN(6,f)       /(VN(3*f,1)       -VN(3*f-2,1)))      -conj(SN(5,f)      /(VN(3*f-1,1)      -VN(3*f,1))))        +conj(SN(4,f)        /VN(3*f,1))]        -[conj((1)   *j*MY2(1:3,3*f-2:3*f))*VN(3*f-2:3*f,1) ];

 %I_824(:,:,cont)=[(conj(S_824_a_delta/(V_824(1,:,cont-1)-V_824(2,:,cont-1)))-conj(S_824_c_delta/(V_824(3,:,cont-1)-V_824(1,:,cont-1))))+conj(S_824_a_estrela/V_824(1,:,cont-1));(conj(S_824_b_delta/(V_824(2,:,cont-1)-V_824(3,:,cont-1)))-conj(S_824_a_delta/(V_824(1,:,cont-1)-V_824(2,:,cont-1))))+conj(S_824_b_estrela/V_824(2,:,cont-1));(conj(S_824_c_delta/(V_824(3,:,cont-1)-V_824(1,:,cont-1)))-conj(S_824_b_delta/(V_824(2,:,cont-1)-V_824(3,:,cont-1))))+conj(S_824_c_estrela/V_824(3,:,cont-1))]-[conj((1e-6)*j*Y_824)*V_824(:,:,cont-1) ];

        %IN(:,f)=[conj(SN(1:3,f)./VN(3*f-2:3*f,1))+conj([SN(4,f)./(VN(3*f-2,1)-VN(3*f-1,1));SN(5,f)./(VN(3*f-1,1)-VN(3*f,1));SN(6,f)./(VN(3*f,1)-VN(3*f-2,1))])-[conj((1)*j*MY2(1:3,3*f-2:3*f))*VN(3*f-2:3*f,1)]];
    end
    
%%
b=1;
M1=M;
amr1=amr;
    
%Tirar linha da matriz "amr" que contem a localização das linhas do sistema
for a1=1:NL
    amr1(2*a1-1,:)=[];
end

%Identifica os nós da extremidade (varrendo a matriz "amr1' e verificando
%se aparece mais de uma vez)
%Dentro dessa função a primeira coluna das matrizes "nos" e "ramo" é
%formada com os nós e ramos das extremidades, respectivamente

%Essas matrizes organizam de maneira hierarquizada os nós e linhas do sistema
for g = 1:3*NL-NL
     a=0;
     
     if amr1(g,1)==1
        a=2;
     end
        
    for h = 1:3*NL-NL
        if (amr1(g,1) == amr1(h,1) && amr1(g,2) == amr1(h,2) )
            a=a+1;
        end
    end 

    if a<2 
        M1(amr1(g,1),amr1(g,2))=3;
        
        nos(2*b-1,1)=amr1(g,1);
        nos(2*b,1)=amr1(g,2);
        
        l=nos(2*b-1,1);
        c=nos(2*b,1);
        
                if l ~= 1 
                    if M(l-1,c) == 2
                        ramo(2*b-1,1)=l-1;
                        ramo(2*b,1)=c;
                        
                        M1(l-1,c)=4;
                    end
                end
                if l ~= T_l
                    if M(l+1,c) == 2
                        ramo(2*b-1,1)=l+1;
                        ramo(2*b,1)=c;
                        M1(l+1,c)=4;
                    end
                end
                if c ~= 1 
                    if M(l,c-1) == 2
                        ramo(2*b-1,1)=l;
                        ramo(2*b,1)=c-1;
                        M1(l,c-1)=4;
                    end
                end
                if c~= T_c
                    if M(l,c+1) == 2
                        ramo(2*b-1,1)=l;
                        ramo(2*b,1)=c+1;
                        M1(l,c+ 1)=4;
                    end
                end     
        b=b+1;
        
    end
end

%NNE corresponde ao número de nós na extremidade
NNE1=b-1;
NNE=b-1;

%inicializa as matrizes "nos" e "ramo". (a partir da segunda coluna, já que
%a segunda já foi preenchida)
for z=2:NL
    nos(:,z)=zeros(2*NNE,1);
    ramo(:,z)=zeros(2*NNE,1);
end

%Identifica a barra de referência (pedir para colocar a barra de referencia
%na primeira linha)
for a1=1:T_c
    if M(1,a1)==1
        ref_c=a1;
    end
end

b=1;
crit_parada=1e3;
parada=0;


% A função "while" para no momento em que chega na barra de referência
% Nessa função as matrizes "nos" e "ramos" são formadas
    while parada~=crit_parada

for a1=1:NNE
        l=ramo(2*a1-1,b);
        c=ramo(2*a1,b);

            if (l-1)==1 && c==ref_c
                        nos(2*a1-1,b+1)=1;
                        nos(2*a1,b+1)=ref_c;
                           parada=crit_parada;
            else
    
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if l~=0 && c~=0                      
    
                if l ~= 1 
                    if M1(l-1,c) == 1 && ( (l-1)~=nos(2*a1-1,b) || (c)~=nos(2*a1,b) )

                        %testar se os ramos adjacentes já foram computados
                        %(caso tenha sido, manter o processo, caso
                        %contrário armazenar 0)  
                        a5=0;
               if l-1 ~= 1    %cima      
                        if M1(l-1-1,c)==2 || M1(l-1-1,c)==b+4 %|| M(l-1-1,c)==0
                            a5=a5+1;
                        end
               end
               if l-1 ~= T_l %baixo         
                        if M1(l-1+1,c)==2 || M1(l-1+1,c)==b+4 %|| M(l-1+1,c)==0
                            a5=a5+1;
                        end
               end
               if c ~= 1 %esquerda             
                        if M1(l-1,c-1)==2 || M1(l-1,c-1)==b+4  %|| M(l-1,c-1)==0
                            a5=a5+1;
                        end
               end
               if c~= T_c %direita
                        if M1(l-1,c+1)==2 || M1(l-1,c+1)==b+4 %|| M(l-1,c+1)==0
                            a5=a5+1;
                        end
               end
               %verificar se já avaliou todos ao redor
                         if a5==1
                        nos(2*a1-1,b+1)=l-1;
                        nos(2*a1,b+1)=c;
                        M1(l-1,c)=3;
                         else
                             nos(2*a1-1,b+1)=0;
                             nos(2*a1,b+1)=0;
                             %NNE=NNE-1;
                         end
                    end
                end
   
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
                if l ~= T_l
                    if M1(l+1,c) == 1 && ( (l+1)~=nos(2*a1-1,b) || (c)~=nos(2*a1,b) )
                        
                        if (l+1)==1 && c==ref_c
                           parada=crit_parada;
                        end
                        
                        a5=0;
                if l+1 ~= 1    %cima      
                        if M1(l+1-1,c)==2 || M1(l+1-1,c)==b+4  %|| M(l-1-1,c)==0
                            a5=a5+1;
                        end
               end
               if l+1 ~= T_l %baixo         
                        if M1(l+1+1,c)==2 || M1(l+1+1,c)==b+4 %|| M(l-1+1,c)==0
                            a5=a5+1;
                        end
               end
               if c ~= 1 %esquerda             
                        if M1(l+1,c-1)==2 || M1(l+1,c-1)==b+4  %|| M(l-1,c-1)==0
                            a5=a5+1;
                        end
               end
               if c~= T_c %direita
                        if M1(l+1,c+1)==2 || M1(l+1,c+1)==b+4 %|| M(l-1,c+1)==0
                            a5=a5+1;
                        end
               end
               
                         if a5==1
                        nos(2*a1-1,b+1)=l+1;
                        nos(2*a1,b+1)=c;
                        M1(l+1,c)=3;
                         else
                             nos(2*a1-1,b+1)=0;
                             nos(2*a1,b+1)=0;
                             %NNE=NNE-1;
                         end                        
                    end
                end
                
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%               
                
                
                
                if c ~= 1 
                    if M1(l,c-1) == 1 && ( (l)~=nos(2*a1-1,b) || (c-1)~=nos(2*a1,b) )
                        
                        if (l)==1 && (c-1)==ref_c
                           parada=crit_parada;
                        end
                        
                        a5=0;
              if l ~= 1    %cima      
                        if M1(l-1,c-1)==2 || M1(l-1,c-1)==b+4 %|| M(l-1-1,c)==0
                            a5=a5+1;
                        end
               end
               if l ~= T_l %baixo         
                        if M1(l+1,c-1)==2 || M1(l+1,c-1)==b+4 %|| M(l-1+1,c)==0
                            a5=a5+1;
                        end
               end
               if c-1 ~= 1 %esquerda             
                        if M1(l,c-1-1)==2 || M1(l,c-1-1)==b+4  %|| M(l-1,c-1)==0
                            a5=a5+1;
                        end
               end
               if c-1~= T_c %direita
                        if M1(l,c-1+1)==2 || M1(l,c-1+1)==b+4 %|| M(l-1,c+1)==0
                            a5=a5+1;
                        end
               end
                        if a5==1
                        nos(2*a1-1,b+1)=l;
                        nos(2*a1,b+1)=c-1;
                        M1(l,c-1)=3;
                         else
                             nos(2*a1-1,b+1)=0;
                             nos(2*a1,b+1)=0;
                             %NNE=NNE-1;
                        end  
    
                    end
                end
                
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
                
                
                if c~= T_c
                    if M1(l,c+1) == 1 && ( (l)~=nos(2*a1-1,b) || (c+1)~=nos(2*a1,b) )
                        
                        if (l)==1 && (c+1)==ref_c
                           parada=crit_parada;
                        end
                        
                        a5=0;
              if l ~= 1    %cima      
                        if M1(l-1,c+1)==2 || M1(l-1,c+1)==b+4  %|| M(l-1-1,c)==0
                            a5=a5+1;
                        end
               end
               if l ~= T_l %baixo         
                        if M1(l+1,c+1)==2 || M1(l+1,c+1)==b+4 %|| M(l-1+1,c)==0
                            a5=a5+1;
                        end
               end
               if c+1 ~= 1 %esquerda             
                        if M1(l,c+1-1)==2 || M1(l,c+1-1)==b+4  %|| M(l-1,c-1)==0
                            a5=a5+1;
                        end
               end
               if c+1~= T_c %direita
                        if M1(l,c+1+1)==2 || M1(l,c+1+1)==b+4 %|| M(l-1,c+1)==0
                            a5=a5+1;
                        end
               end
                        if a5==1
                        nos(2*a1-1,b+1)=l;
                        nos(2*a1,b+1)=c+1;
                        M1(l,c+1)=3;
                         else
                             nos(2*a1-1,b+1)=0;
                             nos(2*a1,b+1)=0;
                             %NNE=NNE-1;
                        end 
              
                    end
                end  

end
                             %if nos(2*a1-1,b+1)==0 && nos(2*a1,b+1)==0
                             %NNE=	-1
                             %end

    
l=nos(2*a1-1,b+1);
c=nos(2*a1,b+1);

if l~=0 && c~=0
                if l ~= 1 
                    if M1(l-1,c) == 2 && ( (l-1)~=ramo(2*a1-1,b) || (c)~=ramo(2*a1,b) ) && M1(l-1,c)~=4
                        ramo(2*a1-1,b+1)=l-1;
                        ramo(2*a1,b+1)=c;
                        M1(l-1,c)=b+4;
                    end
                end
                if l ~= T_l
                    if M1(l+1,c) == 2 && ( (l+1)~=ramo(2*a1-1,b) || (c)~=ramo(2*a1,b) ) && M1(l+1,c)~=4
                        ramo(2*a1-1,b+1)=l+1;
                        ramo(2*a1,b+1)=c;
                        M1(l+1,c)=b+4;
                    end
                end
                if c ~= 1 
                    if M1(l,c-1) == 2 && ( (l)~=ramo(2*a1-1,b) || (c-1)~=ramo(2*a1,b) ) && M1(l,c-1)~=4
                        ramo(2*a1-1,b+1)=l;
                        ramo(2*a1,b+1)=c-1;
                        M1(l,c-1)=b+4;
                    end
                end
                if c~= T_c
                    if M1(l,c+1) == 2 && ( (l)~=ramo(2*a1-1,b) || (c+1)~=ramo(2*a1,b) ) && M1(l,c+1)~=4
                        ramo(2*a1-1,b+1)=l;
                        ramo(2*a1,b+1)=c+1;
                        M1(l,c+1)=b+4; 
                    end
                end 
%%}

end
            end
% NO FINAL TESTAR SE O PRÓXIMO NÓ SÓ TEM UMA LIGAÇÃO

end

b=b+1;
 end

 %Simplificar a matriz "nos" e "ramo" (eliminando as linhas e colunas com
 %zeros)
q=0;
T1 = size(nos);
T_l1=T1(1,1);
T_c1=T1(1,2);

T2 = size(ramo);
T_l2=T2(1,1);
T_c2=T2(1,2);
 
 while q == 0 
    for r = 1:T_l1
        if nos(r,T_c1) ~= 0 
            q = q+1;
         else q = q+0;
        end
    end
    if q == 0
        T_c1 = T_c1-1;
        nos = nos(:,1:T_c1);
        end
 end 

 q=0;
  while q == 0 
    for r = 1:T_l2
        if ramo(r,T_c2) ~= 0 
            q = q+1;
         else q = q+0;
        end
    end
    if q == 0
        T_c2 = T_c2-1;
        ramo = ramo(:,1:T_c2);
        end
  end


%% 
% Forma a matriz I_ramos (I_ij)= (I_j)+ Somatório (I_jM), sendo M os nós a
% jusante do nó j

size1=size(ramo);
l_r=size1(1,1);
c_r=size1(1,2);
cont=1;

%primeira coluna
for x=1:2:l_r
    I_ramos(2*cont-1,1)=ramo(2*cont-1,1);
    I_ramos(2*cont,1)=ramo(2*cont,1);
    
    I_ramos(2*cont-1,2)=nos(2*cont-1,1);
    I_ramos(2*cont,2)=nos(2*cont,1);
    cont=cont+1;
end


for y=2:c_r
    cont1=1;
    for x=1:2:l_r
        cont2=3;
        
        l=ramo(2*cont1-1,y);
        c=ramo(2*cont1,y);
        
   if ramo(2*cont1-1,y)~=0 && ramo(2*cont1,y)~=0
        
        I_ramos(2*cont-1,1)=ramo(2*cont1-1,y);
        I_ramos(2*cont,1)=ramo(2*cont1,y);
    
        I_ramos(2*cont-1,2)=nos(2*cont1-1,y);
        I_ramos(2*cont,2)=nos(2*cont1,y);
        
        
        for y1=1:c_r
            cont4=1;
            
            if y1<y  
        for x1=1:2:l_r
              if c <ref_c
                  if c>ramo(2*cont4,y1) && abs(c-ramo(2*cont4,y1))<=2 && abs(l-ramo(2*cont4-1,y1))<=2 && (c~=ramo(2*cont4,y1)||l~=ramo(2*cont4-1,y1))&& (abs(c-ramo(2*cont4,y1))<2 || abs(l-ramo(2*cont4-1,y1))<2)
                      I_ramos(2*cont-1,cont2)=ramo(2*cont4-1,y1);
                      I_ramos(2*cont,cont2)=ramo(2*cont4,y1);
                      cont2=cont2+1;
                  end
              end
              if c >ref_c
                    if c<ramo(2*cont4,y1) && abs(c-ramo(2*cont4,y1))<=2 && abs(l-ramo(2*cont4-1,y1))<=2 && (c~=ramo(2*cont4,y1)||l~=ramo(2*cont4-1,y1))&& (abs(c-ramo(2*cont4,y1))<2 || abs(l-ramo(2*cont4-1,y1))<2)
                      I_ramos(2*cont-1,cont2)=ramo(2*cont4-1,y1);
                      I_ramos(2*cont,cont2)=ramo(2*cont4,y1);
                      cont2=cont2+1;
                  end
              end
              if c ==ref_c
                  if l<ramo(2*cont4-1,y1) && abs(c-ramo(2*cont4,y1))<=2 && abs(l-ramo(2*cont4-1,y1))<=2 && (c~=ramo(2*cont4,y1)||l~=ramo(2*cont4-1,y1))&& (abs(c-ramo(2*cont4,y1))<2 || abs(l-ramo(2*cont4-1,y1))<2)
                      I_ramos(2*cont-1,cont2)=ramo(2*cont4-1,y1);
                      I_ramos(2*cont,cont2)=ramo(2*cont4,y1);
                      cont2=cont2+1;
                  end
              end

cont4=cont4+1;
        end  
        end
        end
    cont=cont+1;
   end
   cont1=cont1+1;
    end
end
%%
size13=size(I_ramos);
size1=size13(1,1);
size2=size13(1,2);

%%
%Encontra o valor das correntes dos ramos presente na matriz "I_ramos1"
for x=1:size1/2
    for y=1:NN
    
        if I_ramos(2*x-1,2)==posicoes(1,y) && I_ramos(2*x,2)==posicoes(2,y)
            
            I_ramos1(3*x-2:3*x,1)=IN(:,y); % 2ª parcela da eq: (I_ij)=(I_j)+Somatório(I_jM) --> (I_j)
                    
        end
    end
end

for x= (NNE1+1):size1/2 
    for y= 1:size1/2
        for z= 3:size2
     if I_ramos(2*x-1,z)==I_ramos(2*y-1,1) && I_ramos(2*x,z)==I_ramos(2*y,1)

            I_ramos1(3*x-2:3*x,1)=I_ramos1(3*x-2:3*x,1)+I_ramos1(3*y-2:3*y,1); % 2ª parcela da eq: (I_ij)=(I_j)+Somatório(I_jM) --> Somatório(I_jM)
                    
        end
    end
end

end
%% 
% Mesmo procedimento para o processo realizado para a corrente; nesse caso
% foi obtida a matriz "V_nos". Vj=Vi-Zij*Iij 
% 1ª coluna: Vj; 2ª coluna: linha(ij); 3ª coluna: Vi

size_nos=size(nos);
l_n=size_nos(1,1);
c_n=size_nos(1,2);

cont2=1;
for y = (c_n-1):(-1):1
    cont1=1;
    for x = 1:2:l_n

       l2=nos(2*cont1-1,y);
       c2=nos(2*cont1,y); 
       
       if nos(2*cont1-1,y)~=0 && nos(2*cont1,y)~=0
       
           V_nos(2*cont2-1,1)=l2;
           V_nos(2*cont2,1)=c2;
           
           y11=2;

           for y1=(y+1):c_n % ou y+0

               cont3=1;
               for x1= 1:2:l_n
                   
                   l1=nos(2*cont3-1,y1);
                   c1=nos(2*cont3,y1);
                   
                   if ((l1~=0 && c1~=0) && l2>=  l1 && abs(l2-l1)<=2 && abs(c2-c1)<=2 && (l1==l2 || c1==c2) ) && ( l2~=l1 || c2~=c1 ) && y11~=1 && M((l1+l2)/2,(c1+c2)/2)==2;
                   
                       V_nos(2*cont2-1,3)=l1;
                       V_nos(2*cont2,3)=c1;
                       y11=1;
                       
                   end
                   cont3=cont3+1;
               end
           end
        cont2=cont2+1;   
       end
        cont1=cont1+1;
    end
end

%Encontrar os ramos entre os nós
for x=1:2*(NN-1)
    V_nos(x,2)=(V_nos(x,1)+V_nos(x,3))/2;
end

%% identificar posições dos ramos (para relacionar com as impedâncias)
% A matriz "posicoes_r" ordena os ramos da mesma forma que as impedâncias
% estão ordenadas na matriz "MZ"

size_feeder=size(M);
l_feeder=size_feeder(1,1);
c_feeder=size_feeder(1,2);
cont=1;

for x = 1:l_feeder
    for y = 1:c_feeder

if M(x,y)==2
    posicoes_r(1,cont)=x;
    posicoes_r(2,cont)=y;
    cont=cont+1;
end

    end
end

%% Atribuir a tensão de referência na matriz "V_nos1"
% V_nos1=[Vref(1,1)*1 ;Vref(1,2)*(-0.5-j*0.8660254) ;Vref(1,3)*(-0.5+j*0.8660254)];    
V_nos1=1.05*[Vref*1 ;Vref*(-0.5-j*0.8660254) ;Vref*(-0.5+j*0.8660254)];

%%
% Calculo das outras tensões nodais: "V_nos1"

cont_reg=1;
cont_trafo=1;

for aux = 1:(NN-1)
%for aux = 11:11
    r1= V_nos(2*aux-1,2);
    r2=V_nos(2*aux,2);
    
    n1= V_nos(2*aux-1,3);
    n2=V_nos(2*aux,3);
    
    
    
    %Calculo da tensão dos nós diretamente conectados à barra de referência
    if n1==1 && n2==ref_c
        for aux1=1:NL %varrer as posições das impedâncias
            if r1==posicoes_r(1,aux1) && r2==posicoes_r(2,aux1)
                for aux2=1:NL %varrer as posições das correntes dos ramos
                     if r1==I_ramos(2*aux2-1,1) && r2==I_ramos(2*aux2,1)
              
                         if M3(r1,r2)==4
                      V_nos1(3*(aux+1)-2:3*(aux+1),1)=CR(1,cont_reg)*[1 0 0; 0 1 0; 0 0 1]*(V_nos1(3*1-2:3*1,1))- MZ(1:3,3*aux1-2:3*aux1)*I_ramos1(3*aux2-2:3*aux2,1);
                          cont_reg=cont_reg+1;
                         else
                             
                      V_nos1(3*(aux+1)-2:3*(aux+1),1)=V_nos1(3*1-2:3*1,1)- MZ(1:3,3*aux1-2:3*aux1)*I_ramos1(3*aux2-2:3*aux2,1);
 
                         end

                   end
                end
            end
        end
    else
     %Calculo da tensão dos outros nós
        for aux3=1:(NN-1)
             if n1== V_nos(2*aux3-1,1) && n2==V_nos(2*aux3,1)

                NLL=size(posicoes_r);
                NLL1=NLL(1,2);
                  for aux1=1:NLL1 %varrer as posições das impedâncias

                       if r1==posicoes_r(1,aux1) && r2==posicoes_r(2,aux1)

                          for aux2=1:NL %varrer as posições das correntes dos ramos
                               if r1==I_ramos(2*aux2-1,1) && r2==I_ramos(2*aux2,1)
                                   
                                   
                                   if M3(r1,r2)==4
                          V_nos1(3*(aux+1)-2:3*(aux+1),1)=CR(1,cont_reg)*[1 0 0; 0 1 0; 0 0 1]*(V_nos1(3*(aux3+1)-2:3*(aux3+1),1))- MZ(1:3,3*aux1-2:3*aux1)*I_ramos1(3*aux2-2:3*aux2,1);
                          cont_reg=cont_reg+1;
                                  else
                             
                          V_nos1(3*(aux+1)-2:3*(aux+1),1)=V_nos1(3*(aux3+1)-2:3*(aux3+1),1)- MZ(1:3,3*aux1-2:3*aux1)*I_ramos1(3*aux2-2:3*aux2,1);
 
                                   end

                               end
                          end
                       end
                  end
             end
        end
    end
    
end

 VN=V_nos1;
 V_nos2(:,ite)=V_nos1;
 I_ramos2(:,ite)=I_ramos1;
 
 %Cálculo do erro
 for e=1:NN*3
    er(e)= abs(abs(V_nos2(e,ite))-abs(V_nos2(e,ite-1)));
 end
 erro(ite-1)=sum(er);

ite=ite+1;

end
abs(V_nos1)*sqrt(3)/4160;

%% Análise de curto-circuito
ramos_cc(1,1)=posicoes(1,1);
ramos_cc(2,1)=posicoes(2,1);

for z=2:NN
x1=posicoes(1,z);
y1=posicoes(2,z);    

ramos_cc(2*z-1,1)=x1;
ramos_cc(2*z,1)=y1;

%procurar na matriz V_nos
for z1=1:(NN-1)
    
if x1==V_nos(2*z1-1,1) && y1==V_nos(2*z1,1)
    ramos_cc(2*z-1,2)=V_nos(2*z1-1,2);
    ramos_cc(2*z,2)=V_nos(2*z1,2);
    
    x1=ramos_cc(2*z-1,2);
    y1=ramos_cc(2*z,2);
    
    %procurar na matriz I_ramos
    for a=1:NL
        for b=3:size2
            
         if x1==I_ramos(2*a-1,b) && y1==I_ramos(2*a,b) && ( M3(x1,y1)==2 || M3(x1,y1)==3)
            ramos_cc(2*z-1,3)=I_ramos(2*a-1,1);
            ramos_cc(2*z,3)=I_ramos(2*a,1);
            
            x1=ramos_cc(2*z-1,3);
            y1=ramos_cc(2*z,3);
            
            cont6=4;
            
            for a1=(a+1):NL
                for b1=3:size2
                             if x1==I_ramos(2*a1-1,b1) && y1==I_ramos(2*a1,b1) && (M3(x1,y1)==2 || M3(x1,y1)==3)
                                 ramos_cc(2*z-1,cont6)=I_ramos(2*a1-1,1);
                                 ramos_cc(2*z,cont6)=I_ramos(2*a1,1);
            
                                 x1=ramos_cc(2*z-1,cont6);
                                 y1=ramos_cc(2*z,cont6);
            
                                 cont6=cont6+1;
                             end    
                end
            end
            
         end
         
        end
    end
end
end
end

% Somatório das impedâncias
MZ_CC2=zeros(3*NN,3);

size_CC=size(ramos_cc);
size_CC1=size_CC(1,1);
size_CC2=size_CC(1,2);

for a= 2:size_CC1/2
    for b=2:size_CC2
        
        for c=1:NL

        if ramos_cc(2*a-1,b)==posicoes_r(1,c) && ramos_cc(2*a,b)==posicoes_r(2,c);
            MZ_CC2(3*a-2:3*a,1:3)=MZ_CC2(3*a-2:3*a,1:3)+ MZ(1:3,3*c-2:3*c);
        end
        
        %atribuir metade de Z_1_2 para o nó Z_1
        if posicoes_r(1,c)==2 && posicoes_r(2,c)==ref_c
            MZ_CC2(1:3,1:3)=MZ(1:3,3*c-2:3*c)/2;
        end
        
        end
        
    end
end

%%

% Transformar a matriz MZ em "zeros" e "uns"
for x=1:3
    for y=1:3*NL
        if MZ(x,y)>0.010
            MZ1(x,y)=1;
        else
            MZ1(x,y)=0;
        end
    end
end
%Utilizar apenas a diagonal principal
    for y=1:NL
        MZ2(1,y)=MZ1(1,3*y-2);
        MZ2(2,y)=MZ1(2,3*y-1);
        MZ2(3,y)=MZ1(3,3*y-0);
    end


%Identificar as fases com tensão

            V_nos3(1,1)=V_nos1(1,1);
            V_nos3(2,1)=V_nos1(2,1);
            V_nos3(3,1)=V_nos1(3,1);

for x=1:(NN-1)
   
    a1=V_nos(2*x-1,2);
    a2=V_nos(2*x,2);
    
    for y=1:NL
        if a1==posicoes_r(1,y) && a2==posicoes_r(2,y)
            V_nos3(3*(x+1)-2,1)=V_nos1(3*(x+1)-2,1).*MZ2(1,y);
            V_nos3(3*(x+1)-1,1)=V_nos1(3*(x+1)-1,1).*MZ2(2,y);
            V_nos3(3*(x+1)-0,1)=V_nos1(3*(x+1)-0,1).*MZ2(3,y);
            
        end
    end
    
end

VrefT = Vref*1.6495;

res1=abs(V_nos3)/VrefT*sqrt(3);
res2=abs(I_ramos1);

V_S = (abs(V_nos3)/Vref);

c2 = 0;

for c1=1:NL
    ISS(c1,1) = res2(2*c1-1+c2,1);
    ISS(c1,2) = res2((2*c1)+c2,1);
    ISS(c1,3) = res2((2*c1+1)+c2,1);
    c2 = c2 + 1;
end

c2 = 0;

for c1=1:NN
    VSS(c1,1) = V_S(2*c1-1+c2,1);
    VSS(c1,2) = V_S((2*c1)+c2,1);
    VSS(c1,3) = V_S((2*c1+1)+c2,1);
    c2 = c2 + 1;
end

VSS(1,:)=[1 1 1];

%Coloca a tensão na sequência dos dados de entrada
for x = 1:2*(NN-1)
    V_nosX(x,1) = V_nos(x,1);
end
for x = 1:(NN-1)
    V_nosXX(1:2,x) = [V_nosX(2*x-1);V_nosX(2*x)];
end
VSSX = VSS;
for x = 2:NN
    for y = 1:NN-1
        if posicoes(1:2,x) == V_nosXX(1:2,y)
            VSSX(x,1:3) = VSS(y+1,1:3);
        end
    end
end

%Encontra o Modo de Operação da GD
MOD(NN,1:3) = 0;
SN_S_T = SN_S_T'; 
SN_FP = SN_FP';

for x = 1:NN
    for y = 1:3
        if SN_S_T(x,y) ~= 0
            if SN_FP(x,y) > 1
                SN_FP(x,y) = -1.2222*VSSX(x,y)+2.2222;
            end
            if SN_FP(x,y) == 0.9 && abs(VSSX(x,y)) >= 0 && abs(VSSX(x,y)) < 0.1
                MOD(x,y) = 1; %Modo 1
            end
            if SN_FP(x,y) == 0.9 && VSSX(x,y) >= 0.1 && VSSX(x,y) < 0.5
                MOD(x,y) = 2; %Modo 2
            end
            if SN_FP(x,y) == (0.333*VSSX(x,y)+0.733)
                MOD(x,y) = 3; %Modo 3
            end
            if SN_FP(x,y) == 1 && abs(VSSX(x,y)) >= 0.80 && abs(VSSX(x,y)) < 1
                MOD(x,y) = 4; %Modo 4
            end
            if SN_FP(x,y) == -1.2222*VSSX(x,y)+2.2222
                MOD(x,y) = 5; %Modo 5
            end
        end
        if SN_S_T(x,y) == 0 || (PPA(y,x)+PPR(y,x)) == 0
            MOD(x,y) = -1;
        end
        if MOD(x,y) == 0
            MPL = x;
            MPC = y;
        end
    end
end

contM = 0;

for x = 1:NN
    for y = 1:3
        if MOD(x,y) == -1 || MOD(x,y) == 1 || MOD(x,y) == 2 || MOD(x,y) == 3 || MOD(x,y) == 4 || MOD(x,y) == 5
            contM = contM + 1;
        end
    end
end

if contM == NN*3
    Modifica = 0; 
end

if Modifica == 1
    if VSSX(MPL,MPC) >= 0 && VSSX(MPL,MPC) < 0.1
        SN_FP(MPL,MPC) = 0.9000;
    end
    if VSSX(MPL,MPC) > 0.1 && VSSX(MPL,MPC) < 0.5
        SN_FP(MPL,MPC) = 0.9000;
    end
    if VSSX(MPL,MPC) >= 0.5 && VSSX(MPL,MPC) < 0.8
        SN_FP(MPL,MPC) = 0.333*VSSX(MPL,MPC)+0.733;
    end
    if VSSX(MPL,MPC) >= 0.8 && VSSX(MPL,MPC) < 1
        SN_FP(MPL,MPC) = 1.0000;
    end
    if VSSX(MPL,MPC) >= 1 && VSSX(MPL,MPC) < 1.1
       SN_FP(MPL,MPC) = -1.2222*VSSX(MPL,MPC)+2.2222;
    end
    
    PPR(MPC,MPL) = abs((PPA(MPC,MPL).*SN_P_T(MPC,MPL).*(1-SN_FP(MPL,MPC))./(SN_FP(MPL,MPC).*SN_Q_T(MPC,MPL)*j)));
       
%     %Modificação da curva B
%     if PPR(MPC,MPL) > 0.44
%         PPR(MPC,MPL) = 0.44;
%     end
    
    SN(MPL,MPC) = (1-PPA(MPC,MPL)).*SN_P_Y(MPC,MPL)'+(1-PPR(MPC,MPL)).*SN_Q_Y(MPC,MPL)'*j;
    SN(MPL,MPC) = (1-PPA(MPC,MPL)).*SN_P_D(MPC,MPL)'+(1-PPR(MPC,MPL)).*SN_Q_D(MPC,MPL)'*j;
end

SN_S_T = SN_S_T';
SN_FP = SN_FP';

% prompt = 'PARA';
% NC = input(prompt);

end

MOD
VSSX
SN_FP'
PPR'

toc
%ISS

VSS = VSSX;
VSS(1,:)=[1 1*(cosd(-120)+j*sind(-120)) 1*(cosd(120)+j*sind(120))];

MY_CC = 0;

%%%%Obtém a matriz de Admitância
for i1 = 0:NN-1
    i3 = 0;
    for i2 = 1+3*i1:3+3*i1
        i3 = i3 + 1;
        MY_C(i3,1:3) = MZ_CC2(i2,1:3);
    end
    MY_C(1:3,1:3) = inv(MY_C(1:3,1:3));
    YF(1+3*i1:3+3*i1,1:3) = MY_C(1:3,1:3);
    
    YFS(i1+1,1:3) = [YF(1+3*i1,1)+YF(1+3*i1,2)+YF(1+3*i1,3); YF(2+3*i1,1)+YF(2+3*i1,2)+YF(2+3*i1,3); YF(3+3*i1,1)+YF(3+3*i1,2)+YF(3+3*i1,3)];
end

%%%%Falta Monofásica
for i1 = 0:NN-1

%Entrada IEEE
    if VSS(i1+1,1) ~= 0
        VSS(i1+1,1) = [1.3728+0.0000i];
    end
    if VSS(i1+1,2) ~= 0
        VSS(i1+1,2) = [-0.6864-1.1890i];
    end
    if VSS(i1+1,3) ~= 0
        VSS(i1+1,3) = [-0.6864+1.1890i];
    end
    
%Fase A
    Ca = [1 0 0 YF(1+3*i1,1) YF(1+3*i1,2) YF(1+3*i1,3) YFS(i1+1,1);
          0 1 0 YF(2+3*i1,1) YF(2+3*i1,2) YF(2+3*i1,3) YFS(i1+1,2);
          0 0 1 YF(3+3*i1,1) YF(3+3*i1,2) YF(3+3*i1,3) YFS(i1+1,3);
          0 0 0 1         0         0         0       ;
          0 0 1 0         0         0         0       ;
          0 1 0 0         0         0         0       ;
          0 0 0 0         0         0         1       ]; 
    IPa= YF(1+3*i1:3+3*i1,1:3)*(1000)*[VSS(i1+1,1);VSS(i1+1,2);VSS(i1+1,3)];  IPa = [IPa(1);IPa(2);IPa(3);0;0;0;0];
    Xa(1:7,i1+1) = abs(inv(Ca)*IPa);
%Fase B
    Cb = [1 0 0 YF(1+3*i1,1) YF(1+3*i1,2) YF(1+3*i1,3) YFS(i1+1,1);
          0 1 0 YF(2+3*i1,1) YF(2+3*i1,2) YF(2+3*i1,3) YFS(i1+1,2);
          0 0 1 YF(3+3*i1,1) YF(3+3*i1,2) YF(3+3*i1,3) YFS(i1+1,3);
          0 0 0 0         1         0         0       ;
          1 0 0 0         0         0         0       ;
          0 0 1 0         0         0         0       ;
          0 0 0 0         0         0         1       ];
    IPb= YF(1+3*i1:3+3*i1,1:3)*(1000)*[VSS(i1+1,1);VSS(i1+1,2);VSS(i1+1,3)];  IPb = [IPb(1);IPb(2);IPb(3);0;0;0;0];
    Xb(1:7,i1+1) = abs(inv(Cb)*IPb);
%Fase C
    Cc = [1 0 0 YF(1+3*i1,1) YF(1+3*i1,2) YF(1+3*i1,3) YFS(i1+1,1);
          0 1 0 YF(2+3*i1,1) YF(2+3*i1,2) YF(2+3*i1,3) YFS(i1+1,2);
          0 0 1 YF(3+3*i1,1) YF(3+3*i1,2) YF(3+3*i1,3) YFS(i1+1,3);
          0 0 0 0         0         1         0       ;
          1 0 0 0         0         0         0       ;
          0 1 0 0         0         0         0       ;
          0 0 0 0         0         0         1       ]; 
    IPc= YF(1+3*i1:3+3*i1,1:3)*(1000)*[VSS(i1+1,1);VSS(i1+1,2);VSS(i1+1,3)];  IPc = [IPc(1);IPc(2);IPc(3);0;0;0;0];
    Xc(1:7,i1+1) = abs(inv(Cc)*IPc);
end

%%%%Falta Bifásica
for i1 = 0:NN-1

%Entrada IEEE
    if VSS(i1+1,1) ~= 0
        VSS(i1+1,1) = [1.3728+0.0000i];
    end
    if VSS(i1+1,2) ~= 0
        VSS(i1+1,2) = [-0.6864-1.1890i];
    end
    if VSS(i1+1,3) ~= 0
        VSS(i1+1,3) = [-0.6864+1.1890i];
    end
    
%Fases AB
    Cab = [1 0 0 YF(1+3*i1,1) YF(1+3*i1,2) YF(1+3*i1,3) YFS(i1+1,1);
          0 1 0 YF(2+3*i1,1) YF(2+3*i1,2) YF(2+3*i1,3) YFS(i1+1,2);
          0 0 1 YF(3+3*i1,1) YF(3+3*i1,2) YF(3+3*i1,3) YFS(i1+1,3);
          0 0 0 1         0         0         0       ;
          0 0 0 0         1         0         0       ;
          0 0 1 0         0         0         0       ;
          1 1 0 0         0         0         0       ]; 
    IPab= YF(1+3*i1:3+3*i1,1:3)*(1000)*[VSS(i1+1,1);VSS(i1+1,2);VSS(i1+1,3)];  IPab = [IPab(1);IPab(2);IPab(3);0;0;0;0];
    Xab(1:7,i1+1) = abs(inv(Cab)*IPab);
%Fases AC
    Cac = [1 0 0 YF(1+3*i1,1) YF(1+3*i1,2) YF(1+3*i1,3) YFS(i1+1,1);
          0 1 0 YF(2+3*i1,1) YF(2+3*i1,2) YF(2+3*i1,3) YFS(i1+1,2);
          0 0 1 YF(3+3*i1,1) YF(3+3*i1,2) YF(3+3*i1,3) YFS(i1+1,3);
          0 0 0 1         0         0         0       ;
          0 0 0 0         0         1         0       ;
          0 1 0 0         0         0         0       ;
          1 0 1 0         0         0         0       ];
    IPac= YF(1+3*i1:3+3*i1,1:3)*(1000)*[VSS(i1+1,1);VSS(i1+1,2);VSS(i1+1,3)];  IPac = [IPac(1);IPac(2);IPac(3);0;0;0;0];
    Xac(1:7,i1+1) = abs(inv(Cac)*IPac);
%Fases CB
    Ccb = [1 0 0 YF(1+3*i1,1) YF(1+3*i1,2) YF(1+3*i1,3) YFS(i1+1,1);
          0 1 0 YF(2+3*i1,1) YF(2+3*i1,2) YF(2+3*i1,3) YFS(i1+1,2);
          0 0 1 YF(3+3*i1,1) YF(3+3*i1,2) YF(3+3*i1,3) YFS(i1+1,3);
          0 0 0 0         0         1         0       ;
          0 0 0 0         1         0         0       ;
          1 0 0 0         0         0         0       ;
          0 1 1 0         0         0         0       ]; 
    IPcb= YF(1+3*i1:3+3*i1,1:3)*(1000)*[VSS(i1+1,1);VSS(i1+1,2);VSS(i1+1,3)];  IPcb = [IPcb(1);IPcb(2);IPcb(3);0;0;0;0];
    Xcb(1:7,i1+1) = abs(inv(Ccb)*IPcb);
end

%%%%Falta Trifásica
for i1 = 0:NN-1
    if VSS(i1+1,1) == 0 || VSS(i1+1,2) == 0 || VSS(i1+1,3) == 0
       VSS(i1+1,1) = 0; VSS(i1+1,2) = 0; VSS(i1+1,3) = 0;
    end
    
%Entrada IEEE
    if VSS(i1+1,1) ~= 0 || VSS(i1+1,2) ~= 0 || VSS(i1+1,3) ~= 0
        VSS(i1+1,1) = [1.3728+0.0000i];
        VSS(i1+1,2) = [-0.6864-1.1890i];
        VSS(i1+1,3) = [-0.6864+1.1890i];
    end
    
    Ct = [1 0 0 YF(1+3*i1,1) YF(1+3*i1,2) YF(1+3*i1,3) YFS(i1+1,1);
          0 1 0 YF(2+3*i1,1) YF(2+3*i1,2) YF(2+3*i1,3) YFS(i1+1,2);
          0 0 1 YF(3+3*i1,1) YF(3+3*i1,2) YF(3+3*i1,3) YFS(i1+1,3);
          0 0 0 1         0         0         0       ;
          0 0 0 0         1         0         0       ;
          0 0 0 0         0         1         0       ;
          1 1 1 0         0         0         0       ];
    IP = YF(1+3*i1:3+3*i1,1:3)*(1000)*[VSS(i1+1,1);VSS(i1+1,2);VSS(i1+1,3)];  
    IPt= [IP(1);IP(2);IP(3);0;0;0;0];
    Xt(1:7,i1+1) = abs(inv(Ct)*IPt);  
end

for i1 = 1:NN
    IFM(i1,1:3) = [Xa(1,i1) Xb(2,i1) Xc(3,i1)]; %Falta Monofásica (A, B e C)
    IFB(i1,1:3) = [Xab(1,i1) Xac(3,i1) Xcb(2,i1)]; %Falta Bifásica (AB, AC e CB)
    IFT(i1,1:3) = [Xt(1,i1) Xt(2,i1) Xt(3,i1)]; %Falta Trifásica (ABC)
end

M=M3;

