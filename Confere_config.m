
T = size(M);
T_l=T(1,1);
T_c=T(1,2);
a = 0; 
f = 0;
q = 0;
S = 0;
T_Traf = 0;

%Simplificar
while q == 0 
    for r = 1:T_l
        if M(r,T_c) == 1 || M(r,T_c) == 2 || M(r,T_c) == 3 || M(r,T_c) == 4
            q = q+1;
         else q = q+0;
        end
        if M(r,T_c) == 3
            T_Traf = 1;
        end
    end
    if q == 0
        T_c = T_c-1;
        M = M(:,1:T_c);
        end
end 
T = size(M);
T_l=T(1,1);
T_c=T(1,2);
q = 0;
while q == 0 
    for r = 1:T_c
        if M(T_l,r) == 1 || M(T_l,r) == 2 || M(T_l,r) == 3 || M(T_l,r) == 4
            q = q+1;
         else q = q+0;
        end
    end
    if q == 0
        T_l = T_l-1;
        M = M(1:T_l,:);
        end
end 

%Esquerda
q = 0;
T_c=1;
%Simplificar
while q == 0 
    for r = 1:T_l
        if M(r,T_c) == 1 || M(r,T_c) == 2 || M(r,T_c) == 3 || M(r,T_c) == 4
            q = q+1;
         else q = q+0;
        end
    end
    if q == 0
        M(:,1) = [];
        end
end 

T = size(M);
T_l=T(1,1);
T_c=T(1,2);

%Confere linha
for l = 1:T_l
    x=0;
    for c = 1:T_c
        if M(l,c)==1 || M(l,c)==2 || M(l,c)==3 || M(l,c)==4
            x=x+1;
        end
    end
    if x == 0
       S = 1;
    end
end

%Confere coluna
x = 0;
if S ~= 1
    for c = 1:T_c
        x=0;
        for l = 1:T_l
            if M(l,c)==1 || M(l,c)==2 || M(l,c)==3 || M(l,c)==4
                x=x+1;
            end
        end
        if x == 0
           S = 1;
        end
    end
end

f = 0;
NN = 0;
NL = 0;

%Confere vizinhança
if S ~= 1
    for l = 1:T_l
        for c = 1:T_c
            if M(l,c) == 1
                NN = NN + 1;
            elseif M(l,c) == 2 || M(l,c) == 3 || M(l,c) == 4
                NL = NL + 1;    
            end
            
            if M(l,c) == 2 || M(l,c) == 3 || M(l,c) == 4
                if l ~= 1 
                    if M(l-1,c) == 1
                        a = a + 1;
                    end
                end
                if l ~= T_l
                    if M(l+1,c) == 1
                        a = a + 1;
                    end
                end
                if c ~= 1 
                    if M(l,c-1) == 1
                        a = a + 1;
                    end
                end
                if c ~= T_c
                    if M(l,c+1) == 1
                        a = a + 1;
                    end
                end
                if a < 2
                   f = f + 1;
                end
                a=0;
            end
        end
    end
end

if f > 0
    S = 1;
end
N=0;
if S ~= 1
    for l = 1:T_l
        for c = 1:T_c
            N = M(l,c) + N;
        end
    end
    else error('Matriz incorreta')
end

%Encontra a matriz de linhas
r = 1;
aux3=1;
if S ~= 1
    for l = 1:T_l
        for c = 1:T_c
            a=0;
            if M(l,c) == 2 || M(l,c) == 3 || M(l,c) == 4
                amr(aux3,:)=[l,c];
                aux3=aux3+1;
                
                if l ~= 1 
                    if M(l-1,c) == 1
                        a = a + 1;
                        amr(aux3,:)=[l-1,c];
                        aux3=aux3+1;     
                    end
                end
                if l ~= T_l
                    if M(l+1,c) == 1
                        a = a + 1;
                        amr(aux3,:)=[l+1,c];
                        aux3=aux3+1;
                    end
                end
                if c ~= 1 
                    if M(l,c-1) == 1
                        a = a + 1;
                        amr(aux3,:)=[l,c-1];
                        aux3=aux3+1;
                    end
                end
                if c ~= T_c
                    if M(l,c+1) == 1
                        a = a + 1;
                        amr(aux3,:)=[l,c+1];
                        aux3=aux3+1;
                    end
                end
                if a == 0
                    a
                   f = f + 1;
                end
            end
        end
    end
end


