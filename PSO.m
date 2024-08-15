
function PSO   
clear all;
clc;
tic
format long; 
a=0;  
b=1;   
vmax=0.25*(b-a);       
N=100;               
c1=2;%1.2;  
c2=2;%1.2;              
wmax=0.9;           
wmin=0.4;          
D=30;                    
Tmax=100;            
Py=zeros(N,2);    
A=zeros(N,2);      
Xa=[];

for i=1:N                   
    for j=1:D
        x(i,j)=rand(1);  
        v(i,j)=0; 
    end
end



for i=1:N                           
    Py(i,1)=f1(x(i,1:D));             
    Py(i,2)=f2(x(i,1:D));           
    Xy(i,1:D)=x(i,1:D);             
    A(i,1)=Py(i,1);                 
    A(i,2)=Py(i,2);
    Xa=[Xa;x(i,1:D)];               
end
    disp('Xa=')
    Xa
%----------A集合成员的支配关系比较，得到当前的非支配解集C----------%
[C,Xc]=ParetoC(A,Xa);      
    disp('Xc=')
    Xc
%----------从C中随机取一个值，将其位置值赋给Pg----------%
B=C;Xb=Xc;          %存放最终非劣解
%n0=randint(1,1,[1,size(B,1)]);     
%Pg=Xb(n0,2:D);        

Pg=leader(B,Xb);

clear('A');  
%----------主循环，迭代过程----------%
for t=1:Tmax                           
    Xa=[];C=[];Xc=[];
    if ~mod(t,10)
        fprintf('current iter:%d\n',t)
    end
    w=wmax-(wmax-wmin)*t/Tmax; 
    for i=1:N                             
        v(i,2:D)=w*v(i,2:D)+c1*rand*(Xy(i,2:D)-x(i,2:D))+c2*rand*(Pg-x(i,2:D)); 
        for j=2:D
                if v(i,j)>vmax
                   v(i,j)=vmax;
                elseif v(i,j)<-vmax
                       v(i,j)=-vmax;
                end
        end
        x(i,2:D)=x(i,2:D)+v(i,2:D); 
        for j=1:D
            if(x(i,j)>b)   x(i,j)=b;   end   
            if(x(i,j)<a)   x(i,j)=a;   end    
        end

        A(i,1)=f1(x(i,1:D));
        A(i,2)=f2(x(i,1:D));
        Xa=[Xa;x(i,1:D)]; 
        
        if ((A(i,1)<Py(i,1))&&(A(i,2)<=Py(i,2)))        
            Py(i,1)=A(i,1); Py(i,2)=A(i,2);
            Xy(i,1:D)=x(i,1:D);   
        elseif ((A(i,1)<=Py(i,1))&&(A(i,2)<Py(i,2)))    
               Py(i,1)=A(i,1); Py(i,2)=A(i,2);
               Xy(i,1:D)=x(i,1:D); 
        elseif ((Py(i,1)<A(i,1))&&(Py(i,2)<=A(i,2)))       
                Py(i,1)=Py(i,1);Py(i,2)=Py(i,2);Xy(i,1:D)=Xy(i,1:D);
        elseif ((Py(i,1)<=A(i,1))&&(Py(i,2)<A(i,2)))         
                 Py(i,1)=Py(i,1);Py(i,2)=Py(i,2);Xy(i,1:D)=Xy(i,1:D);
        else    
            warning('OFF');
            m=randint(1,1,[1,2]);
            if m==1
               Py(i,1)=A(i,1); Py(i,2)=A(i,2);
               Xy(i,1:D)=x(i,1:D);
            end
        end
    end

    
    [C,Xc]=ParetoC(A,Xa);
            
    
    for i=1:size(B,1)
        C=[C;B(i,:)];    
        Xc=[Xc;Xb(i,:)];
    end
    [B,Xb]=ParetoC(C,Xc);
    if size(B,1)>100  
    n1=randint(100,1,[1,size(B,1)]);
    B=B(n1(:),:);   
    Xb=Xb(n1(:),:);
    end
    disp('Xb=')
    Xb
   
    n=randint(1,1,[1,size(B,1)]);      
    Pg=Xb(n,2:D);                     
    Pg=leader(B,Xb);
    disp('Pg')
    size(Pg)
    clear('A');
    
    if mod(t,20)==0
           plot(B(:,1),B(:,2),'ro','linewidth',1.5) 
           xlabel('f1(x)');ylabel('f2(x)');
           grid on
        else
            plot(B(:,1),B(:,2),'ro','linewidth',1.5) 
            xlabel('f1(x)');ylabel('f2(x)');
            grid on
            title(strcat('u',num2str(t),'P'));
        end
        title(strcat('u',num2str(t),'Pareto'));
        pause(0.2)
end
toc
save solutionWW-ZDT1-MOPSOceshi.txt Xb -ASCII
save solution-ZDT1-MOPSOceshi.txt B -ASCII

plot(B(:,1),B(:,2),'ro','linewidth',1.5); 
xlabel('f1(x)');
ylabel('f2(x)');
title(strcat('u',num2str(t),'A'));
grid on




   
  