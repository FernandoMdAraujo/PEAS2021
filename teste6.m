clear all;
close all;
clc

N=500;
h=[0.25 1 0.25];
VarRuido=0.01;
simb1=2*rand(1,N)-1;
a=1*(simb1>=0)-1*(simb1<0); % sinal binário
yo=filter(h,1,a);
eta=sqrt(VarRuido)*randn(1,N);
ya=yo+eta;

A=[0 0 0; 1 0 0; 0 1 0];
b=[1 0 0]';
H=h;
n=1:N;
p10=10;

x=zeros(3,N+1);
y=zeros(1,N);
ySR=zeros(1,N);

v_x = [1 0 0];
v_y = 0.02;

x=zeros(3,N+1);
vx = (2*randi([0 1], 1, N) - 1);
vy = sqrt(v_y)*randn(1,N);

for i=1:N
    x(:,i+1)=A*x(:,i)+vx(i)*v_x';
    y(:,i)=H*x(:,i)+vy(:,i); %Observacao com ruido
    ySR(:,i)=H*x(:,i); %Obsevacao SEM ruido
end


%Algoritimo de Kalman
%Inicializacao
xpri=zeros(3,N);
xpos=zeros(3,N);
Ppos=zeros(3,3);

Ppri=p10*eye(3);
Sy=v_y;

%loop para varrer todos os valores de n
##for i=1:N-1 %Segunda Versao do Kalman
##  S=H*Ppri*H'+Sy;
##  Kappa=A*Ppri*H'*inv(S);
##  alpha=y(:,i)-H*xpri(:,i);
##  xpri(:,i+1)=A*xpri(:,i)+Kappa*alpha;
##  Ppri=A*Ppri*A'-Kappa*H*Ppri*A'+vx*vx';
##end

for i=1:N-1 %Primeira Versao do Kalman
  S=H*Ppri*H'+Sy;
  K=A*Ppri*H'*inv(S);
  alpha=y(:,i)-H*xpri(:,i);
  xpos(:,i)=xpri(:,i)+K*alpha;
  Ppos=(eye(3)-K*H)*Ppri;
  xpri(:,i+1)=A*xpri(:,i)+K*alpha;
  Ppri=A*Ppri*A'+vx*vx';
end


figure(1)
subplot(311)
plot(n,a,'.b');
title('item a - Sinal Aleatório a(n)');
axis([0 N -1.5 1.5])
subplot(312)
plot(n,a,'.b');
title('item a - Segmento de 100 amostras do Sinal Aleatório a(n)');
axis([0 100 -1.5 1.5])
subplot(313)
plot(n,ya);
title('iten b - Sinal Observado y(n)');

figure(2)
subplot(311);
plot(y(1,:),'r','Linewidth',2)
hold on
plot(ySR(1,:),'s')
hold off; grid
legend('y(n) com Ruido','y(n)Sem Ruido')
title('Comparação entre os sinais y(n) com ruído e sem ruído');

subplot(312);
plot(y(1,:),'r','Linewidth',2)
hold on
plot(ySR(1,:),'s')
plot(xpri(1,:),xpri(2,:),xpri(3,:),'g' ,'Linewidth',2)
hold off; grid
legend('y(n) com Ruido','y(n)Sem Ruido',' ','x a priori')
title('Resultado obtido com Filtro de Kalman a priori');

subplot(313);
plot(y(1,:),'r','Linewidth',2)
hold on
plot(ySR(1,:),'s')
plot(xpos(1,:),xpos(2,:),xpos(3,:),'g','Linewidth',3)
hold off; grid
legend('y(n) com Ruido','y(n)Sem Ruido',' ','x a posteriori')
title('Resultado obtido com Filtro de Kalman a posteriori');

