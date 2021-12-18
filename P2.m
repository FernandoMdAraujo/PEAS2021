%P2 - Parte Computacional
%Fernando Marques de Araujo - 6801112
clc
clear
close all
% Questao 1
% Dados da Questao
rs0 = 20;
rs1 = 19.39;
rs2 = 17.59;

sigma2_eta = 1;

% Modelo
% x(n) = s(n) +eta(n)

Rs = [rs0 rs1 rs2; rs1 rs0 rs1; rs2 rs1 rs0];
Reta = sigma2_eta * eye(3);

Rx = Rs + Reta

[Qa,Da] = eig(Rx)

%item d
load dadosAP2; % O arquivo dadosAP2.mat contém as variáveis x1, x2 e x3.
% i) Estimativa da matriz de covariâcia Cx e PCA com base na decomposição
% de autovalores e autovetores de Cx.
%x1=(x1-mean(x1));
x2=(x2-mean(x2));
x3=(x3-mean(x3));

X = [x1; x2; x3];
T = size(X,2);
Cx = 1/T * (X*X')
[Qi, Di] = eig(Cx)
y_eig = Qi'*X;
##display(var(y_eig,0,3), 'var(y_eig)');
Proporcao=diag(Di)*100/sum(diag(Di));
[diag(Di) Proporcao]

% iii) PCA com base na decomposição
%       de autovalores e autovetores de Cx.

K=8758;
figure(1)
subplot(3,1,1)
plot(X(1,1:K), '.'); grid on
ylabel('x_1')
title('Sinais no tempo Antes do PCA')
subplot(3,1,2)
plot(X(2,1:K), '.'); grid on
ylabel('x_2')
subplot(3,1,3)
plot(X(3,1:K), '.'); grid on
ylabel('x_3')
figure(2)
subplot(3,1,1)
plot(y_eig(1,1:K), '.'); grid on
ylabel('PCA: y_1')
title('Sinais no tempo Depois do PCA')
subplot(3,1,2)
plot(y_eig(2,1:K), '.'); grid on
ylabel('PCA: y_2')
subplot(3,1,3)
plot(y_eig(3,1:K), '.'); grid on
ylabel('PCA: y_3')
xlabel('n')

% iv
figure(3)
subplot(2,2,1)
plot(X(1,:),X(2,:), '.')
hold on
plot(2*Qi(1,3)*[-10 10], 2*Qi(2,3)*[-10 10],'g','Linewidth',2)
plot(2*Qi(1,2)*[-10 10], 2*Qi(2,2)*[-10 10],'m','Linewidth',2)
plot(2*Qi(1,1)*[-10 10], 2*Qi(2,1)*[-10 10],'r','Linewidth',2)
hold off
title('Diagrama de dispersão')
xlabel('x_1(n)'); ylabel('x_2(n)')
grid
axis equal
subplot(2,2,2)
plot(X(1,:),X(3,:), '.')
hold on
plot(2*Qi(1,3)*[-10 10], 2*Qi(2,3)*[-10 10],'g','Linewidth',2)
plot(2*Qi(1,2)*[-10 10], 2*Qi(2,2)*[-10 10],'m','Linewidth',2)
plot(2*Qi(1,1)*[-10 10], 2*Qi(2,1)*[-10 10],'r','Linewidth',2)
title('Diagrama de dispersão')
xlabel('x_1(n)'); ylabel('x_3(n)')
grid
axis equal
subplot(2,2,3)
plot(X(2,:),X(3,:), '.')
hold on
plot(2*Qi(1,3)*[-10 10], 2*Qi(2,3)*[-10 10],'g','Linewidth',2)
plot(2*Qi(1,2)*[-10 10], 2*Qi(2,2)*[-10 10],'m','Linewidth',2)
plot(2*Qi(1,1)*[-10 10], 2*Qi(2,1)*[-10 10],'r','Linewidth',2)
title('Diagrama de dispersão')
xlabel('x_2(n)'); ylabel('x_3(n)')
grid
axis equal



% v) Histogramas
B=50;
figure(4)
subplot(3,2,1)
hist(X(1,:),B); grid
title('Histograma de x_1(n)')
%axis([-10 10 0 500])

subplot(3,2,3)
hist(X(2,:),B); grid
title('Histograma de x_2(n)')
%axis([-5 5 0 500])

subplot(3,2,5)
hist(X(3,:),B); grid
title('Histograma de x_3(n)')
%axis([-5 5 0 500])
% 
subplot(3,2,2)
hist(y_eig(1,:),B);grid
title('Histograma de y_1(n)')
%axis([-5 5 0 500])

subplot(3,2,4)
hist(y_eig(2,:),B);grid
title('Histograma de y_2(n)')
%axis([-10 10 0 500])

subplot(3,2,6)
hist(y_eig(3,:),B);grid
title('Histograma de y_3(n)')
%axis([-10 10 0 500])
% 

##% extra
##figure(5)
##subplot(2,2,1)
##plot(y_eig(2,:),y_eig(1,:), '.')
##hold on
##title('Diagrama de dispersão')
##xlabel('y_2(n)'); ylabel('y_1(n)')
##grid
##%axis equal
##subplot(2,2,2)
##plot(y_eig(3,:),y_eig(1,:), '.')
##hold on
##title('Diagrama de dispersão')
##xlabel('y_3(n)'); ylabel('y_1(n)')
##grid
##%axis equal
##subplot(2,2,3)
##plot(y_eig(3,:),y_eig(2,:), '.')
##hold on
##title('Diagrama de dispersão')
##xlabel('y_3(n)'); ylabel('y_2(n)')
##grid
##%axis equal
##y1=y_eig(1,:);
##y2=y_eig(2,:);
##y3=y_eig(3,:);
##figure(6)
##subplot(322);
##scatter(y1, y2, 1);
##grid;
##xlabel('y1(n)');
##ylabel('y2(n)');
##title('Diagrama de dispersão de y1(n) por y2(n)');
##
##subplot(324);
##scatter(y1, y3, 1);
##grid;
##xlabel('y1(n)');
##ylabel('y3(n)');
##title('Diagrama de dispersão de y1(n) por y3(n)');
##
##subplot(326);
##scatter(y2, y3, 1);
##grid;
##xlabel('y2(n)');
##ylabel('y3(n)');
##title('Diagrama de dispersão de y2(n) por y3(n)');

%Questao 2
%item a
N=K;
L=K; % comprimento da janela 100 %
K=1;%round((N-L)/R);
R=1;
Nf=N; % usar por exemplo Nf=N
Nf2=round(Nf/2);
% 
Ipm1=zeros(1,Nf); Ipm2=zeros(1,Nf);
janela1= rectwin(L)';
%janela2= blackman(L)'; 
%  
for r=0:K-1
 xr1(1:L)=x1(r*R+[1:L]).*janela1; 
 Ipm1=Ipm1+abs(fft(xr1(1:L),Nf)).^2;
end
% 
U1=sum(janela1.^2)/L; Um1=sum(janela1)^2;  
% Estimativa da DSP: Equacao slide 10
Ipm1=Ipm1/(L*K*U1); % 
% % 
figure(6)
plot([0:Nf2-1],10*log10(Ipm1(1:Nf2)),'r');
title('media dos Periodogramas janelados I_{r,Nf}(k) (dB)')
xlabel('k = 0:Nf/2;  (fft com Nf pontos)')
legend('janela retangular')
figure(7)
plot([0:Nf2-1],(Ipm1(1:Nf2)),'r');
grid
title('media dos Periodogramas janelados I_{r,Nf}(k) - Escala linear')
xlabel('k = 0:Nf/2;  (fft com Nf pontos)')
legend('janela retangular')

%item b
varR1E = mean(Ipm1(500:Nf/2))
[vy1,vx1] = max(abs(Ipm1))
[vy2,vx2] = max(abs(Ipm1(Ipm1<max(Ipm1))))
wsinal1 = 2*pi*(vx1-1)/Nf
wsinal2 = 2*pi*(vx2-1)/Nf
A1 = sqrt((vy1-varR1E)*4*U1*L/(Um1))
A2 = sqrt((vy2-varR1E)*4*U1*L/(Um1))
