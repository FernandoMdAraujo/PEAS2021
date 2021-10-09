%Fernando Marques de Araujo - NUP 6801112

%Creditos pelo Codigo:
%Este script foi feito com base na Aula #11 da Professora Maria D. Miranda
%Aula ministrada em 28/set/2021
%Curso PTC3451 – Processamento Estatístico e Adaptativo de Sinais

%Parte Computacional P1

clear
close all

%Sinal de voz
[ou, fs]=audioread('entradaP1.wav');
[do, fs]=audioread('desejadoP1.wav');
%sound(ou, fs); pause; sound(do, fs);

N=length(ou);

Ta=1/fs;

M=800;

deltaNLMS=1e-12;
muNLMS=2;
muLMS=0.01;
[yNLMS, eNLMS, wNLMS, eaNLMS]= NLMSWn(M, muNLMS, ou, do, deltaNLMS);
[yLMS, eLMS, wLMS, eaLMS]=LMSWn(M, muLMS, ou, do);

v1=[0.06 8.7 -1.2 1.2];
v2=[0.06 8.7 -0.6 0.6];
v3=[0.06 8.7 -150 0]

%Plot do sinal Eco
%figure(1)
title('erro - sinal sem eco');
subplot(411); plot ([0:N-1]/fs,ou(1:N), 'y'); title('Entrada do FA')
axis(v1);grid
subplot(412); plot ([0:N-1]/fs,do(1:N), 'r'); title('Desejado')
axis(v2);grid
subplot(413); plot ([0:N-1]/fs,10*log10(eNLMS(1:N).^2)); title('Erro^2 (dB) NLMS')
axis(v3);grid
subplot(414); plot ([0:N-1]/fs,10*log10(eLMS(1:N).^2)); title('Erro^2 (dB) LMS')
axis(v3);grid

%Foi deixado em aquivo separado, por isso aparece comentado aqui
% function [yNLMS, eNLMS, wNLMS, eaNLMS]= NLMSWn(M, muNLMS, u, d, deltaNLMS)
% NNLMS=length(u); 
% wNLMS=zeros(M,1); 
% yNLMS=zeros(NNLMS, 1);
% eNLMS=zeros(NNLMS, 1); 
% uvNLMS=zeros(M,1);
% eaNLMS=zeros(NNLMS, 1);
% 
% for i=1:NNLMS
% uvNLMS=[u(i); uvNLMS(1:M-1)];
% yNLMS(i)=uvNLMS'*wNLMS(:,i);
% eNLMS(i)=d(i)-yNLMS(i);
% muNLMSusado=muNLMS/(uvNLMS'*uvNLMS+deltaNLMS);
% wNLMS(:,i+1)=wNLMS(:,i)+muNLMSusado*eNLMS(i)*uvNLMS;
% 
% end
% end

% function [yLMS, eLMS, wLMS, eaLMS]=LMSWn(M, muLMS, u, d)
% NLMS=length(u); 
% wLMS=zeros(M,1); 
% yLMS=zeros(NLMS, 1);
% eLMS=zeros(NLMS, 1); 
% uvLMS=zeros(M,1);
% eaLMS=zeros(NLMS, 1);
% 
% for i=1:NLMS
%    uvLMS=[u(i); uvLMS(1:M-1)];
%    yLMS(i)=uvLMS'*wLMS(:,i);
%    eLMS(i)=d(i)-yLMS(i);
%    wLMS(:,i+1)=wLMS(:,i)+muLMS*eLMS(i)*uvLMS; 
% end
% end