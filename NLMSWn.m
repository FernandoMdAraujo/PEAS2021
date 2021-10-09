%Função - NLMS
function [yNLMS, eNLMS, wNLMS, eaNLMS]= NLMSWn(M, muNLMS, u, d, deltaNLMS)
NNLMS=length(u); 
wNLMS=zeros(M,1); 
yNLMS=zeros(NNLMS, 1);
eNLMS=zeros(NNLMS, 1); 
uvNLMS=zeros(M,1);
eaNLMS=zeros(NNLMS, 1);

for i=1:NNLMS
uvNLMS=[u(i); uvNLMS(1:M-1)];
yNLMS(i)=uvNLMS'*wNLMS(:,i);
eNLMS(i)=d(i)-yNLMS(i);
muNLMSusado=muNLMS/(uvNLMS'*uvNLMS+deltaNLMS);
wNLMS(:,i+1)=wNLMS(:,i)+muNLMSusado*eNLMS(i)*uvNLMS;

end
end