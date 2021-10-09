function [yLMS, eLMS, wLMS, eaLMS]=LMSWn(M, muLMS, u, d)
NLMS=length(u); 
wLMS=zeros(M,1); 
yLMS=zeros(NLMS, 1);
eLMS=zeros(NLMS, 1); 
uvLMS=zeros(M,1);
eaLMS=zeros(NLMS, 1);

for i=1:NLMS
   uvLMS=[u(i); uvLMS(1:M-1)];
   yLMS(i)=uvLMS'*wLMS(:,i);
   eLMS(i)=d(i)-yLMS(i);
   wLMS(:,i+1)=wLMS(:,i)+muLMS*eLMS(i)*uvLMS; 
end
end