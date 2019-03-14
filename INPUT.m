function [AFA,HBETA,SHETA0,a,b,JR,COOR,MEL,TOM,WDE,N]=INPUT(ft,NP,NE,NM)
%%
%读取结点坐标；
COOR=zeros(2,NP);
for I=1:NP
    C=fscanf(ft,'%d %f %f',3);
    IP=C(1);COOR(1,IP)=C(2);COOR(2,IP)=C(3);
end
%%
%读取每个单元包含的结点号码；
MEL=zeros(5,NE);
for I=1:NE
    C=fscanf(ft,'%d %d %d %d %d %d',6);
    NEE=C(1);NME=C(2);MEL(1,NEE)=C(3);MEL(2,NEE)=C(4);MEL(3,NEE)=C(5);MEL(4,NEE)=C(6);
    MEL(5,NEE)=NME;
end
%%
%结点总的自由度；
JR=zeros(NP,1);
for I=1:NP
    JR(I)=1;
end
N=0;
for I=1:NP
    if (JR(I)==1)
        N=N+1;
        JR(I)=N;
    end
end
%%
%获取材料参数和绝热温升的参数；
TOM=zeros(4,NM);
WDE=zeros(3,NM);
for J=1:NM
    C=fscanf(ft,'%d %f %f %f %f %f %f %f',8);
    JJ=C(1);TOM(1,JJ)=C(2);TOM(2,JJ)=C(3);TOM(3,JJ)=C(4);TOM(4,JJ)=C(5);WDE(1,JJ)=C(6);WDE(2,JJ)=C(7);WDE(3,JJ)=C(8);
    AFA=24*TOM(4,JJ)/(TOM(1,JJ)*TOM(2,JJ));%以天为单位的导温系数；
    HBETA=24*TOM(3,JJ)/(TOM(1,JJ)*TOM(2,JJ));%以天为单位的放热系数；
    SHETA0=WDE(1,JJ);%单位是（C和d);
    a=WDE(2,JJ);
    b=WDE(3,JJ);
    fprintf('MATERIAL PROPERTIES:%d %f %f %f %f\n',JJ,TOM(1,JJ),TOM(2,JJ),TOM(3,JJ),TOM(4,JJ));
end
        