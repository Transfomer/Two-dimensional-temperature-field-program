function [AFA,HBETA,SHETA0,a,b,JR,COOR,MEL,TOM,WDE,N]=INPUT(ft,NP,NE,NM)
%%
%��ȡ������ꣻ
COOR=zeros(2,NP);
for I=1:NP
    C=fscanf(ft,'%d %f %f',3);
    IP=C(1);COOR(1,IP)=C(2);COOR(2,IP)=C(3);
end
%%
%��ȡÿ����Ԫ�����Ľ����룻
MEL=zeros(5,NE);
for I=1:NE
    C=fscanf(ft,'%d %d %d %d %d %d',6);
    NEE=C(1);NME=C(2);MEL(1,NEE)=C(3);MEL(2,NEE)=C(4);MEL(3,NEE)=C(5);MEL(4,NEE)=C(6);
    MEL(5,NEE)=NME;
end
%%
%����ܵ����ɶȣ�
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
%��ȡ���ϲ����;��������Ĳ�����
TOM=zeros(4,NM);
WDE=zeros(3,NM);
for J=1:NM
    C=fscanf(ft,'%d %f %f %f %f %f %f %f',8);
    JJ=C(1);TOM(1,JJ)=C(2);TOM(2,JJ)=C(3);TOM(3,JJ)=C(4);TOM(4,JJ)=C(5);WDE(1,JJ)=C(6);WDE(2,JJ)=C(7);WDE(3,JJ)=C(8);
    AFA=24*TOM(4,JJ)/(TOM(1,JJ)*TOM(2,JJ));%����Ϊ��λ�ĵ���ϵ����
    HBETA=24*TOM(3,JJ)/(TOM(1,JJ)*TOM(2,JJ));%����Ϊ��λ�ķ���ϵ����
    SHETA0=WDE(1,JJ);%��λ�ǣ�C��d);
    a=WDE(2,JJ);
    b=WDE(3,JJ);
    fprintf('MATERIAL PROPERTIES:%d %f %f %f %f\n',JJ,TOM(1,JJ),TOM(2,JJ),TOM(3,JJ),TOM(4,JJ));
end
        