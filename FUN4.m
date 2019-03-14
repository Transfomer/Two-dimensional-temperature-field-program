function [NC,PN,DET,XJAC]=FUN4(XY,R,S)
%%
%�����κ����Լ��Ծֲ������ƫ������
JBZB=[-1.0 1.0 1.0 -1.0;-1.0 -1.0 1.0 1.0];
NC=zeros(4,1);
PN=zeros(2,4);
for I=1:4
    NC(I)=0.25*(1+R*JBZB(1,I))*(1+S*JBZB(2,I));
    PN(1,I)=0.25*JBZB(1,I)*(1+S*JBZB(2,I));
    PN(2,I)=0.25*JBZB(2,I)*(1+R*JBZB(1,I));
end
%%
%�����ſɱȾ����Լ�������ʽ��ֵ��
XJAC=zeros(2,2);
for I=1:2
    for J=1:2
        XJAC(I,J)=0;
        for K=1:4
            XJAC(I,J)=XJAC(I,J)+PN(I,K)*XY(J,K);
        end
    end
end
DET=det(XJAC);
%fprintf('%f\t',DET);
    