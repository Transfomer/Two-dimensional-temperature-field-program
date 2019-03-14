function [XH]=COUNTXH(XY,AFA)
%%
%XH-单元热传导矩阵；
RSTG=zeros(3,1);
RSTG(1)=-0.7745966692414830;
RSTG(2)=0.00;
RSTG(3)=-RSTG(1);
H=zeros(3,1);
H(1)=0.5555555555555560;
H(2)=0.8888888888888890;
H(3)=H(1);
XH=zeros(4,4);
%%
NI=zeros(4,4);
for I=1:4
    for J=1:4
        for K=1:3
            S=RSTG(K);
            SH=H(K);
            for M=1:3
                R=RSTG(M);
                RH=H(M);
                [DNX,RJAC]=FDNX(XY,R,S);
                [NC,PN,DET,XJAC]=FUN4(XY,R,S);
                NI(I,J)=DNX(1,I)*DNX(1,J)+DNX(2,I)*DNX(2,J);
                XH(I,J)=XH(I,J)+AFA*NI(I,J)*SH*RH*DET;
            end
        end
    end
end
            
