function [XR]=COUNTXR(XY)
%%
%XR-单元热容矩阵；
%高斯积分点和权系数
RSTG=zeros(3,1);
RSTG(1)=-0.7745966692414830;
RSTG(2)=0.00;
RSTG(3)=-RSTG(1);
H=zeros(3,1);
H(1)=0.5555555555555560;
H(2)=0.8888888888888890;
H(3)=H(1);
XR=zeros(4,4);
for I=1:4
    for J=1:4
        for K=1:3
            S=RSTG(K);
            SH=H(K);
            for M=1:3
                R=RSTG(M);
                RH=H(M);
                [NC,PN,DET,XJAC]=FUN4(XY,S,R);
                XR(I,J)=XR(I,J)+NC(I)*NC(J)*SH*RH*DET;
            end
        end
    end
end

            
