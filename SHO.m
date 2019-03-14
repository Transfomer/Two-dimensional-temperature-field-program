function [SK1,SK2]=SHO(SK1,SK2,JR,COOR,MEL,AFA,S,NE,MA)
%%
%在方程组的两边形成用差分的解法形成的关于热传导矩阵的部分；
NN=zeros(4,1);
XY=zeros(2,4);
for IE=1:NE
    NEE=IE;
    for K=1:4
        IEK=MEL(K,NEE);
        NN(K)=JR(IEK);
        for M=1:2
            XY(M,K)=COOR(M,IEK);
        end
    end
    [XH]=COUNTXH(XY,AFA);
    for I=1:4
        for J=1:4
            II=NN(I);
            JJ=NN(J);
            if (JJ>0) && (II>=JJ)
                SK1(MA(II)-II+JJ)=SK1(MA(II)-II+JJ)+S*XH(I,J);
                SK2(MA(II)-II+JJ)=SK2(MA(II)-II+JJ)-(1-S)*XH(I,J);
            end
        end
    end
end
