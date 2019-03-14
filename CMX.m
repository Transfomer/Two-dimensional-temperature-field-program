function [F]=CMX(F,iew,JR,NSE,MEL,COOR,WG,HBETA,Ta)
%%
%求方程荷载由于放热边界引起的一部分；
RSTG=zeros(3,1);
RSTG(1)=-0.7745966692414830;
RSTG(2)=0.00;
RSTG(3)=-RSTG(1);
H=zeros(3,1);
H(1)=0.5555555555555560;
H(2)=0.8888888888888890;
H(3)=H(1);
KCRD=[1 1 2 2];
KFACE=[1 2 1 4;4 3 2 3];
FVAL=[-1 1 -1 1];
NODES=zeros(2,1);
XF=zeros(4,1);
NN=zeros(4,1);
XY=zeros(2,4);
RST=zeros(2,1);
RSH=zeros(2,1);
%%
for IE=1:NSE
    NEE=iew(IE);
     for I=1:4
        XF(I)=0;
    end 
    for K=1:4
        IEK=MEL(K,NEE);
        NN(K)=JR(IEK);
        for M=1:2
            XY(M,K)=COOR(M,IEK);
        end
    end
    ML=KCRD(WG);
    if (ML==1)
        MM=2;
    elseif(ML==2)
        MM=1;
    end
    RST(ML)=FVAL(WG);
    for I=1:2
        NODES(I)=KFACE(I,WG);
        M=NODES(I);
        P=NN(M);
        for K=1:3
            RST(MM)=RSTG(K);
            RSH(MM)=H(K);
            [NC,PN,DET,XJAC]=FUN4(XY,RST(1),RST(2));
            LX=sqrt(XJAC(MM,1)^2+XJAC(MM,2)^2);
            XF(M)=XF(M)+HBETA*Ta*RSH(MM)*NC(M)*LX;
        end
        if (P>0)
            F(P)=F(P)+XF(M);
            %fprintf('%f\t', F(P));
        end
    end
end
            
    