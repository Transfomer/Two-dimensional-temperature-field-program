function [SK1,SK2]=SKO(iew,JR,SK1,SK2,MEL,COOR,MA,WG,NSE,HBETA,S)
%%
%求方程的K,将G叠加到K上面（放热边界对热传导矩阵的贡献矩阵）；
RSTG=zeros(2,1);
RSTG(1)=-0.5773502691896260;
RSTG(2)=-RSTG(1);
H=zeros(2,1);
H(1)=1.0;
H(2)=1.0;
NN=zeros(4,1);
XY=zeros(2,4);
RST=zeros(2,1);
RSH=zeros(2,1);
NODES=zeros(10,1);
KCRD=[1 1 2 2];
KFACE=[1 2 1 4;4 3 2 3];
FVAL=[-1 1 -1 1];
XG=zeros(4,4);
%%
for IE=1:NSE
    NEE=iew(IE);
    for I=1:4
        for J=1:4
            XG(I,J)=0;
        end
    end
    for K=1:4
        IEK=MEL(K,NEE);
        NN(K)=JR(IEK);
        for M=1:2
            XY(M,K)=COOR(M,IEK);
        end
    end
        ML=KCRD(WG);
        %fprintf('%d\n',ML);
        if (ML==1)
            MM=2;
        elseif(ML==2)
            MM=1;
        end
        RST(ML)=FVAL(WG);
        for I=1:2
            NODES(I)=KFACE(I,WG);
            M=NODES(I);
            for J=1:2
                NODES(J)=KFACE(J,WG);
                T=NODES(J);
                for K=1:2
                    RST(MM)=RSTG(K);
                    RSH(MM)=H(K);
                    [NC,PN,DET,XJAC]=FUN4(XY,RST(1),RST(2));
                    LX=sqrt(XJAC(MM,1)^2+XJAC(MM,2)^2);
                     %fprintf('%f\t',LX);
                    XG(M,T)=XG(M,T)+HBETA*RSH(MM)*NC(M)*NC(T)*LX;
                    %fprintf('%f\t',XG(M,T));
                    %fprintf('\n');
                end
                II=NN(M);
                JJ=NN(T);
                if (JJ>0)&&(II>=JJ)
                    SK1(MA(II)-II+JJ)=SK1(MA(II)-II+JJ)+S*XG(M,T);
                    SK2(MA(II)-II+JJ)=SK2(MA(II)-II+JJ)-(1-S)*XG(M,T);
                end
            end
        end
end
