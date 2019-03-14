function [FT]=CMT(MEL,COOR,JR,SHETA0,a,b,time,N,NE)
%%
%绝热温升的温度叠加到结点上；
RSTG=zeros(3,1);
RSTG(1)=-0.7745966692414830;
RSTG(2)=0.00;
RSTG(3)=-RSTG(1);
H=zeros(3,1);
H(1)=0.5555555555555560;
H(2)=0.8888888888888890;
H(3)=H(1);
FT=zeros(N,1);
XF1=zeros(4,1);
NN=zeros(4,1);
XY=zeros(2,4);
for IE=1:NE
    NEE=IE;
    [WEN_RESULT]=WEN(SHETA0,a,b,time);
    for I=1:4
        XF1(I)=0;
    end   
    for K=1:4
        IEK=MEL(K,NEE);
        NN(K)=JR(IEK);
        for M=1:2
            XY(M,K)=COOR(M,IEK);
        end
    end
    for I=1:4
        J=NN(I);
        for K=1:3
            SS=RSTG(K);
            SH=H(K);
            for M=1:3
                R=RSTG(M);
                RH=H(M);
                [NC,PN,DET,XJAC]=FUN4(XY,R,SS);
                %fprintf('%f\t%d\n',WEN_RESULT,time);
                %fprintf('%f\t',SH);
                XF1(I)=XF1(I)+WEN_RESULT*NC(I)*SH*RH*DET;
                %fprintf('%f\t',XF1(I));
            end
        end
        if (J>0)
            FT(J)=FT(J)+XF1(I);
            %fprintf('%f\t',FT(J));
        end
    end
end
        
    