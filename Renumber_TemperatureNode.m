function [XY,TF,MEL]=Renumber_TemperatureNode(Mel_PhP,xyPhP,T1)
global EleTopo nnel 
n=size(Mel_PhP,1);
P=zeros(4,1);
TT=zeros(4,1);
x=zeros(4,1);
y=zeros(4,1);
xyEle=zeros(4,2);
e=0;
TF=zeros(n*121,1);
MEL=zeros(n*100,4);
for I=1:n
    for J=1:4
        P(J)=Mel_PhP(I,J);
        x(J)=xyPhP(P(J),1);
        y(J)=xyPhP(P(J),2);
        xyEle(J,1:2)=xyPhP(P(J),1:2);
        TT(J)=T1(P(J));
    end
    for K=0:10
        for L=1:11
            t=L+K*11+(I-1)*121;
            XY(t,1)=x(1)+K*(x(2)-x(1))/10;
            XY(t,2)=y(1)+(L-1)*(y(4)-y(1))/10;
            xyGP(1)=XY(t,1);
            xyGP(2)=XY(t,2);
            [xylocGP] = LocCoordFromGloCoordIsoPara(xyGP,xyEle,nnel,2,EleTopo);%确定物理边界在数学网格的哪个点；
            %fprintf('%f\t',xylocGP);
            xlocGP = xylocGP(1);  ylocGP = xylocGP(2); 
            [fw, dwdr, dwds] = get_shape_QUAD(xlocGP, ylocGP,nnel,0); %求物理边界上的点在对应的数学网格上的形函数；
            %fprintf('%f\t',TT);
            for M=1:4
                TF(t,1)=TF(t,1)+fw(M)*TT(M);
            end
        end
    end
    for P=0:9
        for H=1:10
            e=e+1;
            MEL(e,1)=H+P*11+(I-1)*121;
            MEL(e,4)=H+P*11+1+(I-1)*121;
            MEL(e,2)=H+(P+1)*11+(I-1)*121;
            MEL(e,3)=H+(P+1)*11+(I-1)*121+1;
        end
    end
end

    
    
    
        
    
    