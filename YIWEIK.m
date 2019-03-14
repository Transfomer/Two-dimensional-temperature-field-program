function [MA,NH]=YIWEIK(fp,JR,MEL,N,NE)
%%
%求MA;
MA=zeros(N,1);
NN=zeros(4,1);
for I=1:N
    MA(I)=0;
end
for IE=1:NE
    for K=1:4
        IEK=MEL(K,IE);
        NN(K)=JR(IEK);
    end
    L=N;
    %找到最小半带宽；
    for I=1:4
        if(NN(I)>0)
            if (NN(I)<L)
                L=NN(I);
            end
        end
    end
    for I=1:4
        JP=NN(I);
        if(JP>0)
            MA(JP)=JP;
        end
    end
end
for I=2:N
    MA(I)=MA(I)+MA(I-1);
end
NH=MA(N);
fprintf('FREEDOM N=%d\t STORAGE NH=%d\n',N,NH);
fprintf(fp,'FREEDOM N=%d\t STORAGE NH=%d\r\n',N,NH);
            