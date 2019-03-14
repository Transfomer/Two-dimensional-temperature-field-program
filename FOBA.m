function [T1]=FOBA(FF,SK,MA,N)
%%
for I=2:N
    L=I-MA(I)+MA(I-1)+1;
    K=I-1;
    if(L<=K)
        for LP=L:K
            IP=MA(I)-I+LP;
            FF(I)=FF(I)-SK(IP)*FF(LP);
        end
    end
end
for I=1:N
    II=MA(I);
    FF(I)=FF(I)/SK(II);
end
for I=N:-1:2
    L=I-MA(I)+MA(I-1)+1;
    K=I-1;
    if(L<=K)
        for J=L:K
            IJ=MA(I)-I+J;
            FF(J)=FF(J)-SK(IJ)*FF(I);
        end
    end
end
T1=zeros(N,1);
for I=1:N
    T1(I)=FF(I);
end
%fprintf('%f\t',T1);