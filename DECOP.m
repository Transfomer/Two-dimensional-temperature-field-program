function [SK]=DECOP(SK,MA,SK1,NH,N)
%%
for I=1:NH
    SK(I)=SK1(I);
end
for I=2:N
    L=I-MA(I)+MA(I-1)+1;
    L1=L+1;
    K=I-1;
    if(L1<=K)
        for J=L1:K
            IJ=MA(I)-I+J;
            M=J-MA(J)+MA(J-1)+1;
            if(L>M)
                M=L; 
            end
            MP=J-1;
            if(M<=MP)
                for LP=M:MP
                    IP=MA(I)-I+LP;
                    JP=MA(J)-J+LP;
                    SK(IJ)=SK(IJ)-SK(IP)*SK(JP);
                end
            end
        end
    end
    if(L<=K)
        for LP=L:K
            IP=MA(I)-I+LP;
            LPP=MA(LP);
            SK(IP)=SK(IP)/SK(LPP);
            II=MA(I);
            SK(II)=SK(II)-SK(IP)*SK(IP)*SK(LPP);
        end
    end
end
            
            