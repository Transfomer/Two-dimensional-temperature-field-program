function [DNX,RJAC]=FDNX(XY,R,S)
%%
%PN-对局部坐标的偏导；
%DNX-对整体坐标的偏导；
[NC,PN,DET,XJAC]=FUN4(XY,R,S);
RJAC=zeros(2,2);
if (DET<1.0E-5)
    fprintf('ERROR:NEGATIVE OR ZERO JACOBIAN DETERMINANT FOR ELEMENT');
else
    RJAC=XJAC^-1;
end
DNX=zeros(2,4);
for I=1:4
    for J=1:2
        DNX(J,I)=0;
        for K=1:2
            DNX(J,I)=DNX(J,I)+RJAC(J,K)*PN(K,I);
        end
    end
end

