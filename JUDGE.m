function [FT,FF,F0,F1,T1,SK,SK1,SK2]=JUDGE(MEL,COOR,JR,SK1,SK2,MA,T0,SHETA0,a,b,S,F,DT,N,NE,NH,fp,D,ht)
%%
%求解方程组；
time=0;
FF=zeros(N,1);
F0=zeros(N,1);
F1=zeros(N,1);
SK=zeros(1000,1);
fprintf('天数\t   结点号码\t   温度值\n');
fprintf(fp,'天数\t结点号码\t  温度值\r\n');
for M=1:D
    for J=1:N
        FF(J)=0;
    end
    [FT]=CMT(MEL,COOR,JR,SHETA0,a,b,time,N,NE);
    %fprintf('%f\t',FT);
    %fprintf('\n');
    for J=1:N
        F0(J)=FT(J);
    end
    time=time+DT;
    [FT]=CMT(MEL,COOR,JR,SHETA0,a,b,time,N,NE);
    %fprintf('%f\t',FT);
    %fprintf('\n');
    for J=1:N
        F1(J)=FT(J);
    end
    for J=1:N
        for K=1:J
            FF(J)=FF(J)+SK2(MA(J)-J+K)*T0(K);
            %fprintf('%d',T0);
        end
        for K=J+1:N
            FF(J)=FF(J)+SK2(MA(K)-K+J)*T0(K);
        end
        FF(J)=FF(J)+F(J)+(1-S)*F0(J)+S*F1(J);
    end
    %fprintf('%f\t',FF);
    %fprintf('\n');
    [SK]=DECOP(SK,MA,SK1,NH,N);
    %fprintf('%f\t',SK);
    %fprintf('\n');
    [T1]=FOBA(FF,SK,MA,N);
    TMAX=0;
    PP=0;
    for I=1:N
        if(T1(I)>TMAX)
            TMAX=T1(I);
            PP=I;
        end
    end
    if(ht==1)
        fprintf(' %d\t      %d\t      %6.2f\n',M,PP,TMAX);
        fprintf(fp, '%d\t    %d\t  %6.2f\r\n',M,PP,TMAX);
    elseif(ht==2)
        for I=1:N
            fprintf('第%d个：%f\t   ',I,T1(I));
        end
        fprintf('\n');
    end
    iz=1;
    if(iz==1)
        [XY,TF,MEL1]=Renumber_TemperatureNode(MEL',COOR',T1);
        plot_temperature(XY(:,1), XY(:,2), TF, MEL1,M);
        fmat(:,M)=getframe;
    end
    for I=1:N
        T0(I)=T1(I);
    end
end
movie(fmat,10);
    

        