%% PROGRAM TTF
clear;clc;close all;format long;
%% 指定全局变量；
global NP NE NM NR %结点数;单元数;材料数;已知温度数；
global NH N %按一维存储的容量;自由度总数;
%global NEE NME %材料号;
%% 加载路径;
addpath(fullfile(pwd,'\postplot\'));
addpath(fullfile(pwd,'\calculate\'));
%% 定义变量的大小;
KPT=zeros(2,5000);%已知温度结点矩阵;
XH=zeros(4,4);%h;
XG=zeros(4,4);%g;
XY=zeros(2,4);
XR=zeros(4,4);%r;
%% 主程序；
filepath=strcat(pwd,'\prep\');
fileoutpath=strcat(pwd,'\post\');
%这里可以修改输入文件的名称；
FILENAME='INPUTFILE2.txt';
%这里可以修改输出文件的名称；
OUTPUT='OUT1.txt';
fp=fopen(strcat(fileoutpath,OUTPUT),'w');
fclose(fp);
fp=fopen(strcat(fileoutpath,OUTPUT),'a');
ft=fopen(strcat(filepath,FILENAME),'r');
C=fscanf(ft,'%d %d %d %d',4);
NP=C(1);NE=C(2);NM=C(3);NR=C(4);
fprintf('NP=%d\t NE=%d\t NM=%d\t NR=%d\n',NP,NE,NM,NR);
fprintf(fp,'NP=%d\t NE=%d\t NM=%d\t NR=%d\r\n',NP,NE,NM,NR);
[AFA,HBETA,SHETA0,a,b,JR,COOR,MEL,TOM,WDE,N]=INPUT(ft,NP,NE,NM); %Ta-周围介质温度;AFA-导温系数;HBETA-放热系数；
%IT-已知一结点的温度；S-差分参数；DT-时间步长；MEL-单元信息;COOR-整体结点坐标矩阵;TOM-材料常数;WDE-绝热温升参数;
[MA,NH]=YIWEIK(fp,JR,MEL,N,NE);%JR-结点自由度序号矩阵;形成MA
C=fscanf(ft,'%f %f %d %f',4);
S=C(1);DT=C(2);D=C(3);ST=C(4);
fprintf('S=%f\t DT=%f\t D=%d\t ST=%f\n',S,DT,D,ST);
fprintf(fp,'S=%4.1f\t DT=%4.1f\t D=%d\t ST=%4.1f\r\n',S,DT,D,ST);
C=fscanf(ft,'%f %f',2);
Ta=C(1);iz=C(2);
fprintf('Ta=%d\t iz=%d\n',Ta,iz);%iz-散热边界批号;
fprintf(fp,'Ta=%d\t iz=%d\r\n',Ta,iz);
[SK1,SK2]=ARE(JR,MEL,COOR,DT,NE,MA);
[SK1,SK2]=SHO(SK1,SK2,JR,COOR,MEL,AFA,S,NE,MA);%形成H ;
F=zeros(N,1);
if(iz>0)
    for JJ=1:iz
        C=fscanf(ft,'%d %d %d',3);%JS-散热边界批号;NSE-第JS批放热边界的单元个数;WG-散热边界单元面号;
        JS=C(1);NSE=C(2);WG=C(3);
        iew=zeros(NSE,1);
        for M=1:NSE
            iew(M)=fscanf(ft,'%d',1);
        end
        [SK1,SK2]=SKO(iew,JR,SK1,SK2,MEL,COOR,MA,WG,NSE,HBETA,S);%循环,将G累加到H上
        [F]=CMX(F,iew,JR,NSE,MEL,COOR,WG,HBETA,Ta);%形成F;
        %fprintf('%f\t',SK2);
    end
end
T0=zeros(N,1);
for I=1:N
    T0(I)=ST;
end
ht=1;
[FT,FF,F0,F1,T1,SK,SK1,SK2]=JUDGE(MEL,COOR,JR,SK1,SK2,MA,T0,SHETA0,a,b,S,F,DT,N,NE,NH,fp,D,ht);
%%
fprintf('PROGRAM HAS BEEN ENDEN');
fprintf(fp,'PROGRAM HAS BEEN ENDEN');
fclose(fp);
%% plot单元的结点和网格划分情况（figure1)；
tt=0;
if(tt==1);
    %% 求计算区域网格划分在x和y的份数；
    for I=1:NP
        x1=max(COOR(1,:));
        x2=min(COOR(1,:));
        y1=max(COOR(2,:));
        y2=min(COOR(1,:));
        if (COOR(1,I)==x1)&&(COOR(2,I)==y2)
            t1=I;
        end
        if (COOR(2,I)==y1)&&(COOR(1,I)==x1)
            t2=I/t1;
        end
    end    
    %% 求计算区域的边界的结点号；
    [C,n1]=fscanf(ft,'%d',t1);
    BC_S(:)=C(:);
    IEN_BC_s=zeros(2,n1-1);
    for I=1:n1-1
        IEN_BC_s(1,I)=BC_S(I);
        IEN_BC_s(2,I)=BC_S(I+1);
    end
    [C]=fscanf(ft,'%d',t2);
    BC_E(:)=C(:);
    IEN_BC_e=zeros(2,n1-1);
    for I=1:n1-1
        IEN_BC_e(1,I)=BC_E(I);
        IEN_BC_e(2,I)=BC_E(I+1);
    end
    [C]=fscanf(ft,'%d',t1);
    BC_N(:)=C(:);
    IEN_BC_n=zeros(2,n1-1);
    for I=1:n1-1
        IEN_BC_n(1,I)=BC_N(I);
        IEN_BC_n(2,I)=BC_N(I+1);
    end
    [C]=fscanf(ft,'%d',t2);
    BC_W(:)=C(:);
    IEN_BC_w=zeros(2,n1-1);
    for I=1:n1-1
        IEN_BC_w(1,I)=BC_W(I);
        IEN_BC_w(2,I)=BC_W(I+1);
    end
    [fig,ax] = init_plot_figure();
    Xmin=min(COOR(1,:))-0.1;
    Xmax=max(COOR(1,:))+0.1;
    Ymin=min(COOR(2,:))-0.1;
    Ymax=max(COOR(2,:))+0.1;
    ax.XLim = [Xmin, Xmax];
    ax.YLim = [Ymin, Ymax];
    %plot结点号码；
    x=COOR(1,:)';
    y=COOR(2,:)';
    z=zeros(NP,1);
    plot_node_labels(ax, 1:NP, x,y,z);
    IEN=MEL(1:4,:);
    eltype(1:NE)=3;
    %plot单元的透明度，也就是给单元着色；
    plot_element(ax, 1:NE, IEN, eltype, x,y,z, 'surf_alpha', 0.2);
    %plot单元的号码；
    plot_element_labels(ax, IEN, eltype, x, y, z);
    %% 区分四个方向边界的结点；
    % Plot the BC nodes
    %南方向结点颜色为红；
    plot_node(ax, unique(IEN_BC_s(:)), x,y,z, 'r');
    %东方向结点颜色为绿；
    plot_node(ax, unique(IEN_BC_e(:)), x,y,z, 'g');
    %北方向结点颜色为蓝；
    plot_node(ax, unique(IEN_BC_n(:)), x,y,z, 'b');
    %西方向结点颜色为品红；
    plot_node(ax, unique(IEN_BC_w(:)), x,y,z, 'm');
    %%
    count=0;
    neq=0;
    ID=zeros(NP,1);
    ned=1;
    for  I = 1:NP
        for J = 1:ned
            % find Essential BC
            if( fix(I) == 1 )
                ID(J,I) = NP - count;
                count = count + 1;
            else
                % list no gg's
                neq = neq+1;
                ID(J,I) = neq;
            end
        end
    end
    %% plot计算结果的云图；
    [~,ax2] = init_plot_figure();
    Xmin=min(COOR(1,:))-0.1;
    Xmax=max(COOR(1,:))+0.1;
    Ymin=min(COOR(2,:))-0.1;
    Ymax=max(COOR(2,:))+0.1;
    ax2.XLim = [Xmin, Xmax];
    ax2.YLim = [Ymin, Ymax];
    plot_element_solution(ax2, IEN, ID, eltype, x, y, z, T1);
end


    

