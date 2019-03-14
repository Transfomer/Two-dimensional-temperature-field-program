function [WEN_RESULT]=WEN(SHETA0,a,b,time)
%%
%求绝热温升对于时间的导数；
if (time==0)
    WEN_RESULT=0;
end
if (time>0)
    WEN_RESULT=SHETA0*a*b*time^(b-1)*exp(-a*time^b);
end