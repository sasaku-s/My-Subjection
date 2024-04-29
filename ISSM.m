clc
clear all
close all
%%%%%%%%%%%%%%%% 阵列模型 %%%%%%%%%%%%%
M=12;                                     %阵元数
N=300;                                    %快拍数
ts=0.00125;                               %时域采样间隔
f0=100;                                   %入射信号中心频率
f1=80;                                    %入射信号最低频率
f2=120;                                   %入射信号最高频率
c=3e8;                                    %波速
lambda=c/f0;                              %波长
d=lambda/2;                               %阵元间距
SNR=12;                                   %信噪比
b=pi/180;
theat1=0*b;                               %入射信号波束角1
theat2=40*b;                              %入射信号波束角2
theat3=45*b;
theat4=-25*b;
n=ts:ts:N*ts;
theat=[theat1 theat2 theat3 theat4]';

%%%%%%%%%%%%%%%% 信号处理 %%%%%%%%%%%%%%%%
s1=chirp(n,80,1,120);                      %线性调制
sa=fft(s1,512);                           %FFT变换

s2=chirp(n+0.125,80,1,120);                
sb=fft(s2,512);                           

s3=chirp(n+0.25,80,1,120);
sc=fft(s3,512);

s4=chirp(n+0.375,80,1,120);
sd=fft(s4,512);

%%%%%%%%%%%%%%%%%%%%% ISM算法 %%%%%%%%%%%%%%%%%%
P=1:4;%%%信号数
a=zeros(M,4);
sump=zeros(1,181);
for i=1:N
    f=80+(i-1)*1.25;
    s=[sa(i) sb(i) sc(i) sd(i)]';%%接收信号矩阵
    for m=1:M
        a(m,P)=exp(-1i*2*pi*f*d/c*sin(theat(P))*(m-1))';
    end
    R=a*(s*s')*a';
    [em,zm]=eig(R);
    [zm1,pos1]=max(zm);
    for l=1:2
        [zm2,pos2]=max(zm1);
        zm1(:,pos2)=[];
        em(:,pos2)=[];
    end
    k=1;
    for ii=-90:1:90
        arfa=sin(ii*b)*d/c;
        for iii=1:M
            tao(1,iii)=(iii-1)*arfa;
        end
        A=[exp(-1i*2*pi*f*tao)]';
        p(k)=A'*em*em'*A;
        k=k+1;
    end
    sump=sump+abs(p);
end
pmusic=(1/N)*sump;
pm=1./pmusic;
thetaesti=-90:1:90;
plot(thetaesti,20*log(abs(pm)));
xlabel('入射角/度');
ylabel('空间谱/dB');
grid on
