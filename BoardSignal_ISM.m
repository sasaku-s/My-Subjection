clc
clear
close
M=12;
N=200;
ts=0.01;
f0=100;
f1=80;
f2=120;
c=1500;
lambda=c/f0;
d=lambda/2;
SNR=15;
b=pi/180;
theat1=30*b;%%%信号1
theat2=25*b;%%%信号2
n=ts:ts:N*ts;
theat=[theat1 theat2]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%信号处理%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s1=chirp(n,80,1,120);%%线性调频1
sa=fft(s1,2048);
figure
specgram(s1,256,1E3,256,250);%频谱
s2=chirp(n+0.100,80,1,120);%%线性调频2
sb=fft(s2,2048);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%ISM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P=1:2;
a=zeros(M,2);
sump=zeros(1,181);
for i=1:N
    f=80+(i-1)*1.0;
    s=[sa(i) sb(i)]';
    for m=1:M
        a(m,P)=exp(-1i*2*pi*d/c*sin(theat(P))*(m-1))';
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
pmusic=1/33*sump;
pm=1./pmusic;
thetaesti=-90:1:90;
plot(thetaesti,20*log(abs(pm)));
xlabel('入射角/度');
ylabel('空间谱/db');
grid on