% Code for the paper "Computational adaptive optimal control with an
% application to a car engine control problem",  Yu Jiang and Zhong-Ping
% Jiang,vol. 48, no. 10, pp. 2699-2704, Oct. 2012.
% \copyright Copyright 2011-2014 Yu Jiang, New York University.
function []=engine_main()
clc;
x_save=[];
t_save=[];
flag=1;  % 1: learning is on. 0: learning is off.

% System matrices used for simulation purpose
A =[-20 0 0
    -0.0363 -0.9497 1.0000
    -16.1492 -5.6966 -0.6241]

B=[20
    0
    0]
[xn,un]=size(B);%size of B. un-column #, xn row #
% Set the weighting matrices for the cost function
Q=diag([6.8 8 1.5]);
R=[2];

% Initialize the feedback gain matrix
K=zeros(un,xn); % Only if A is Hurwitz, K can be set as zero.
N=100;  %Length of the window, should be at least greater than xn^2
% NN=10;  %Max iteration times
T=.01;  %Duration of time for each integration
%x0=[10;2;100;2;-1;-2]; %Initial condition
x0=[0.001;0.003;0.003];
i1=(rand(1,100)-.5)*1000;
%i1是一个100个元素的行向量 ，范围（-500，500） rand是一个取（0，1）的随机数
i2=(rand(1,100)-.5)*1000;

Dxx=[];XX=[];XU=[];  % Data matrices
X=[x0;kron(x0',x0')';kron(x0,zeros(un,1))]';
% Run the simulation and obtain the data matrices \delta_{xx}, I_{xx},
% and I_{xu}

for i=1:N
    % Simulation the system and at the same time collect online info.
    [t,X]=ode45(@mysys, [(i-1)*T,i*T],X(end,:));%ode45- 0.01S内采集45次
    %Append new data to the data matrices
    Dxx=[Dxx;kron(X(end,1:xn),X(end,1:xn))-kron(X(1,1:xn),X(1,1:xn))];
    % X的最后一行的第一列到第6列克劳内特积-X的第一行的第一列到第六列克劳内特积
    XX=[XX;X(end,xn+1:xn+xn^2)-X(1,xn+1:xn+xn^2)];
    % X的最后一行的第7列到第42列-X的第一行的第7列到第42列                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
    XU=[XU;X(end,xn+xn^2+1:end)-X(1,xn+xn^2+1:end)];
    %X的最后一行的第43列到最后一列-X的第一行的第43列到最后一列
    % Keep track of the system trajectories
    x_save=[x_save;X(end,1:3)];
    t_save=[t_save;t];
end

Dxx=processing_Dxx(Dxx); % Only the distinct columns left  仅留下不同的列
% K=zeros(un,xn);  % Initial stabilizing feedback gain matrix
P_old=zeros(xn);P=eye(xn)*10 ;% Initialize the previous cost matrix
it=0;
p_save=[];
k_save=[];

% Counter for iterations
% Track the cost matrices in all the iterations
% Track the feedback gain matrix in each iterations
[K0,P0]=lqr(A,B,Q,R) ;
P0% Calculate the ideal solution for comparion purpose
k_save=[norm(K-K0)];

while norm(P-P_old)>1e-10 & it<16   % Stopping criterion for learning
    it=it+1;
    P_old=P;
    QK=Q+K'*R*K;
    X2=XX*kron(eye(xn),K') ; %1*12
    X1=[Dxx,-X2-XU];
    Y=-XX*QK(:);
    pp=X1\Y;
    P=reshape_p(pp)
    p_save=[p_save,norm(P-P0)];
    BPv=pp(end-(xn*un-1):end);
    K=inv(R)*reshape(BPv,un,xn)/2
    k_save=[k_save,norm(K-K0)];
end

P
K

for i=N+1 :800
        [t,X]=ode45(@mysys, [(i-1)*T,i*T],X(end,:));
        x_save=[x_save;X(end,1:3)];
     t_save=[t_save;t];
end
% Plot the trajectories
figure(1)
subplot(2,1,1)
plot([0:length(p_save)-1],p_save,'o',[0:length(p_save)-1],p_save)
axis([0,10,-0.5,5])
legend('||P_k-P^*||')
xlabel('Number of iterations')
% 
% 
% figure(2)
subplot(2,1,2)
plot([0:length(k_save)-1],k_save,'^',[0:length(k_save)-1],k_save)
axis([-0.5,10,-0.5,10])
legend('||K_k-K^*||')
xlabel('Number of iterations')
% % Post-learning simulation
% [tt,xx]=ode23(@mysys,[t(end) 200],X(end,:)');
% % Keep track of the post-learning trajectories
% t_final=[t_save;tt];
% x_final=[x_save;xx];
% 
 figure(2)
 plot([(0:length(x_save)-1)]*0.01,x_save(:,1),'g-');hold on;
 plot([(0:length(x_save)-1)]*0.01,x_save(:,2),'b-');hold on;
 plot([(0:length(x_save)-1)]*0.01,x_save(:,3),'r-');hold on;
legend('x_1','x_2','x_3')
xlabel('Time (sec)')
% 
% figure(4)
% plot(t_final,sqrt(sum(x_final(:,1:6).^2,2)),'Linewidth',2)
% axis([0,200,-50,200])
% legend('||x||')
% xlabel('Time (sec)')
% 
% figure(5)
% plot(t_final,3.6*x_final(:,1),'k-.', ...
%     t_final, x_final(:,6),'-','Linewidth',2)
% axis([0,10,-80,50])
% legend('y_1 (MAF)','y_2 (MAP)')
% xlabel('Time (sec)')

% The following nested function gives the dynamics of the sytem. Also,
% integraters are included for the purpose of data collection.
    function dX=mysys(t,X)
        %global A B xn un i1 i2 K flag
        x=X(1:xn);
        if t>=2;   % See if learning is stopped     N=200 T=.01s
            flag=0;
        end
        if flag==1
            u=zeros(un,1);
            for i=i1
                u(1)=u(1)+0.001*sin(i*t)/length(i1); % constructing the
                % exploration noise
            end
           
            u=10000*u;
        else
            u=-K*x;
        end
        dx=A*x+B*u;
        dxx=kron(x',x')';
        dux=kron(x',u')';
        dX=[dx;dxx;dux];
    end
% This nested function reconstruct the P matrix from its distinct elements
% 此嵌套函数根据其不同元素重建P矩阵
    function P=reshape_p(p)
        P=zeros(xn);
        ij=0;
        for i=1:xn
            
            for j=1:i
                ij=ij+1;
                P(i,j)=p(ij);
                P(j,i)=P(i,j);
            end
        end
    end

% The following nested function removes the repeated columns from Dxx

    function Dxx=processing_Dxx(Dxx)
        ij=[]; ii=[];
        for i=1:xn
            ii=[ii (i-1)*xn+i];
        end
        for i=1:xn-1
            for j=i+1:xn
                ij=[ij (i-1)*xn+j];
            end
        end
        Dxx(:,ii)=Dxx(:,ii)/2;
        Dxx(:,ij)=[];
        Dxx=Dxx*2;
    end
end





