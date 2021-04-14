function [xbest]=BAS_Algorithm_up(y,smp)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%初始化部分

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eta=0.95;

c=5;      %ratio between step and d0

step=1;   %initial step set as the

% largest input range

n=100;%iterations

%k=20;%space dimension

grid_numN =96; %网格点数
f_grid=0:1/grid_numN:1-1/grid_numN;
f=f_grid'; %初始频率网络
z=ones(grid_numN,1); %初始 t+1次 估计矩阵

%x=rands(k,1);%intial value

x=zeros(grid_numN,2);
x(:,1)=f;
x(:,2)=z;

xbest=x;

fbest=func(xbest,y,smp);

fbest_store=fbest;

f_store=[0;x(:,1);fbest];

display(['0:','xbest1=[',num2str(xbest(:,1)'),'],fbest1=',num2str(fbest)])
display(['0:','xbest2=[',num2str(xbest(:,2)'),'],fbest2=',num2str(fbest)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%迭代部分

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:n

d0=step/c;

dir0=rands(grid_numN,1);

dir0=dir0/(eps+norm(dir0));

dir=zeros(grid_numN,2);
dir(:,1)=dir0;
dir(:,2)=dir0;

xleft=x+dir*d0;

fleft=func(xleft,y,smp);

xright=x-dir*d0;

fright=func(xright,y,smp);

x=x-step*dir*sign(fleft-fright);

func=func(x,y,smp);

%%%%%%%%%%%

if func<fbest

xbest=x;

fbest=func;

end

%%%%%%%%%%%

f_store=cat(2,f_store,[i;x(:,1);func]);

fbest_store=[fbest_store;fbest];

display([num2str(i),':xbest1=[',num2str(xbest(:,1)'),'],fbest1=',num2str(fbest)])
display([num2str(i),':xbest2=[',num2str(xbest(:,2)'),'],fbest2=',num2str(fbest)])
%%%%%%%%%%%

step=step*eta;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%数据显示部分

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1),clf(1),

plot(f_store(1,:),f_store(end,:),'r-o')

hold on,

plot(f_store(1,:),fbest_store,'b-.')

xlabel('iteration')

ylabel('minimum value')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%被优化的函数，这部分需要换用你自己的被优化函数

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function func=func(x,y,smp)

gamma=1;
dd=1./(abs(x(:,2)).^2+gamma); 
D=diag(dd); %生成预测对角阵

A=exp(-2i*pi*smp'*x(:,2)');
L=0;
for j=1:length(x(:,2))
    L=L+log(abs(x(j,2))^2+gamma);
end
func=L+(norm(y-A*x(:,2)))^2;
end
