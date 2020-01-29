function [S]   =  Periodical_Simulated1(Ori,Perio,rate,mean)
%% 
%case1:
%%%%%%%%%%%%%%%%%%%%-----------Chang creat stripe noise---------%%%%%%
%��1����10�������������ʼ�У�Ȼ����10Ϊ���ڼ���������
%ǰ�����ͬʱ����һ���̶�������mean
%�������ͬʱ��ȥһ���̶�������mean
%������Ч���Ϳ��Դﵽ�ӵ���������Ϊ1
rand('seed',2);
Location = randperm(Perio,floor(rate*Perio));
[Row, Col] = size(Ori);Row=2*Row;Col=2*Col;
S=zeros(Row,Col);
for i=0:1:floor(Col/Perio)-1   
    S(:,Location(1:floor(length(Location)/2)) + i * Perio) = S(:,Location(1:floor(length(Location)/2))...
      + i * Perio) + mean;
    S(:,Location(floor(length(Location)/2)+1:length(Location))+ i * Perio) = ...
    S(:,Location(floor(length(Location)/2)+1:length(Location))+ i * Perio) - mean;
end

%%
%case2  ��   �����ϵ�ֵ��ͬ�����ǲ�ͬ�����ϵ�ֵ��ͬ(�����и���������Ҳ���Դﵽ��������Ϊ1
% rand('seed',3);
% Location = randperm(Perio,floor(rate*Perio));
% [Row, Col] = size(Ori);
% u= randi([0 50],2,floor(Col/Perio));
% S=zeros(Row,Col);
% for i=0:1:floor(Col/Perio)-1   
%     S(:,Location(1:floor(length(Location)/2)) + i * Perio) = S(:,Location(1:floor(length(Location)/2))...
%     + i * Perio) + u(1,i+1);
%     S(:,Location(floor(length(Location)/2)+1:length(Location))+ i * Perio) = ...
%     S(:,Location(floor(length(Location)/2)+1:length(Location))+ i * Perio) - u(2,i+1);   %+ u(2,i+1)����ֵȫΪ����
% end
% P_Stripe=Ori+S;

%%
%case3  ��   �����ϵ�ֵ��ͬ�����Ҳ�ͬ�����ϵ�ֵ��ͬ(�����и���������������ʹ����������Ϊ1
% rand('seed',0);
% Location = randperm(Perio,floor(rate*Perio));
% [Row, Col] = size(Ori);
% u= randi([0 60],Row,2*floor(Col/Perio));
% S=zeros(Row,Col);
% for i=0:1:floor(Col/Perio)-1   
%     S(:,Location(1:floor(length(Location)/2)) + i * Perio) = S(:,Location(1:floor(length(Location)/2))...
%     + i * Perio) + repmat(u(:,i+1),1,floor(length(Location)/2));
%     num=length(Location)-floor(length(Location)/2);
%     S(:,Location(floor(length(Location)/2)+1:length(Location))+ i * Perio) = ...
%     S(:,Location(floor(length(Location)/2)+1:length(Location))+ i * Perio) - repmat(u(:,floor(Col/Perio)+i+1),1,num);   %+ u(2,i+1)����ֵȫΪ����
% end
% P_Stripe=Ori+S;
%%
%case4:    ���������ɷֵ��ȵ���
% rand('seed',0);
% Location = randperm(Perio,floor(rate*Perio));
% [Row, Col] = size(Ori);
% u= randi([0 60],2,floor(Col/Perio));
% S=zeros(Row,Col);
% for i=0:1:floor(Col/Perio)-1   
%     S(:,Location(1:floor(length(Location)/2)) + i * Perio) = S(:,Location(1:floor(length(Location)/2))...
%     + i * Perio) + u(1,i+1);
%     S(:,Location(floor(length(Location)/2)+1:length(Location))+ i * Perio) = ...
%     S(:,Location(floor(length(Location)/2)+1:length(Location))+ i * Perio) - u(2,i+1);   %+ u(2,i+1)����ֵȫΪ����
% end
% Location = randperm(Row*Col,1000);
% S(Location)=randi([0 60],1,1000);
% P_Stripe=Ori+S;
%%
% rand('seed',2);
% % Location = randperm(Perio,floor(rate*Perio));
% [Row, Col] = size(Ori);
% S=zeros(Row,Col);
% for i=0:1:floor(Col/Perio)-1   
%     S(:,4 + i * Perio) = S(:,4 ...
%       + i * Perio) + mean;
%     S(:,5+ i * Perio) = ...
%     S(:,5+ i * Perio) - mean;
% end
% P_Stripe=Ori+S;
%%
%case2  ��   �����ϵ�ֵ��ͬ�����ǲ�ͬ�����ϵ�ֵ��ͬ(�����и���������Ҳ���Դﵽ��������Ϊ1
% rand('seed',2);
% [Row, Col] = size(Ori);
% u= randi([0 50],2,floor(Col/Perio));
% S=zeros(Row,Col);
% for i=0:1:floor(Col/Perio)-1   
%     S(:,4 + i * Perio) = S(:,4 ...
%     + i * Perio) + u(1,i+1);
%     S(:,5+ i * Perio) = ...
%     S(:,5+ i * Perio) - u(2,i+1);   %+ u(2,i+1)����ֵȫΪ����
% end
% P_Stripe=Ori+S;