function [NonP_Stripe,S]   =  NonPeriodical_Simulated0_S(Ori,rate,mean)
%% 
%case1:
%%%%%%%%%%%%%%%%%%%%-----------Chang creat stripe noise---------%%%%%%
%��1����10�������������ʼ�У�Ȼ����10Ϊ���ڼ���������
%ǰ�����ͬʱ����һ���̶�������mean
%�������ͬʱ��ȥһ���̶�������mean
%������Ч���Ϳ��Դﵽ�ӵ���������Ϊ1
% rand('seed',1);
% [Row, Col] = size(Ori);
% S=zeros(Row,Col);
% Location = randperm(Col,round(rate*Col));
% S(:,Location(1:round(rate*Col/2)))=S(:,Location(1:round(rate*Col/2))) + mean;
% S(:,Location(round(rate*Col/2)+1:round(rate*Col)))=...
%            S(:,Location(round(rate*Col/2)+1:round(rate*Col))) - mean;
% NonP_Stripe=Ori+S;
 %%
 %case2  ��   �����ϵ�ֵ��ͬ�����ǲ�ͬ�����ϵ�ֵ��ͬ(�����и���������Ҳ���Դﵽ��������Ϊ1
 rand('seed',1);
[Row, Col] = size(Ori);
Location1 = randperm(Col,round(rate*Col));%���������Ҫ���������е�λ��
%u=randi([0 100],1,length(Location1));  %�������һ����Location������ͬ��0��100֮������������



S=zeros(Row,Col);
% Location(1:round(rate*Col/2)) ȡlocation��ǰһ�� 
%Location(round(rate*Col/2)+1:round(rate*Col)) ȡLocation�ĺ�һ��
S(1:100,Location1(1:round(rate*Col/2)))= ...
S(1:100,Location1(1:round(rate*Col/2))) +mean;

%S(1:200,Location1(1:round(rate*Col/2))) + repmat(u(1:round(rate*Col/2)),200,1);
%S(201:400,Location1(round(rate*Col/2)+1:round(rate*Col)))= ...
%S(201:400,Location1(round(rate*Col/2)+1:round(rate*Col))) - repmat(u(round(rate*Col/2)+1:length(Location1)),200,1);

%S(:,Location(1:round(rate*Col/2)))=S(:,Location(1:round(rate*Col/2))) + repmat(u(1:round(rate*Col/2)),Row,1);
%S(:,Location(round(rate*Col/2)+1:round(rate*Col)))= S(:,Location(round(rate*Col/2)+1:round(rate*Col))) - repmat(u(round(rate*Col/2)+1:length(Location)),Row,1);
 NonP_Stripe=Ori+S;