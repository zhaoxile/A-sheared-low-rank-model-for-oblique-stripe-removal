function [u,s,ii,relchg]=adm_groupsparse(I,opts,theta,direction)
%% ---the right version!!!!
%---date: 2016-5-10
[m,n]=size(I);
Dt= defDDt;
C=getC(I);
D1=defDD1t;  %y-direction finite function
D2=defDD2t;  %x-direction finite function
%%%%--------------initialization----------%%%%
f=I;
s=zeros(m,n);  
u=zeros(m,n);
w=zeros(m,n);
p1=zeros(m,n);  %%!!!
p2=p1;
p3=p1;


up = u;
lamda1=opts.lamda1;
lamda2=opts.lamda2;
lamda3=opts.lamda3;%lamda1=0;

beta1=opts.beta1;
beta2=opts.beta2;
beta3=opts.beta3;%beta1=0;


tol=opts.tol;
maxitr=opts.maxitr;
%Denom=beta2*ones(m,n)+beta3*C.eigsD2tD2;%% eyes(m,n)
[s1] = Shear(s, theta, direction);
D1u=D1(u);


%%%%--------------finite diff-------------%%%%
%%%%--------------Main loop---------------%%%%
ii=0;
relchg1=1;

relchg=[];
while relchg1> tol && ii<maxitr 
    V1=D1u+p3/beta3;

    %%%%----------d-subproblem------------%%%%
    v=sign(V1).*max(0,abs(V1)-lamda3/beta3);

    %%%%----------group sparse l_{2,1}----%%%%
    
    for i=1:n
        r=s1(:,i)+p2(:,i)/beta2;
        q(:,i)=r.*max(norm(r)-lamda2/beta2,0)/(norm(r)+eps);
    end
   %q=zeros(m,n); 
     %%%-------global sparse l_0--------%%%%%
%      V3=s+p2/beta2;
%      jj=find(abs(V3)>=sqrt(2*lamda1/beta2));
%      y(jj)=V3(jj);
 %%%%----------d-subproblem------------%%%%
      [U,a,V] = svd (s1+p1/beta1);
    [c,e]=size(a);
    if c>e
        c=e;
    end
    for i=1:c
        a(i,i)=a(i,i)-lamda1/beta1;
        if a(i,i)<0
            a(i,i)=0;
        end
    end
    d=U*a*V';
 %d=zeros(m,n);
    %%%%----------w-subproblem------------%%%%
    [q3] = ShearTrans(q-p2/beta2, theta, direction);% q3=zeros(m,n);
    [d1] = ShearTrans(d-p1/beta1, theta, direction);%d1=zeros(m,n);
    
    a21=ones(m,n);
    a12=ones(m,n);
    a11=ones(m,n)+beta3*(abs(psf2otf([1,-1],[m n])).^2);
    a22=ones(m,n)+beta2*ones(m,n)+beta1*ones(m,n);%%

    b1 = fft2(f+Dt((v - p3/beta3)*beta3,ones(m,n)));

    b2 = fft2(f+beta1*d1+beta2*q3);
    %%%solve the linear system ax=b by Cramer's Rule
    Det=a11.*a22-a12.*a21;
   % Det(Det==0) = 1;
    u1=(b1.*a22-b2.*a21)./(Det);
    s1=(a11.*b2-a12.*b1)./(Det);
    s=real(ifft2(s1));
    u=real(ifft2(u1));
    %  Update Lam
    % =================  
    %%%%%%%%%%%%%%%%%%%%%%%
    relchg1=norm(u-up,'fro')/norm(u,'fro');
    relchg = [relchg,relchg1];
    ii=ii+1;
    %%%%------------update p--------------%%%%
%    [s1] = Shear(s, theta, direction);
     D1u=D1(u);
     [s1]= Shear(s, theta, direction);
    p1=p1+beta1*(s1-d);
    p2=p2+beta2*(s1-q);
    
    p3=p3+beta3*(D1u-v);

    up = u;
end
%%%%--------------Subfunction-------------%%%%
function C=getC(I)
sizeI=size(I);
C.eigsD1tD1=abs(psf2otf([1,-1],sizeI)).^2;
C.eigsD2tD2=abs(psf2otf([1;-1],sizeI)).^2;
end
%%%%--------------Subfunction-------------%%%%
function D1=defDD1t
D1=@(U)ForwardD1(U);
end

function D2=defDD2t
D2=@(U)ForwardD2(U);
end

function Dt=defDDt
Dt= @(X,Y) Dive(X,Y);
end

function Dux=ForwardD1(U)
Dux=[diff(U,1,2),U(:,1)-U(:,end)];
end

function Duy=ForwardD2(U)
Duy=[diff(U,1,1);U(1,:)-U(end,:)];
end

function DtXY = Dive(X,Y)
DtXY = [X(:,end) - X(:, 1), -diff(X,1,2)];
DtXY = DtXY + [Y(end,:) - Y(1, :); -diff(Y,1,1)];
end
end