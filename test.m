 clear all; close all
path(path,genpath(pwd));

load('MODIS.mat');load('feizhouqi250.210.mat')

J=(I+S)/255;

 psnr_com = 0;
 psnr=[];
for i=0.0055%[0.0022 0.003 0.0035 0.004 0.0045 0.005 0.0055]
     for j=0.0008%[0.0006 0.0007 0.0008 0.0009 0.00095]
         for v=0.001
             for k=0.00025%[0.00015 0.0002 0.00025 0.0003 0.00035]

opts.lamda1=i;   opts.lamda2=j;   opts.lamda3=k;
opts.beta1=v;   opts.beta2=v;   opts.beta3=v;

opts.tol=1.e-4;   opts.maxitr=2500;

theta=21; direction='r';
%[x] = Shear(J, theta, direction);
[u,s,ii,relchg]=adm_groupsparse(J,opts,theta,direction);
psnr1=psnr_fun(I/255,u);
ssim_g = ssim_index(I,u*255);
psnr=[psnr;i,j,k,v,psnr1,ssim_g];
             end
         end
     end
end
figure,subplot(221),imshow(I,[]);subplot(222),imshow(J,[]);
subplot(223),imshow(u,[]);subplot(224),imshow(s,[]);

