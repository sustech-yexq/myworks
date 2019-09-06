%  Copyright (c)2017 Ries Lab, European Molecular Biology Laboratory,
%  Heidelberg.
%  
%  This file is part of GPUmleFit_LM Fitter.
%  
%  GPUmleFit_LM Fitter is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%  
%  GPUmleFit_LM Fitter is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%  
%  You should have received a copy of the GNU General Public License
%  along with GPUmleFit_LM Fitter.  If not, see <http://www.gnu.org/licenses/>.
%  
%  
%  Additional permission under GNU GPL version 3 section 7

%% example script for fitting of 3D data
%  Installation instructions
%  Requires Matlab 2016a or newer including the following tool boxes:      
%  image analysis, curve fitting, optimization, statistics                                                                                                   
%  GPU support for CUDA 8               不知道这个要求是否达到                                                                                                                                                                                         
clc;


%% add path to helper functions
addpath('shared')        %%% 添加“shared”文件夹作为路径                                                

%% make bead calibration
%  run 3D calibration GUI and make 3D calibration
%  e.g. using bead files: in example_data/beadstacks_3D_astig                                                                                                                                                                                                           
%  save e.g. as example_data/bead_astig_3dcal.mat
% or generate calibration files programmatically:

if ~exist([pwd filesep 'example_data' filesep 'bead_astig_3dcal.mat'],'file')   %   only run if no calibration files have been generated, as this takes time.
    example_calibration     %%%   如果文件夹example_data中不存在bead_astig_3dcal.mat,则运行这个脚本。这个脚本的目的是产生这个bead_astig_3dcal.mat数据。
                            %%%   但这个函数运行不了，报错，因为其中一个调用函数缺失。
end

%% load bead calibration
cal=load([pwd filesep 'example_data' filesep 'bead_astig_3dcal.mat']); % load bead calibration    % 加载这些数据给cal。 这些都是校准的数据。


%% either simulate data or load experimental tiff file      或者模拟数据或者加载实验tiff文件，以下程序分这两种情况进行
mode=2;
if mode ==1   % simulate data    模拟数据
    offset=0;       %%% 这个变量干嘛用的     
    conversion=1;   %%% 这个变量干嘛用的
    %numlocs=50000;    %  numlocs: number of simulated molecules. For maximum fitting speeds on the GPU,
                      %  it should be >10^5. For smaller GPUs, adust to avoid memory overflow
                      %  error to 10^5-10^5. For CPU fitter, use 10^3-10^4.
     
    RoiPixelsize=13;  % ROI in pixels for simulation    %%%  ROI（region of interest）是从被处理的图像以方框、圆、椭圆等等方式勾勒出需要处理的区域，称为感兴趣区。
    dz=cal.cspline.dz;  % coordinate system of spline PSF is corner based and in units pixels / planes 完全不理解这句解释的意思
    z0=cal.cspline.z0; % distance and midpoint of stack in spline PSF, needed to translate into nm coordinates  z0是指什么
    dx=floor(RoiPixelsize/2); % distance of center from edge
    ground_truth.z=linspace(-500,500,50000)'; % define some coordinates. Alternatively, use rand    %%% linspace（x1,x2,N)用于产生x1，x2之间的N个值的行矢量，相数数值间距相等。 
    ground_truth.x=linspace(-0.5,0.5,50000)';
    ground_truth.y=sin(ground_truth.x*4*pi);      %%% 为什么y的值是x的sin函数。
    coordinates=horzcat(ground_truth.x+dx,ground_truth.y+dx,ground_truth.z/dz+z0);
    Intensity=2000; % number of photons / localization
    background=10;  % number of background photons per pixel
    imstack = simSplinePSF(RoiPixelsize,cal.cspline.coeff,Intensity,background,coordinates); % simulate images   调用函数，赋值给imstack。目前还不理解imstack的数据
               %%% 因为不理解imstack这个数据，所以不这个函数的作用是什么
               %%% coeff是什么，因为打不开它的数据，所以不理解它是什么。
               
    imstack=single(imstack); % imstack needs to be in photons. The fitters require the stacks in single-format;
    
else % experimental data  实验数据
    file=[pwd filesep 'example_data' filesep 'single_bead.tif']; % experimental astgmatic bead stack.    这个astgmatic怎么翻译
    imstackadu=readfile_tif(file); % Stack of ROIs in photons.   
    offset=400;        %%% 这个变量的意思是什么呢？
    conversion=0.1;     %%% 同样，这个变量的意思是什么呢？
    imstack=(single(imstackadu)-offset)/conversion;  % if necessary, convert ADU into photons.  
    ground_truth.z=((1:size(imstack,3))'-size(imstack,3)/2+1)*10; %  dz=10 nm; convert frame to nm   %%% 为什么真实值是这样表示？
end


%% fit image stacks

% load sCMOS varmap or set sCMOSvarmap to zero
% sCMOSvarmap=ones(size(imstack),'single'); %the variance map, same size as imstack, single format. 
sCMOSvarmap=0; % if scalar : use EMCCD fitter;    为什么这个要设为0
numlocs=size(imstack,3); 


if isfield(cal,'gauss_zfit')&&~isempty(cal.gauss_zfit) % only if Gaussian PSF model was calibrated in the calibrate3D_GUI, we can test the Gaussian fit
    % fit Gaussian model, direct z fit, emCCD mode                             %%%  isfield(cal,'gauss_zfit')语句是判断 gauss_zfit是否为cal的成员。
    
    tic      %%%%记录此刻的时间，即开始时间，与toc搭配使用。
    [P_Gau,CRLB3,tianjia1]=mleFit_LM(imstack,3,50,single(cal.gauss_zfit),sCMOSvarmap,1);    %%% 添加了一个输出变量tianjia1，为了对应mleFit_LM函数的输出个数
                      %%% 这个函数是最大似然估计函数吗？返回的参数P是什么？
    tgz=toc; %%% 保存运行完的时间，并赋值给tgz
    disp(['Gauss z: ' num2str(numlocs/tgz) ' fits/s']);
    gaussz.x=P_Gau(:,1);
    gaussz.y=P_Gau(:,2); %x,y in pixels       %%%   P有7列，那3，4，6，7分别是什么？为什么1和2分别给x和y，而不是其他列？
    gaussz.z=P_Gau(:,5)*1000;            %%% 为什么要乘1000呢？

    % fit elliptical Gaussian model, extract z from sx, sy. emCCD mode
    tic
    [P_ell_Gau,CRLB4,tianjia2]=mleFit_LM(imstack,4,50,1,sCMOSvarmap,1);                     %%% 添加了一个变量tianjia2，为了对应mleFit_LM函数的输出个数
    tgsxsy=toc;
    disp(['Gaussxy: ' num2str(numlocs/tgsxsy) ' fits/s']);
    gausssxsy.x=P_ell_Gau(:,1);
    gausssxsy.y=P_ell_Gau(:,2); %x,y in pixels 
    sx=P_ell_Gau(:,5);
    sy=P_ell_Gau(:,6);         %%% sx和sy分别是什么呢？
    gausssxsy.z=sxsy2z(sx,sy,cal.gauss_sx2_sy2); % z from sx, sy      %%%  调用函数sxsy2z，获取z的值

else
    gaussz.x=zeros(size(imstack,3),1);gaussz.y=zeros(size(imstack,3),1); gaussz.z=zeros(size(imstack,3),1);
    gausssxsy.x=zeros(size(imstack,3),1);gausssxsy.y=zeros(size(imstack,3),1); gausssxsy.z=zeros(size(imstack,3),1);
end

% fit elliptical Gaussian model, cspline, emCCD mode
tic
[P_cspline,CRLB,tianjia3]=mleFit_LM(imstack,5,50,single(cal.cspline.coeff),sCMOSvarmap,1);        %%% 添加了一个变量tianjia3，为了对应mleFit_LM函数的输出个数
tspline=toc;
disp(['cspline: ' num2str(numlocs/tspline) 'fits/s']);
dx=floor(size(imstack,1)/2);   %%% 所以dx是什么意思呢？
cspline.x=P_cspline(:,1)-dx;   %%% 为什么要减去这个dx呢？
cspline.y=P_cspline(:,2)-dx; %x,y in pixels 
cspline.z=(P_cspline(:,5)-cal.cspline.z0)*cal.cspline.dz;

z_err=sqrt(CRLB(:,5))*cal.cspline.dz; %CRLB is the variance of the parameters, to obtain a localization precision we have to take the square root

% calculate error for all fits. 
zrange=400; % Only take into account central part. Focus +/- zrange
inz=abs(ground_truth.z)<zrange;

% difference between fitted z and ground truth
gaussz.zrel=gaussz.z-ground_truth.z;
gaussz.zrel(isnan(gaussz.zrel))=0;    %%% isnan判断数组中的元素是否为 NaN。该语句将gaussz.zrel数组中为NaN的元素设为0。

gausssxsy.zrel=gausssxsy.z-ground_truth.z;
gausssxsy.zrel(isnan(gausssxsy.zrel))=0;

cspline.zrel=cspline.z-ground_truth.z;

% numpoints
numpoints=min(2000,size(imstack,3));
range=1:round(length(ground_truth.z)/numpoints):length(ground_truth.z);
%plot fitted z vs ground-truth z
    figure(1);
    % hold off
    plot(ground_truth.z(range),gaussz.zrel(range)-mean(gaussz.zrel(inz)),'r.');
    hold on
    plot(ground_truth.z(range),gausssxsy.zrel(range)-mean(gausssxsy.zrel(inz)),'.','Color',[0.7,0.7,0]);     %%% [0.7,0.7,0] 这个是什么的范围？
    hold on
    plot(ground_truth.z(range),cspline.zrel(range)-mean(cspline.zrel(inz)),'b.');     
    hold on 
    
    gs=fit(ground_truth.z(range),double(gaussz.zrel(range)-mean(gaussz.zrel(inz))),'smoothingspline','SmoothingParam',.00001);
    gss=fit(ground_truth.z(range),double(gausssxsy.zrel(range)-mean(gausssxsy.zrel(inz))),'smoothingspline','SmoothingParam',.00001);
    gz=fit(ground_truth.z(range),double(cspline.zrel(range)-mean(cspline.zrel(inz))),'smoothingspline','SmoothingParam',.00001);

    plot(ground_truth.z(range),ground_truth.z(range)*0,'k','LineWidth',3); %%% 这里是要画z=0的位置吗？
    hold on

    plot(ground_truth.z(range),gs(ground_truth.z(range)),'r','LineWidth',3);    %%% 作出拟合曲线
    hold on
    plot(ground_truth.z(range),gss(ground_truth.z(range)),'Color',[0.7,0.7,0],'LineWidth',3);
    hold on
    plot(ground_truth.z(range),gz(ground_truth.z(range)),'b','LineWidth',3);
    hold on
    plot([-zrange -zrange],[-1000 1000],'k');
    hold on
    plot([zrange zrange],[-1000 1000],'k');
    
if mode==1
   tt=10;
else
   tt=1000;
end
    ss=std(cspline.zrel(inz))*tt;         %%% std(x)算出x的标准偏差
    ylim([-ss ss]);
    xlabel('ground truth z (nm)');
    ylabel('z (fitted) - z (ground truth) (nm)');
    title('Accuracy');

    cspline.dz=std(ground_truth.z(inz)-cspline.z(inz),'omitnan');   %%% std函数用来计算标准偏差, "omitnan" 从结果中忽略 NaN 值
    zgauss.dz=std(ground_truth.z(inz)-gaussz.z(inz),'omitnan');
    gausssxsy.dz=std(ground_truth.z(inz)-gausssxsy.z(inz),'omitnan');

    legendtxt{4}='smoothing spline';
    legendtxt{3}=['spline fit. error: ' num2str(cspline.dz, '%3.1f') ' nm'];
    legendtxt{1}=['Gaussian z fit. error: ' num2str(zgauss.dz, '%3.1f') ' nm'];
    legendtxt{2}=['Gaussian sx, sy fit. error: ' num2str(gausssxsy.dz, '%3.1f') ' nm'];
    legend(legendtxt)           %%% legend 函数，作图标，显示（）里的内容

 %calculate lateral error
 if mode == 1
     cspline.dx=std(ground_truth.x(inz)-cspline.x(inz),'omitnan'); %in pixels
     cspline.dy=std(ground_truth.y(inz)-cspline.y(inz),'omitnan');
     dx_gaussz=std(ground_truth.x(inz)-gaussz.x(inz),'omitnan');
     dy_gaussz=std(ground_truth.y(inz)-gaussz.y(inz),'omitnan');
     gausssxsy.dx=std(ground_truth.x(inz)-gausssxsy.x(inz),'omitnan');
     gausssxsy.dy=std(ground_truth.y(inz)-gausssxsy.y(inz),'omitnan');
     %plot 3D scatter plot
    
     figure(2);
     hold off
     scatter3(ground_truth.x(range),ground_truth.y(range),ground_truth.z(range),3);      %%% 绘制三维散点图，scatter3(x,y,z,'.',c) % c 为颜色，x,y,z长度需相同
     hold on
     scatter3(cspline.x(range),cspline.y(range),cspline.z(range),5);
     xlabel('x (nm)');ylabel('y (nm)');zlabel('z (nm)');
     legend('ground truth','cspline fit')
     title('fitted vs ground truth positions')
 else
     cspline.dx=std(cspline.x(inz),'omitnan'); %in pixels
     cspline.dy=std(cspline.y(inz),'omitnan');
     dx_gaussz=std(gaussz.x(inz),'omitnan');
     dy_gaussz=std(gaussz.y(inz),'omitnan');
     gausssxsy.dx=std(gausssxsy.x(inz),'omitnan');
     gausssxsy.dy=std(gausssxsy.y(inz),'omitnan');
 end

