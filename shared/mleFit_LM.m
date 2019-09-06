function [P,CRLB,LogL]=mleFit_LM(varargin)
% varargin:
%    imstack, startpsf/coeff, iterations, fitmode, isemccd, hidereport,initZ
%        1. imagestack (single)
%        2. fitmode
%           1 fix PSF
%           2 free PSF
%           3 Gauss fit z
%           4 fit PSFx, PSFy elliptical
%           5 cspline
%
%    optional:
%       3. iterations (default=30)
%       4. paramaters for fitters:
%          1. fix PSF: PSFxy sigma 
%          2. free PSF: start PSFxy
%          3. Gauss fit z: parameters: PSFx, Ax, Ay, Bx, By, gamma, d, PSFy
%               (single)
%          4. fit PSFx, PSFy elliptical: start PSFx, PSFy
%          5. cspline: cspline coefficients (single)

%       5. varmap: Variance map for sCMOS. If size(varmap) ~= size(imagestack):
%          no sCMOS correction is used. Default= emCCD
%       6. silent (suppress output if 1)
%       7. z start parameter (more than one: return solution with maximum
%          LIkelihood). Units: distance of stack calibration, center based

% Output:
% P
%   1. X, Y, Photons, Background, Iterations
%   2. X, Y, Photons, Background, PSFxy, Iterations
%   3. X, Y, Photons, Background, Z, Iterations
%   4. X, Y, Photons, Background, PSFx, PSFy, Iterations
%   5. X, Y, Photons, Background, Z, Iterations
%   6. X, Y, Photons, Background, Z, Iterations
% CRLB: cramer-rao lower bounds, as in P
% LogL: log-likelihood.

% Only for fitmode 6: P1 etc: results with z-startparameter<0, P2 etc:
% results with z-startparameter>0

P=[];CRLB=[];LogL=[]; % in case the function exits early
narginh=nargin;       %%% nargin是用来判断输入变量个数的函数。用法，a=nargin或者a=nargin(fx)。
% determine of it runs on GPU, otherwise use CPU as default
persistent fitter                                      %%% persistent 用于定义持久性变量。
allfitters={@GPUmleFit_LM,@CPUmleFit_LM};              %%% “@”符号 函数句柄，可以理解为一个函数的代号，这样在调用时可以调用函数句柄而不用调用该函数。比如sin是matlab中的一个函数，但sin只是函数名
                                                       %%% f = @sin;这行代码定义了一个函数句柄，变量名是f。这样就可以当做参数传递了，而且还可以跟sin函数按相同的语法规则使用：a = f(pi); %可以得到a=0
allfittersnames={'GPUmleFit_LM','CPUmleFit_LM'};       %%% 这句的意思是什么呢？
if isempty(fitter)        %%% isempty 判断数组的元素是否为空。
    testim=single(varargin{1}(:,:,1));                       %%% single函数把一个矩阵中所有元素都变为单精度的      varargin是一个元胞数组，用于代替所有变量。调用该函数时根据需要来改变输入参数的个数
    for k=1:length(allfitters)
        try                                                  %%% 程序首先运行try和catch之间的程序代码，如果没有发生错误则不执行catch和end之间的程序代码，而是执行end后的程序；
            allfitters{k}(testim,1);                         %%% 如果在执行try和catch之间的程序代码时产生错误，则立即执行catch和end之间的程序代码，然后继续执行end后的程序     
            fitter=k;
            break
        catch err
            % fitter did not work
        end
    end
    % disp(['using: ',char(allfitters(fitter))]);程序报错：元胞元素必须为字符数组。所以加%不运行     
end

%convert all parameters to the correct format
imagestack=single(varargin{1});

if narginh>1 && ~isempty(varargin{2})
    fitmode=varargin{2};
else
    fitmode=2;
end

if narginh>2 && ~isempty(varargin{3})
    iterations=varargin{3};
else
    iterations=30;
end

if narginh>3 && ~isempty(varargin{4})
    fitpar=single(varargin{4});
else
    if ~(fitmode==3) &&  ~(fitmode==5) 
        fitpar=1; %for Gaussian fit
    else
        disp('for z fitting (Gauss or spline) fitting parameters (e.g. spline coefficients) are required');
        return
    end
end

if narginh>4 && ~isempty(varargin{5})
    varmap=single(varargin{5});
else
    varmap=0;
end

if narginh>5 && ~isempty(varargin{6})
    silent=single(varargin{6});
else
    silent=1;
end

coeffsize=size(fitpar);
%backward compatibility
if fitmode==6
    fitmode=5;
    varargin{7}=[-coeffsize(3)/6, coeffsize(3)/6];
    narginh=max(narginh,7);
end

if narginh>6 && ~isempty(varargin{7}) && fitmode==5 %only for spline fitting
    z0=(varargin{7}); %in units of calibration stack distance
    zstart=single(z0+coeffsize(3)/2);
elseif fitmode==5
    zstart=single(coeffsize(3)/2);
else
    zstart=0;
end

[P,CRLB,LogL]=allfitters{fitter}(imagestack,fitmode,iterations,fitpar,varmap,silent,zstart(1));          %%% 因为不知道allfitters的内容，所以很难理解这条语句的内容？ 

if length(zstart)>1
    for k=2:length(zstart)
        [Ph,CRLBh,LogLh]=allfitters{fitter}(imagestack,fitmode,iterations,fitpar,varmap,silent,zstart(k));
   %      indbettero=LogLh<LogL;
        indbetter=LogLh-LogL>1e-4; %copy only everything if LogLh increases by more than rounding error.
        P(indbetter,:)=Ph(indbetter,:);
        CRLB(indbetter,:)=CRLBh(indbetter,:);
        LogL(indbetter)=LogLh(indbetter);
    end
end
 if varargin{2}==6 %2D fit: find proper results
     [P1,CRLB1,LogL1,P2,CRLB2,LogL2]=allfitters{fitter}(varargin{:});
     ind1=LogL1>=LogL2;
     ind2=LogL1<LogL2;
     P=zeros(size(P1),'single');CRLB=zeros(size(CRLB1),'single');LogL=zeros(size(LogL1),'single');
     P(ind1,:)=P1(ind1,:);P(ind2,:)=P2(ind2,:);
     CRLB(ind1,:)=CRLB1(ind1,:);CRLB(ind2,:)=CRLB2(ind2,:);
    LogL(ind1,:)=LogL1(ind1,:);LogL(ind2,:)=LogL2(ind2,:);
 else
     [P,CRLB,LogL]=allfitters{fitter}(varargin{:});
     P1=[];CRLB1=[];LogL1=[];P2=[];CRLB2=[];LogL2=[];
 end
%%
% 
% <latex>
% \begin{tabular}{|c|c|} \hline
% $n$ & $n!$ \\ \hline
% 1 & 1 \\
% 2 & 2 \\
% 3 & 6 \\ \hline
% \end{tabular}
% </latex>
% 
 clear(allfittersnames{fitter})

