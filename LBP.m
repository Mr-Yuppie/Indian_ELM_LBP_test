%LBP returns the local binary pattern image or LBP histogram of an image.
% 返回图像的局部二值图或LBP直方图

% 1
%  J = LBP(I,R,N,MAPPING,MODE) returns either a local binary pattern coded image or the local binary pattern histogram of an intensity image I. 
% The LBP codes are computed using N sampling points on a circle of radius R and using mapping table defined by MAPPING. 
% 该LBP在半径为R的圆使用N个采样点并使用由MAPPING定义的映射表来计算。
% See the getmapping function for different mappings and use 0 for no mapping. Possible values for MODE are
%  请参阅getmapping函数的不同映射，并使用0表示无映射。 MODE的可能值为
% 'h' or 'hist'  to get a histogram of LBP codes LBP直方图
%       'nh'           to get a normalized histogram 归一化的直方图
%  Otherwise an LBP code image is returned. 或者返回LBP图像。

%  2
%  J = LBP(I) returns the original (basic) LBP histogram of image I  返回图像I的原始（基本）LBP直方图

%  3
%  J = LBP(I,SP,MAPPING,MODE) computes the LBP codes using n sampling points defined in (n * 2) matrix SP. The sampling points should be defined around the origin (coordinates (0,0)).
%  使用在（n * 2）矩阵SP中定义的n个采样点来计算LBP码。 采样点应围绕原点（坐标（0,0））进行定义。

%  Examples
%  --------
%       I=imread('rice.png');
%       mapping=getmapping(8,'u2'); %using uniform patterns
%       H1=LBP(I,1,8,mapping,'h'); %LBP histogram in (8,1) neighborhood
%                                  
%       subplot(2,1,1),stem(H1);
%
%       H2=LBP(I);
%       subplot(2,1,2),stem(H2);
%
%       SP=[-1 -1; -1 0; -1 1; 0 -1; -0 1; 1 -1; 1 0; 1 1];
%       I2=LBP(I,SP,0,'i'); %LBP code image using sampling points in SP
%                           %and no mapping. Now H2 is equal to histogram
%                           %of I2.

function [result,lbp_map] = lbp(varargin) % image,radius,neighbors,mapping,mode)
% Version 0.3.2
% Authors: Marko Heikkil?and Timo Ahonen

% Changelog
% Version 0.3.2: A bug fix to enable using mappings together with a predefined spoints array
% Version 0.3.1: Changed MAPPING input to be a struct containing the mapping
% table and the number of bins to make the function run faster with high number
% of sampling points. Lauge Sorensen is acknowledged for spotting this problem.
% 将MAPPING输入更改为包含映射表和数据块的结构，bins数目使函数在具有具有大量采样点时运行速度更快，

% Check number of input arguments.
error(nargchk(1,5,nargin));

image=varargin{1};

if nargin==1
    spoints=[-1 -1; -1 0; -1 1; 0 -1; -0 1; 1 -1; 1 0; 1 1];
    neighbors=8;
    mapping=0;
    mode='h';
end

if (nargin == 2) && (length(varargin{2}) == 1)
    error('Input arguments');
end

if (nargin > 2) && (length(varargin{2}) == 1)
    radius=varargin{2};
    neighbors=varargin{3};
    
    spoints=zeros(neighbors,2);

    % Angle step.
    a = 2*pi/neighbors;
    
    for i = 1:neighbors
        spoints(i,1) = -radius*sin((i-1)*a);
        spoints(i,2) = radius*cos((i-1)*a);
    end
    
    if(nargin >= 4)
        mapping=varargin{4};
        if(isstruct(mapping) && mapping.samples ~= neighbors)
            error('Incompatible mapping');
        end
    else
        mapping=0;
    end
    
    if(nargin >= 5)
        mode=varargin{5};
    else
        mode='h';
    end
end

if (nargin > 1) && (length(varargin{2}) > 1)
    spoints=varargin{2};
    neighbors=size(spoints,1);
    
    if(nargin >= 3)
        mapping=varargin{3};
        if(isstruct(mapping) && mapping.samples ~= neighbors)
            error('Incompatible mapping');
        end
    else
        mapping=0;
    end
    
    if(nargin >= 4)
        mode=varargin{4};
    else
        mode='h';
    end   
end

% Determine the dimensions of the input image. 输入图像的大小
[m,n,d] = size(image);

% padding elements avoid edges 填充元素避免边缘
W = radius;
X = zeros(m+2*W, n+2*W, d);
X(W+1:m+W, W+1:n+W, :) = image;
X(W+1:m+W, 1:W, :) = image(:, W:-1:1, :);
X(W+1:m+W, n+W+1:n+2*W, :) = image(:, n:-1:n-(W-1), :);
X(1:W, :, :) = X(2*W:-1:(W+1), :, :);
X(m+(W+1):m+2*W, :, :) = X(m+W:-1:(m+1), :, :);


image=X;
d_image=double(image);
[ysize xsize] = size(image);



miny=min(spoints(:,1));
maxy=max(spoints(:,1));
minx=min(spoints(:,2));
maxx=max(spoints(:,2));
%块大小，在块内计算LBP
% Block size, each LBP code is computed within a block of size bsizey*bsizex
bsizey=ceil(max(maxy,0))-floor(min(miny,0))+1;
bsizex=ceil(max(maxx,0))-floor(min(minx,0))+1;

% Coordinates of origin (0,0) in the block （0，0）坐标在块中
origy=1-floor(min(miny,0));
origx=1-floor(min(minx,0));

% Minimum allowed size for the input image depends on the radius of the used LBP operator.
% 输入图像的最小允许尺寸取决于所使用的LBP操作符的半径。
if(xsize < bsizex || ysize < bsizey)
  error('LBP.m--Too small input image. Should be at least (2*radius+1) x (2*radius+1)');
end

% Calculate dx and dy;
dx = xsize - bsizex;
dy = ysize - bsizey;

% Fill the center pixel matrix C.
C = image(origy:origy+dy,origx:origx+dx);
d_C = double(C);

bins = 2^neighbors;

% Initialize the result matrix with zeros.
result=zeros(dy+1,dx+1);

%Compute the LBP code image

for i = 1:neighbors
  y = spoints(i,1)+origy;
  x = spoints(i,2)+origx;
  % Calculate floors, ceils and rounds for the x and y.
  fy = floor(y); cy = ceil(y); ry = round(y);
  fx = floor(x); cx = ceil(x); rx = round(x);
  % Check if interpolation is needed.
  if (abs(x - rx) < 1e-6) && (abs(y - ry) < 1e-6)
    % Interpolation is not needed, use original datatypes
    N = image(ry:ry+dy,rx:rx+dx);
    D = N >= C; 
  else
    % Interpolation needed, use double type images 
    ty = y - fy;
    tx = x - fx;

    % Calculate the interpolation weights.
    w1 = (1 - tx) * (1 - ty);
    w2 =      tx  * (1 - ty);
    w3 = (1 - tx) *      ty ;
    w4 =      tx  *      ty ;
    % Compute interpolated pixel values
    N = w1*d_image(fy:fy+dy,fx:fx+dx) + w2*d_image(fy:fy+dy,cx:cx+dx) + ...
        w3*d_image(cy:cy+dy,fx:fx+dx) + w4*d_image(cy:cy+dy,cx:cx+dx);
    D = N >= d_C; 
  end  
  % Update the result matrix.
  v = 2^(i-1);
  result = result + v*D;
end
lbp_map=result;

%Apply mapping if it is defined
if isstruct(mapping)
    bins = mapping.num;
    for i = 1:size(result,1)
        for j = 1:size(result,2)
            result(i,j) = mapping.table(result(i,j)+1);
        end
    end
end

if (strcmp(mode,'h') || strcmp(mode,'hist') || strcmp(mode,'nh'))
    % Return with LBP histogram if mode equals 'hist'.
    result=hist(result(:),0:(bins-1));
    if (strcmp(mode,'nh'))
        result=result/sum(result);
    end
else
    %Otherwise return a matrix of unsigned integers
    if ((bins-1)<=intmax('uint8'))
        result=uint8(result);
    elseif ((bins-1)<=intmax('uint16'))
        result=uint16(result);
    else
        result=uint32(result);
    end
end

end




