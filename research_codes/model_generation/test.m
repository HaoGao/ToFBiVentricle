
theta_new = linspace(-pi, pi-0.005, 100);

figure
hold on
% 
  a=SA_INTERP4(or_data,theta_new);
 plot(a(1,:),a(2,:))

%  b=SA_INTERP2(or_data,theta_new);
%  plot(b(1,:),b(2,:))

%SA_INTERP3(or_data);



function  inter_data=SA_INTERP(or_data,theta_new)

points=[];theta=[];rho=[];z=[];
points = fnplt(cscvn(or_data));
%[theta,rho,z] = cart2pol(points(1,:),points(2,:),points(3,:));
center=mean(points,2);
center(1)=max(points(1,:))-(max(points(1,:))-min(points(1,:)))/3;
center(2)=max(points(2,:))-(max(points(2,:))-min(points(2,:)))/3;

for i=1:size(points,2)
    vec=points(1:2,i)-center(1:2);
    theta(i)=atan2(vec(2),vec(1));
    rho(i)=norm(vec);
end

[theta_uni,m1,n1] = unique(theta,'stable');
rho_uni=rho(m1);


%%% more cycles

uniqueNumbers = unique(n1);
counts = histcounts(n1, uniqueNumbers);  % combine the last two values together
maxcyle=max(counts);

r_new = interp1(theta_uni, rho_uni, theta_new, 'linear', 'extrap');

for i=1:length(theta_new)
    x(i)=r_new(i)*cos(theta_new(i))+center(1);
    y(i)=r_new(i)*sin(theta_new(i))+center(2);
end

z_new(1:length(theta_new))=center(3);
inter_data=[x; y; z_new];

return
end

%%%%%%%%%%%%%%%%
%%%   Updated function 

function  inter_data=SA_INTERP2(or_data,theta_new)

%or_data(:,length(or_data)+1)=or_data(:,1);
points_in=[];theta=[];rho=[];z=[];
points_in = fnplt(cscvn(or_data));
x=points_in(1,:);
y=points_in(2,:);
center=mean(points_in,2);

distances=[];
distances = sqrt(diff(x).^2 + diff(y).^2);
curveLength = sum(distances);


numPoints = 100;
segmentLength = curveLength / numPoints;
points = zeros(numPoints, 2);

currentLength = 0;
segmentIndex = 1;

for i = 1:length(distances)
    currentLength = currentLength + distances(i);

    if currentLength >= segmentIndex * segmentLength

        fraction = (segmentIndex * segmentLength - (currentLength - distances(i))) / distances(i);
        points(segmentIndex, 1) = x(i) + fraction * (x(i+1) - x(i));
        points(segmentIndex, 2) = y(i) + fraction * (y(i+1) - y(i));

        segmentIndex = segmentIndex + 1;

        if segmentIndex > numPoints
            break;
        end
    end
end

 
z_new(1:length(points))=center(3);
inter_data=[points';z_new];
%inter_data=points;


return
end

function inter_data=SA_INTERP3(or_data)

% 随机生成封闭曲线环的点
 numPoints = 100; % 曲线环上的点数量
% theta = linspace(0, 2*pi, numPoints+1); % 在0到2*pi之间均匀分布的角度向量
% theta(end) = []; % 去除最后一个角度，使曲线环首尾相接
% radius = 1 + 0.1*rand(1, numPoints); % 随机生成曲线环半径，范围为1到1.1
% x = radius .* cos(theta); % 曲线环的x坐标
% y = radius .* sin(theta); % 曲线环的y坐标

% % 示例封闭曲线环的坐标
% x = [0.5*cos(linspace(0, 2*pi, 200)) 0.5*cos(linspace(0, 2*pi, 200))];
% y = [0.5*sin(linspace(0, 2*pi, 200)) 1.5*sin(linspace(0, 2*pi, 200))];

% 计算曲线总长度
curveLength = sum(sqrt(diff(x).^2 + diff(y).^2));

% 计算每份的长度
segmentLength = curveLength / 100;

% 进行等长度插值
numPoints = numel(x);
points = zeros(100, 2); % 存储插值点的数组，每个点包含x和y坐标

totalLength = 0;
currentIndex = 1;

for i = 1:100
    targetLength = i * segmentLength;
    
    % 在曲线上找到最接近目标长度的位置
    while totalLength < targetLength
        segmentLength = sqrt((x(currentIndex+1)-x(currentIndex))^2 + (y(currentIndex+1)-y(currentIndex))^2);
        totalLength = totalLength + segmentLength;
        currentIndex = currentIndex + 1;
    end
    
    % 计算插值点的坐标
    alpha = (targetLength - totalLength + segmentLength) / segmentLength;
    points(i, 1) = (1-alpha) * x(currentIndex-1) + alpha * x(currentIndex);
    points(i, 2) = (1-alpha) * y(currentIndex-1) + alpha * y(currentIndex);
end

% 绘制等长度分割结果
figure;
plot(points(:, 1), points(:, 2), 'bo-', 'LineWidth', 1.5);
axis equal;

end




%%%   Updated function 

function  inter_data=SA_INTERP4(or_data,theta_new)

lot=length(or_data);
or_data(1,lot+1)=0.01*or_data(1,lot)+0.99*or_data(1,1);
or_data(2,lot+1)=0.01*or_data(2,lot)+0.99*or_data(2,1);
or_data(3,lot+1)=or_data(3,1);
%or_data(:,length(or_data)+1)=or_data(:,1);

points_in=[];theta=[];rho=[];z=[];
points_in = fnplt(cscvn(or_data));
x=points_in(1,:);
y=points_in(2,:);
center=mean(points_in,2);

distances=[];
distances = sqrt(diff(x).^2 + diff(y).^2);
curveLength = sum(distances);


numPoints = 1000;
segmentLength = curveLength / numPoints;
points = zeros(numPoints, 2);

currentLength1 = 0;
segmentIndex = 1;
cu=1;
points(1, 1) = x(cu);
points(1, 2) = y(cu);

for i = 1:numPoints-1
    
    aLength = i*segmentLength;


    currentLength2=currentLength1+distances(cu);


    if currentLength1< aLength & aLength<currentLength2
        if distances(cu) > 0.0001;
            fraction=(aLength-currentLength1)/distances(cu);
            points(i+1, 1) = x(cu) + fraction * (x(cu+1) - x(cu));
            points(i+1, 2) = y(cu) + fraction * (y(cu+1) - y(cu));
        else
            cu=cu+1;
                        fraction=(aLength-currentLength1)/distances(cu);
            points(i+1, 1) = x(cu) + fraction * (x(cu+1) - x(cu));
            points(i+1, 2) = y(cu) + fraction * (y(cu+1) - y(cu));
        end
    else
        currentLength1=currentLength2;
        cu=cu+1;
        if distances(cu) > 0.0001;
            fraction=(aLength-currentLength1)/distances(cu);
            points(i+1, 1) = x(cu) + fraction * (x(cu+1) - x(cu));
            points(i+1, 2) = y(cu) + fraction * (y(cu+1) - y(cu));
        else
            cu=cu+1;
                        fraction=(aLength-currentLength1)/distances(cu);
            points(i+1, 1) = x(cu) + fraction * (x(cu+1) - x(cu));
            points(i+1, 2) = y(cu) + fraction * (y(cu+1) - y(cu));
        end
    end

   
end

%%%%%% find 
f2=points(:,2)-mean(points(:,2));
f1=unique(abs(f2));
v1=1000;
k1=0;
while v1>center(1)
    k1=k1+1;
    p1=find(abs(f2)==f1(k1));
    v1=points(p1,1);
end

tarnum=100;
points_new(1,1)=points(p1,1);
points_new(1,2)=points(p1,2);

if points(p1+numPoints/tarnum,2)<points(p1,2)
    for i=2:tarnum
        p2=p1+numPoints/tarnum*(i-1);
        if p2>numPoints
            p2=p2-numPoints;
        end
        points_new(i,1)=points(p2,1);
        points_new(i,2)=points(p2,2);
    end
else
    for i=2:tarnum
        p2=p1-numPoints/tarnum*(i-1);
        if p2<=0
            p2=p2+numPoints;
        end
        points_new(i,1)=points(p2,1);
        points_new(i,2)=points(p2,2);
    end
end


z_new(1:length(points_new))=center(3);
inter_data=[points_new';z_new];



return
end
    
