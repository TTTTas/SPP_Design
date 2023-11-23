clear
close
data1=importfile('双频双系统.pos');
data2=importfile('单频双系统.pos');
data3=importfile('单频GPS系统.pos');
data4=importfile('单频BDS系统.pos');
data5=importfile('双频GPS系统.pos');
data6=importfile('双频BDS系统.pos');
data7=importfile('无对流层改正.pos');
Ref=[-2267805.437,5009342.139,3220992.182];
%% 
enu1=XYZ2ENU(Ref,data1(:,2:4),0);
enu2=XYZ2ENU(Ref,data2(:,2:4),0);
enu3=XYZ2ENU(Ref,data3(:,2:4),0);
enu4=XYZ2ENU(Ref,data4(:,2:4),0);
enu5=XYZ2ENU(Ref,data5(:,2:4),0);
enu6=XYZ2ENU(Ref,data6(:,2:4),0);
enu7=XYZ2ENU(Ref,data7(:,2:4),0);
%% 
figure(1)
subplot(1,3,1)
scatter(enu1(:,1),enu1(:,2),4,[0.635294117647059 0.0784313725490196 0.184313725490196],'filled','o','DisplayName','双频双系统');
hold on;
scatter(enu2(:,1),enu2(:,2),4,[0 0.447058823529412 0.741176470588235],'+','DisplayName','单频双系统');
scatter(0,0,40,'red','filled','square','DisplayName','Ref');
grid on
legend
subplot(1,3,2)
scatter(enu3(:,1),enu3(:,2),4,[1 0.411764705882353 0.16078431372549],'^','DisplayName','单频GPS系统');
hold on;
scatter(enu5(:,1),enu5(:,2),4,[0.552941176470588 0.756862745098039 0.988235294117647],'x','DisplayName','双频GPS系统');
scatter(0,0,40,'red','filled','square','DisplayName','Ref');
grid on
legend
subplot(1,3,3)
scatter(enu4(:,1),enu4(:,2),4,[0.717647058823529 0.274509803921569 1],'*','DisplayName','单频BDS系统');
hold on
scatter(enu6(:,1),enu6(:,2),4,[0.752941176470588 0.909803921568627 0.572549019607843],'<','DisplayName','双频BDS系统');
scatter(0,0,40,'red','filled','square','DisplayName','Ref');
grid on
legend
figure(2)
scatter(enu7(:,1),enu7(:,2),4,[0.96078431372549 0.717647058823529 0.231372549019608],'v','DisplayName','无对流层组合');
hold on
scatter(enu1(:,1),enu1(:,2),4,[0.635294117647059 0.0784313725490196 0.184313725490196],'filled','o','DisplayName','双频双系统');
scatter(0,0,40,'red','filled','square','DisplayName','Ref');
grid on
axis equal
legend
%% 
figure(3)
scatter(enu1(:,1),enu1(:,2),4,[0.635294117647059 0.0784313725490196 0.184313725490196],'filled','o','DisplayName','双频双系统');
hold on;
scatter(enu2(:,1),enu2(:,2),4,[0 0.447058823529412 0.741176470588235],'+','DisplayName','单频双系统');
scatter(enu3(:,1),enu3(:,2),4,[1 0.411764705882353 0.16078431372549],'^','DisplayName','单频GPS系统');
scatter(enu5(:,1),enu5(:,2),4,[0.552941176470588 0.756862745098039 0.988235294117647],'x','DisplayName','双频GPS系统');
scatter(enu4(:,1),enu4(:,2),4,[0.717647058823529 0.274509803921569 1],'*','DisplayName','单频BDS系统');
scatter(enu6(:,1),enu6(:,2),4,[0.752941176470588 0.909803921568627 0.572549019607843],'<','DisplayName','双频BDS系统');
scatter(enu7(:,1),enu7(:,2),4,[0.96078431372549 0.717647058823529 0.231372549019608],'v','DisplayName','无对流层组合');
scatter(0,0,40,'red','filled','square','DisplayName','Ref');
grid on
legend
axis equal
%% 
function blh = XYZ2BLH(xyz, e2, a)
if nargin < 3
    error('XYZ2BLH: 输入参数不足!');
end
if nargin > 3
    error('XYZ2BLH: 输入参数过多!');
end

% 计算经度
L = 180 * atan2(xyz(2), xyz(1)) / pi;
while (L < 0 || L > 360)
    if (L < 0)
        L = L+360;
    end
    if (L > 360)
        L = L-360;
    end
end
blh(1) = L;

% 计算纬度
B = atan(xyz(3) / sqrt(xyz(1) * xyz(1) + xyz(2) * xyz(2)));
B0 = 0;
count = 0;
while (abs(B - B0) > 1e-10 && count < 100)
    B0 = B;
    n = a / sqrt(1 - e2 * sin(B0) * sin(B0));
    B = atan((xyz(3) + n * e2 * sin(B0)) / sqrt(xyz(1) * xyz(1) + xyz(2) * xyz(2)));
    count=count+1;
end
n = a / sqrt(1 - e2 * sin(B) * sin(B));
blh(2) = 180 * B / pi;
blh(3) = xyz(3) / sin(B) - n * (1 - e2);

% 返回结果
end
%% 
function enu = XYZ2ENU(xyz1, xyz2, sys)
WGS84_e2=0.0066943799013;
WGS84_a=6378137.0;
CGCS2000_e2=0.00669438002290;
CGCS2000_a=6378137.0;
if nargin < 3
    error('XYZ2ENU: 输入参数不足!');
end
if nargin > 3
    error('XYZ2ENU: 输入参数过多!');
end

% 转换坐标系
switch sys
    case 0
        blh = XYZ2BLH(xyz1, WGS84_e2, WGS84_a);
    case 1
        blh = XYZ2BLH(xyz1, CGCS2000_e2, CGCS2000_a);
    otherwise
        error('XYZ2ENU: 不支持的坐标系!');
end
enu=zeros(height(xyz2),3);
% 计算ENU坐标
R = [-sin(blh(2)),cos(blh(2)),0; 
    -sin(blh(1)) * cos(blh(2)),-sin(blh(1)) * sin(blh(2)),cos(blh(1)); 
    cos(blh(1)) * cos(blh(2)),cos(blh(1)) * sin(blh(2)),cos(blh(1))];
for i=1:height(xyz2)
    x0 = xyz2(i,1:3)-xyz1;
    enu(i,1:3) = R * x0.';
end

% 返回结果
end