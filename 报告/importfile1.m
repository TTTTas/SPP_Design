function [Sow, X, Y, Z, B, L, H, E, N, U, Gclock, Cclock, m_H, m_V, Vx, Vy, Vz, Vc, sigma0, PDOP, GS, BS] = importfile1(filename, dataLines)
%IMPORTFILE2 从文本文件中导入数据
%  [SOW, X, Y, Z, B, L, H, E, N, U, GCLOCK, CCLOCK, M_H, M_V, VX, VY,
%  VZ, VC, SIGMA0, PDOP, GS, BS] = IMPORTFILE2(FILENAME)读取文本文件 FILENAME
%  中默认选定范围的数据。  以列向量形式返回数据。
%
%  [SOW, X, Y, Z, B, L, H, E, N, U, GCLOCK, CCLOCK, M_H, M_V, VX, VY,
%  VZ, VC, SIGMA0, PDOP, GS, BS] = IMPORTFILE2(FILE,
%  DATALINES)按指定行间隔读取文本文件 FILENAME 中的数据。对于不连续的行间隔，请将 DATALINES 指定为正整数标量或
%  N×2 正整数标量数组。
%
%  示例:
%  [Sow, X, Y, Z, B, L, H, E, N, U, Gclock, Cclock, m_H, m_V, Vx, Vy, Vz, Vc, sigma0, PDOP, GS, BS] = importfile2("D:\GitHub\SPP_Design\报告\单频双系统.pos", [1, Inf]);
%
%  另请参阅 READTABLE。
%
% 由 MATLAB 于 2023-11-20 10:49:11 自动生成

%% 输入处理

% 如果不指定 dataLines，请定义默认范围
if nargin < 2
    dataLines = [1, Inf];
end

%% 设置导入选项并导入数据
opts = delimitedTextImportOptions("NumVariables", 39);

% 指定范围和分隔符
opts.DataLines = dataLines;
opts.Delimiter = ["\t", ":"];

% 指定列名称和类型
opts.VariableNames = ["Var1", "Var2", "Sow", "Var4", "X", "Y", "Z", "Var8", "B", "L", "H", "Var12", "E", "N", "U", "Var16", "Gclock", "Var18", "Cclock", "Var20", "m_H", "Var22", "m_V", "Var24", "Vx", "Vy", "Vz", "Vc", "Var29", "sigma0", "Var31", "Var32", "Var33", "PDOP", "Var35", "GS", "Var37", "BS", "Var39"];
opts.SelectedVariableNames = ["Sow", "X", "Y", "Z", "B", "L", "H", "E", "N", "U", "Gclock", "Cclock", "m_H", "m_V", "Vx", "Vy", "Vz", "Vc", "sigma0", "PDOP", "GS", "BS"];
opts.VariableTypes = ["string", "string", "double", "string", "double", "double", "double", "string", "double", "double", "double", "string", "double", "double", "double", "string", "double", "string", "double", "string", "double", "string", "double", "string", "double", "double", "double", "double", "string", "double", "string", "string", "string", "double", "string", "double", "string", "double", "string"];

% 指定文件级属性
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% 指定变量属性
opts = setvaropts(opts, ["Var1", "Var2", "Var4", "Var8", "Var12", "Var16", "Var18", "Var20", "Var22", "Var24", "Var29", "Var31", "Var32", "Var33", "Var35", "Var37", "Var39"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var1", "Var2", "Var4", "Var8", "Var12", "Var16", "Var18", "Var20", "Var22", "Var24", "Var29", "Var31", "Var32", "Var33", "Var35", "Var37", "Var39"], "EmptyFieldRule", "auto");

% 导入数据
tbl = readtable(filename, opts);

%% 转换为输出类型
Sow = tbl.Sow;
X = tbl.X;
Y = tbl.Y;
Z = tbl.Z;
B = tbl.B;
L = tbl.L;
H = tbl.H;
E = tbl.E;
N = tbl.N;
U = tbl.U;
Gclock = tbl.Gclock;
Cclock = tbl.Cclock;
m_H = tbl.m_H;
m_V = tbl.m_V;
Vx = tbl.Vx;
Vy = tbl.Vy;
Vz = tbl.Vz;
Vc = tbl.Vc;
sigma0 = tbl.sigma0;
PDOP = tbl.PDOP;
GS = tbl.GS;
BS = tbl.BS;
end