function data = importfile(filename, dataLines)
%IMPORTFILE 从文本文件中导入数据
%  DATA = IMPORTFILE(FILENAME)读取文本文件 FILENAME 中默认选定范围的数据。  返回数值数据。
%
%  DATA = IMPORTFILE(FILE, DATALINES)按指定行间隔读取文本文件 FILENAME
%  中的数据。对于不连续的行间隔，请将 DATALINES 指定为正整数标量或 N×2 正整数标量数组。
%
%  示例:
%  data = importfile("D:\GitHub\SPP_Design\报告\单频双系统.pos", [1, Inf]);
%
%  另请参阅 READTABLE。
%
% 由 MATLAB 于 2023-11-20 10:36:44 自动生成

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
opts.VariableNames = ["Var1", "Var2", "VarName3", "Var4", "VarName5", "VarName6", "VarName7", "Var8", "VarName9", "VarName10", "VarName11", "Var12", "Var13", "Var14", "Var15", "Var16", "VarName17", "Var18", "VarName19", "Var20", "Var21", "Var22", "Var23", "Var24", "Var25", "Var26", "Var27", "Var28", "Var29", "VarName30", "Var31", "Var32", "Var33", "VarName34", "Var35", "VarName36", "Var37", "VarName38", "Var39"];
opts.SelectedVariableNames = ["VarName3", "VarName5", "VarName6", "VarName7", "VarName9", "VarName10", "VarName11", "VarName17", "VarName19", "VarName30", "VarName34", "VarName36", "VarName38"];
opts.VariableTypes = ["string", "string", "double", "string", "double", "double", "double", "string", "double", "double", "double", "string", "string", "string", "string", "string", "double", "string", "double", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "double", "string", "string", "string", "double", "string", "double", "string", "double", "string"];

% 指定文件级属性
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% 指定变量属性
opts = setvaropts(opts, ["Var1", "Var2", "Var4", "Var8", "Var12", "Var13", "Var14", "Var15", "Var16", "Var18", "Var20", "Var21", "Var22", "Var23", "Var24", "Var25", "Var26", "Var27", "Var28", "Var29", "Var31", "Var32", "Var33", "Var35", "Var37", "Var39"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var1", "Var2", "Var4", "Var8", "Var12", "Var13", "Var14", "Var15", "Var16", "Var18", "Var20", "Var21", "Var22", "Var23", "Var24", "Var25", "Var26", "Var27", "Var28", "Var29", "Var31", "Var32", "Var33", "Var35", "Var37", "Var39"], "EmptyFieldRule", "auto");

% 导入数据
data = readtable(filename, opts);

%% 转换为输出类型
data = table2array(data);
end