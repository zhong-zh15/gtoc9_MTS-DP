
% /****************************************************************************
% * Copyright (MATLAB), 2020-2031 清华大学航天航空学院动力学与控制实验室
% * 作者: 张众
% * 文件名: FigureTemplate.m
% * 内容简述：
% * 绘图模板，仿照清华大学学报格式
% * 文件历史：
% * 版本号     日期         作者       说明
% * 01a       2020-07-13    张众     创建该文件
% */
 

clc
clear all
close all

x=[1 2 3];
y=[2 3 4];
plot(x,y)

%常见线型和颜色
% 　　1.颜色字符串有'c', 'm', 'y', 'r', 'g', 'b', 'w',和'k'。分别表示青，红紫，黄，红，绿，白和黑。
% 　　2.线型字符串有：'-' 为实线, '--' 为虚线, ':' 为点线, '-.' 为点虚线, 及'none' 表示不用线型。
% 　　3.标记形式有'+', 'o', '*',和'x' ，填入's' 代表正方形, 'd' 代表菱形, '^' 为上三角形, 'v' 为下三角形,
%        '>' 为右三角形, '<' 为左三角形, 'p' 为五角星形, 'h' 为六角星形, none 为不用标记。

%标题
title('\fontname{宋体}\fontsize{8}时间\fontname{Times New Roman}\fontsize{8}(s)')
%位置及大小
set(gcf,'unit','centimeters','position',[10 5 7 5]);
set(gca,'Position',[.13 .18 .8 .70]);
%图例
h= legend('data1','location','NorthEast');                                            
set(h,'FontName','Times New Roman','FontSize',8,'LineWidth',0.3,'FontWeight','normal')
%x轴
xlabel('\fontname{宋体}\fontsize{8}时间\fontname{Times New Roman}\fontsize{8}(s)');   
%y轴
ylabel(['\fontname{宋体}\fontsize{8}加速度\fontname{Times new roman}\fontsize{8}(g)']);
%坐标轴
set(gca,'FontName','Times New Roman','FontSize',8,'LineWidth',0.3);
axis([0,10,0,5]);