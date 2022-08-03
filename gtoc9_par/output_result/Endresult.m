
% /****************************************************************************
% * Copyright (MATLAB), 2020-2031 �廪��ѧ���캽��ѧԺ����ѧ�����ʵ����
% * ����: ����
% * �ļ���: FigureTemplate.m
% * ���ݼ�����
% * ��ͼģ�壬�����廪��ѧѧ����ʽ
% * �ļ���ʷ��
% * �汾��     ����         ����       ˵��
% * 01a       2020-07-13    ����     �������ļ�
% */
 

clc
clear all
close all

x=[1 2 3];
y=[2 3 4];
plot(x,y)

%�������ͺ���ɫ
% ����1.��ɫ�ַ�����'c', 'm', 'y', 'r', 'g', 'b', 'w',��'k'���ֱ��ʾ�࣬���ϣ��ƣ��죬�̣��׺ͺڡ�
% ����2.�����ַ����У�'-' Ϊʵ��, '--' Ϊ����, ':' Ϊ����, '-.' Ϊ������, ��'none' ��ʾ�������͡�
% ����3.�����ʽ��'+', 'o', '*',��'x' ������'s' ����������, 'd' ��������, '^' Ϊ��������, 'v' Ϊ��������,
%        '>' Ϊ��������, '<' Ϊ��������, 'p' Ϊ�������, 'h' Ϊ��������, none Ϊ���ñ�ǡ�

%����
title('\fontname{����}\fontsize{8}ʱ��\fontname{Times New Roman}\fontsize{8}(s)')
%λ�ü���С
set(gcf,'unit','centimeters','position',[10 5 7 5]);
set(gca,'Position',[.13 .18 .8 .70]);
%ͼ��
h= legend('data1','location','NorthEast');                                            
set(h,'FontName','Times New Roman','FontSize',8,'LineWidth',0.3,'FontWeight','normal')
%x��
xlabel('\fontname{����}\fontsize{8}ʱ��\fontname{Times New Roman}\fontsize{8}(s)');   
%y��
ylabel(['\fontname{����}\fontsize{8}���ٶ�\fontname{Times new roman}\fontsize{8}(g)']);
%������
set(gca,'FontName','Times New Roman','FontSize',8,'LineWidth',0.3);
axis([0,10,0,5]);