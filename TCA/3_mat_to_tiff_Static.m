%% 
%% 输入基于TCA_caculation的每个像元的静态TCA分析结果，转换为tif格式
clc
clear

file=dir("D:\mxx\Graduation\Fusion\SMAP_ASCAT_GLDAS\Output\Static\FMSE\*.mat"); %以FMSE指标为例
filename=struct2cell(file);
filenames=filename(1,:);
m=length(filenames);

str_name=cell(1,2);
id_tiff=nan*zeros(m,1);
num=nan*zeros(m,1);
data_metric=nan*zeros(m,3);   

parfor i=1:m
    str_name=strsplit(filenames{1,i}(1:length(filenames{1,i})-4),'_');
    id_tiff(i)=str2double(str_name{1});  %像元id
    num(i)=str2double(str_name{2});      %可用样本数目
    fullpath_boot=strcat('D:\mxx\Graduation\Fusion\SMAP_ASCAT_GLDAS\Output\Static\FMSE\',filenames{1,i});
    ID=load(fullpath_boot);
    data=ID.DATA;
    data_metric(i,:)=data;   %TCA指标

end    
mat2tiff_ForStatic(data_metric,id_tiff,num);  %调用tif转换函数


function [] = mat2tiff_ForStatic(stat,tiffid,sample)
      %% 读取tif模板
      [Templa,T]=geotiffread('D:\mxx\Rescale_sm_tiff_ByGLDAS_TemporalInterplo\asddesunion\smap\20151201.tif');
      Templa(:,:)=nan;
      
      %% 设置输出tif的路径
      str="D:\mxx\Graduation\Fusion\SMAP_ASCAT_GLDAS\toTiff\Static_TCall";
      
      %% 输出每个像元的可用样本数目
      Samples=Templa;
      Samples(tiffid)=sample;
      savename_samples=strcat(str,"\samples.tif");
      geotiffwrite(savename_samples,Samples,T); 
      
 %% 输出每个产品的fMSE指标
 
      R_unrescale=Templa;
      R_unrescale(tiffid)=stat(:,1); 
      savename_CI_R_IC_unrescale=strcat(str,"\fMSE_SMAP.tif");
      geotiffwrite(savename_CI_R_IC_unrescale,R_unrescale,T); 
      
      R_unrescale=Templa;
      R_unrescale(tiffid)=stat(:,2); 
      savename_CI_R_IC_unrescale=strcat(str,"\fMSE_ASCAT.tif");
      geotiffwrite(savename_CI_R_IC_unrescale,R_unrescale,T); 
      
      R_unrescale=Templa;
      R_unrescale(tiffid)=stat(:,3); 
      savename_CI_R_IC_unrescale=strcat(str,"\fMSE_GLDAS.tif");
      geotiffwrite(savename_CI_R_IC_unrescale,R_unrescale,T); 

end

