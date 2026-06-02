%% 
%% 输入基于TCA_caculation的每个像元的动态TCA分析结果，转换为tif格式
clc
clear

file=dir("D:\mxx\Graduation\Fusion\SMAP_ASCAT_GLDAS\Output\Dynamic_W100\fMSE\*.mat");  %以FMSE指标为例
filename=struct2cell(file);
filenames=filename(1,:);
m=length(filenames);
%% 根据研究时间范围设置时间索引
startDate = datenum(2015,3,31);
endDate = datenum(2020,12,31);
dateSequence = startDate:endDate;
vec = datevec(dateSequence);
dateDouble = vec(:,1)*10000 + vec(:,2)*100 + vec(:,3);

id_tiff=[];
date=[];
num=[];
data_median=[];  

parfor i=1:m
    str_name=filenames{1,i}(1:length(filenames{1,i})-4);
    fullpath_boot=strcat('D:\mxx\Graduation\Fusion\SMAP_ASCAT_GLDAS\Output\Dynamic_W100\fMSE\',filenames{1,i});
    ID=load(fullpath_boot);
    data=ID.DATA;

    id=size(data,1); %像元id
    data_median=[data_median;data(:,3:5)];  %TCA指标
    num=[num;data(:,2)];
    date=[date;dateDouble(data(:,1))];%时序id
    id_tiff=[id_tiff;repmat(str2double(str_name),id,1)];
    
end    

[date_sort,I]=sort(date);
clear date file filename filenames str_name m 
id_tiff=id_tiff(I,:);
num=num(I,:);
data_median=data_median(I,:);                            
[~,row,~]=unique(date_sort);
id=[diff(row);length(date_sort)+1-row(length(row))];
data_median_cell=mat2cell(data_median,id,1);        
clear data_median row
id_tiff_cell=mat2cell(id_tiff,id,1);
clear id_tiff
num_cell=mat2cell(num,id,1);
clear num
date_name_cell=mat2cell(date_sort,id,1);
clear date_sort
cellfun(@(x,y,n,p) mattotiff_MedianMetrics_ForStepSize(x,y,n,p),data_median_cell,id_tiff_cell,num_cell,date_name_cell,'UniformOutput',false);

%tif转换函数
function [] = mattotiff_MedianMetrics_ForStepSize(stat,tiffid,sample,date_name)
      %% 读取tif模板
      [Templa,T]=geotiffreadgeotiffread('D:\mxx\Rescale_sm_tiff_ByGLDAS_TemporalInterplo\asddesunion\smap\20151201.tif');
      Templa(:,:)=nan;
      %% 设置输出tif的路径
      str="D:\mxx\Graduation\Fusion\SMAP_ASCAT_GLDAS\toTiff\Dynamic_W100";
      
      %% 输出每个像元的可用样本数目，需要预先在输出路径下新建samples文件夹
      Samples=Templa;
      Samples(tiffid)=sample;
      savename_samples=strcat(str,"\samples\",num2str(date_name(1)),".tif");
      geotiffwrite(savename_samples,Samples,T);     
    
      %% 输出每个产品的fMSE指标，需要预先在输出路径下新建fMSE_SMAP、fMSE_ASCAT和fMSE_GLDAS文件夹
      R_unrescale=Templa;
      R_unrescale(tiffid)=stat(:,1);
      savename_R_Lc1_unrescale=strcat(str,"\fMSE_SMAP\",num2str(date_name(1)),".tif");
      geotiffwrite(savename_R_Lc1_unrescale,R_unrescale,T);

      R_unrescale=Templa;
      R_unrescale(tiffid)=stat(:,2); 
      savename_CI_R_IC_unrescale=strcat(str,"\fMSE_ASCAT\",num2str(date_name(1)),".tif");
      geotiffwrite(savename_CI_R_IC_unrescale,R_unrescale,T); 

      R_unrescale=Templa;
      R_unrescale(tiffid)=stat(:,3);
      savename_R_ASCAT_unrescale=strcat(str,"\fMSE_GLDAS\",num2str(date_name(1)),".tif");
      geotiffwrite(savename_R_ASCAT_unrescale,R_unrescale,T);
      

end

