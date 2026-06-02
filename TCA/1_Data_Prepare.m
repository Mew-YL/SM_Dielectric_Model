%% 同时读取20150331-20201231时间范围内的被动-主动-模型三个土壤水分产品，产品格式为tif，将每个像元的三个产品的时间序列整理成一个mat文件并输出
%% 要求三个产品的时序相同，tif文件均为WGS84坐标系且行列数相同
%% 考虑到SMAP数据在2019年部分时间存在缺失，需事先准备空值tif以补全时间序列
clc;
clear;

file=dir("D:\mxx\Rescale_sm_tiff_ByGLDAS_TemporalInterplo\asddesunion\gldas\*.tif");
filename=struct2cell(file);clear file
filenames=filename(1,:);clear filename
m=length(filenames); %时间序列长度

%% 读取一个tif模板，提取像元数目
[Templa,~]=geotiffread('D:\mxx\Rescale_sm_tiff_ByGLDAS_TemporalInterplo\asddesunion\smap\20151201.tif');
len=size(Templa,1)*size(Templa,2); %像元数目
DATA1=nan*zeros(len,m);
DATA2=nan*zeros(len,m);
DATA3=nan*zeros(len,m);

%% 开启并行计算
if isempty(gcp('nocreate'))
    parpool local;
end   

%% 针对每一时刻i，以SMAP-ASCAT-GLDAS三联体为例，读取三个产品的原始tif文件，并分别整理成三个维度为len*m的矩阵
parfor i=1:m
        fullpath=strcat('D:\mxx\Rescale_sm_tiff_ByGLDAS_TemporalInterplo\asddesunion\smap\',filenames{1,i});
        [smap,~]=geotiffread(fullpath);
        DATA1(:,i)=reshape(smap,len,1); 
        %% ASCAT已经过单位转换(% to cm3/cm3)，取值范围为0-1
        fullpath=strcat('D:\mxx\Rescale_sm_tiff_ByGLDAS_TemporalInterplo\asddesunion\ascat\',filenames{1,i});
        [ascat,~]=geotiffread(fullpath);
        DATA2(:,i)=reshape(ascat,len,1); 
        %% GLDAS已经过单位转换(kg/m2 to cm3/cm3)，取值范围为0-1
        fullpath=strcat('D:\mxx\Rescale_sm_tiff_ByGLDAS_TemporalInterplo\asddesunion\gldas\',filenames{1,i});
        [gldas,~]=geotiffread(fullpath);
        DATA3(:,i)=reshape(gldas,len,1); 

end

%% 针对每一个像元j，提取三个产品的时间序列土壤水分值，以像元id为命名，整理成mat文件输出
for j=1:len
     data=nan*zeros(3,m);
     data(1,:)=DATA1(j,:);
     data(2,:)=DATA2(j,:);
     data(3,:)=DATA3(j,:);
   savename=strcat("D:\mxx\Graduation\Fusion\SMAP_ASCAT_GLDAS\all_union_mat\",num2str(j),".mat");
   parsave(savename,data); 
end
