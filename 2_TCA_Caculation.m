%% 对每个像元的三个产品时间序列土壤水分值，执行动态（101天的滑动窗口，步长为1）和静态TCA分析，并将每个产品的误差方差、R和fMSE计算结果输出
clc
clear

%% 读取一个tif模板，提取像元数目
[Templa,~]=geotiffread('D:\mxx\Rescale_sm_tiff_ByGLDAS_TemporalInterplo\asddesunion\smap\20151201.tif');
len=size(Templa,1)*size(Templa,2); %像元数目

%% 针对每个像元i,执行动态&静态TCA分析
parfor i=1:len
    fullpath=strcat('D:\mxx\Graduation\Fusion\SMAP_ASCAT_GLDAS\all_union_mat\',num2str(i),'.mat'); %% Data_prepare输出结果
    ID=load(fullpath);
    data_inter=ID.data;

    %% 剔除0-1范围外的异常值
    
    id_1=find(0<=data_inter(1,:));
    id_11=find(data_inter(1,:)<=1);    %smap
    id_c1=intersect(id_1,id_11);

    id_2=find(data_inter(2,:)>=0);      %ascat
    id_22=find(data_inter(2,:)<=1);
    id_c2=intersect(id_2,id_22);
       
    id_3=find(0<=data_inter(3,:));
    id_33=find(data_inter(3,:)<=1);   %gldas
    id_c3=intersect(id_3,id_33);

   %% 时间序列id，20150331-20201231依次对应1-m
   id=intersect(intersect(id_c1,id_c2),id_c3);     
   DATA=data_inter(:,id);
   DATA(4,:)=id;   

   %% Dynamic TCA， window size=101, step=1
   
    if size(DATA,2)>=30   %% 只有在时序长度大于等于30时，即至少存在30天土壤水分观测值时，才执行TCA计算，以确保结果的可靠性
         count=0;
         %% FMSE R ERROR初始化，第1列为时间序列id，第2列为滑动窗口内的可用样本数目，后3列依次为SMAP、ASCAT和GLDAS的TCA指标
         
         FMSE=nan*zeros(length(91:1:DATA(4,size(DATA,2))-90),5);
         R=nan*zeros(length(91:1:DATA(4,size(DATA,2))-90),5);
         ERROR=nan*zeros(length(91:1:DATA(4,size(DATA,2))-90),5);
         
         %% 对每个窗口中心的时刻，选取其前面50天和后面50天的土壤水分时序值，执行TCA
        for j=51:1:DATA(4,size(DATA,2))-50
            id_1=find(DATA(4,:)>=j-50,1,'first');
            id_2=find(DATA(4,:)<=j+50,1,'last');
            id_middle=find(DATA(4,:)==j);
            
            if isempty(id_middle)
               id_middle=nan;
            end
            IDp=id_middle-id_1+1;
            num=id_2-id_1+1;
            if num>=30
               count=count+1;
               x=DATA(1:3,id_1:id_2)'; 
               FMSE(count,1)=j;   %时间序列id
               R(count,1)=j;
               ERROR(count,1)=j;
               FMSE(count,2)=num; %滑动窗口内的可用样本数目
               R(count,2)=num;
               ERROR(count,2)=num;

               stats=TC(x);   % 调用TC函数
               FMSE(count,3:end)=stats(7:9)';   %TCA FMSE指标
               R(count,3:end)=stats(4:5)';      %TCA R指标
               ERROR(count,3:end)=stats(1:3)';  %TCA 误差方差指标
            end  
        end
        if count>=1   % 确保每个像元i存在合理的TCA计算结果才进行输出
            FMSE(count+1:end,:)=[]; 
            R(count+1:end,:)=[]; 
            ERROR(count+1:end,:)=[];
            %% 输出每个像元i的时序TCA计算结果
            savenames1=strcat("D:\mxx\Graduation\Fusion\SMAP_ASCAT_GLDAS\Output\Dynamic_W100\fMSE\",num2str(i),'.mat');
            parsave(savenames1,FMSE); 
            savename2=strcat("D:\mxx\Graduation\Fusion\SMAP_ASCAT_GLDAS\Output\Dynamic_W100\SNR\",num2str(i),'.mat');
            parsave(savename2,R); 
            savename3=strcat("D:\mxx\Graduation\Fusion\SMAP_ASCAT_GLDAS\Output\Dynamic_W100\Error\",num2str(i),'.mat');
            parsave(savename3,ERROR); 
        end
    end
    
%% Static TCA

   if size(DATA,2)>=30
     FMSE_Static=nan*zeros(1,3);
     R_Static=nan*zeros(1,3);
     ERROR_Static=nan*zeros(1,3);

      xx=DATA(1:3,:)'; 
      stats_static=TC(x);    
      FMSE_Static(1,:)=stats_static(7:9)';
      R_Static(1,:)=stats_static(4:6)';
      ERROR_Static(1,:)=stats_static(1:3)';
         
      savenames4=strcat("D:\mxx\Graduation\Fusion\SMAP_ASCAT_GLDAS\Output\Static\FMSE\",num2str(i),"_",num2str(size(DATA,2)),".mat");
      parsave(savenames4,FMSE_Static);
      savenames5=strcat("D:\mxx\Graduation\Fusion\SMAP_ASCAT_GLDAS\Output\Static\SNR\",num2str(i),"_",num2str(size(DATA,2)),".mat");
      parsave(savenames5,R_Static);
      savenames6=strcat("D:\mxx\Graduation\Fusion\SMAP_ASCAT_GLDAS\Output\Static\ERROR\",num2str(i),"_",num2str(size(DATA,2)),".mat");
      parsave(savenames6,ERROR_Static);
   end

end    
    
