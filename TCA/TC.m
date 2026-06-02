%% TCA分析
function stats= TC(x)
%%% 1:SMAP 2:ASCAT 3:GLDAS
a=x(:,1); % SMAP时间序列
b=x(:,2); % ASCAT时间序列
c=x(:,3); % GLDAS时间序列
%% 以GLDAS为基准，对SMAP和ASCAT进行尺度缩放，将三个产品的时间序列土壤水分统一至相同数据空间
rescal_a=cov_calc(c,b)/cov_calc(a,b);
rescal_b=cov_calc(c,a)/cov_calc(b,a);
%% 尺度缩放后的SMAP和ASCAT
a_c=rescal_a*(a-mean(a))+mean(c);
b_c=rescal_b*(b-mean(b))+mean(c);

%% TCA分析，后三个元素分别为SMAP、ASCAT和GLDAS的误差方差
y_rescale=[cov_calc(a_c,b_c)*cov_calc(a_c,c)/cov_calc(b_c,c),cov_calc(b_c,a_c)*cov_calc(b_c,c)/cov_calc(a_c,c),cov_calc(c,a_c)*cov_calc(c,b_c)/cov_calc(a_c,b_c),var(a_c)-cov_calc(a_c,b_c)*cov_calc(a_c,c)/cov_calc(b_c,c),var(b_c)-cov_calc(b_c,a_c)*cov_calc(b_c,c)/cov_calc(a_c,c),var(c)-cov_calc(c,a_c)*cov_calc(c,b_c)/cov_calc(a_c,b_c)];

%% 基于y_rescale进一步计算R和fMSE指标

R1=sqrt_define(y_rescale(1)/var(a));
R2=sqrt_define(y_rescale(2)/var(b));
R3=sqrt_define(y_rescale(3)/var(c));

fmse1=sqrt_define(y_rescale(4)/var(a));
fmse2=sqrt_define(y_rescale(5)/var(b));
fmse3=sqrt_define(y_rescale(6)/var(c));

%% 输出TCA计算结果：[ErrorVariance_SMAP, Error   Variance_ASCAT, ErrorVariance_GLDAS, R_SMAP, R_ASCAT, R_GLDAS,fMSE_SMAP, fMSE_ASCAT, fMSE_GLDAS]
stats=[y_rescale(4:6),R1,R2,R3,fmse1,fmse2,fmse3]';
end
function stat= cov_calc(x,y)
stat_inter=cov(x,y);
stat=stat_inter(1,2);
end
function y = sqrt_define(x)
         if x>=0
             y=sqrt(x);
         else
             y=sqrt(-x);
             y=-y;
         end
end
