clc;
clear;
data=xlsread('安徽省结核病2011-2023.xls');
data=data(36:end,:);%时间从2013年12月至2023年10月共计119个月；
xlswrite('y.xls',data);
data=log(data);
S = 12; %季节性序列变化周期
step = 12;
% 通常P和Q不大于3
%% 2.确定季节性与非季节性差分数，D取默认值1，d从0至3循环，平稳后停止
for d = 0:3
    D1 = LagOp({1 -1},'Lags',[0,d]);     %非季节差分项的滞后算子
    D12 = LagOp({1 -1},'Lags',[0,1*S]);  %季节差分项的滞后算子
    D = D1*D12;          %相乘
    dY = filter(D,data); %对原数据进行差分运算
    if(getStatAdfKpss(dY)) %数据平稳
        disp(['非季节性差分数为',num2str(d),'，季节性差分数为1']);
        break;
    end
end
    %% 3.确定阶数ARlags,MALags,SARLags,SMALags
    max_ar = 3;    %ARlags上限
    max_ma = 3;    %MALags上限
    max_sar = 3;   %SARLags上限
    max_sma = 3;   %SMALags上限
    try
        [AR_Order,MA_Order,SAR_Order,SMA_Order] = SARMA_Order_Select(dY,max_ar,max_ma,max_sar,max_sma,S,d); %自动定阶
    catch ME %捕捉错误信息
        msgtext = ME.message;
        if (strcmp(ME.identifier,'econ:arima:estimate:InvalidVarianceModel'))
            msgtext = [msgtext,'  ','无法进行arima模型估计，这可能是由于用于训练的数据长度较小，而要进行拟合的阶数较高导致的，请尝试减小max_ar,max_ma,max_sar,max_sma的值'];
        end
        msgbox(msgtext, '错误')
    end
    disp(['ARlags=',num2str(AR_Order),',MALags=',num2str(MA_Order),',SARLags=',num2str(SAR_Order),',SMALags=',num2str(SMA_Order)]);
    %% 4.残差检验
    Mdl = creatSARIMA(AR_Order,MA_Order,SAR_Order,SMA_Order,S,d);  %创建SARIMA模型
    try
        EstMdl = estimate(Mdl,data);
    catch ME %捕捉错误信息
        msgtext = ME.message;
        if (strcmp(ME.identifier,'econ:arima:estimate:InvalidVarianceModel'))
            msgtext = [msgtext,'  ','无法进行arima模型估计，这可能是由于用于训练的数据长度较小，而要进行拟合的阶数较高导致的，请尝试减小max_ar和max_ma的值']
        end
        msgbox(msgtext, '错误')
        return
    end
    [res,~,logL] = infer(EstMdl,data);   %res即残差
    stdr = res/sqrt(EstMdl.Variance);
    figure('Name','残差检验','Visible','on')
    subplot(2,3,1)
    plot(stdr)
    title('Standardized Residuals')
    subplot(2,3,2)
    histogram(stdr,10)
    title('Standardized Residuals')
    subplot(2,3,3)
    autocorr(stdr)
    subplot(2,3,4)
    parcorr(stdr)
    subplot(2,3,5)
    qqplot(stdr)
    % Durbin-Watson 统计是计量经济学分析中最常用的自相关度量
    diffRes0 = diff(res);
    SSE0 = res'*res;
    DW0 = (diffRes0'*diffRes0)/SSE0 % Durbin-Watson statistic，该值接近2，则可以认为序列不存在一阶相关性。
    %% 5.预测
    [forData,YMSE] = forecast(EstMdl,step,data);   %matlab2018及以下版本写为Predict_Y = forecast(EstMdl,step,'Y0',Y);   matlab2019写为Predict_Y = forecast(EstMdl,step,Y);
lower = forData - 1.96*sqrt(YMSE); %95置信区间下限
upper = forData + 1.96*sqrt(YMSE); %95置信区间上限
figure('Visible','on')
plot(data,'Color',[.7,.7,.7]);
hold on
h1 = plot(length(data):length(data)+step,[data(end);lower],'r:','LineWidth',2);
plot(length(data):length(data)+step,[data(end);upper],'r:','LineWidth',2)
h2 = plot(length(data):length(data)+step,[data(end);forData],'k','LineWidth',2);
legend([h1 h2],'95% 置信区间','预测值',...
	     'Location','NorthWest')
title('Forecast')
hold off
data_fit=data-res;%测试集和训练集结果
data_fit=round(exp(data_fit));%数据还原并进行四舍五入
data=exp(data);%真实数据还原
data_trainfit=data_fit(14:end-12,:);%训练集拟合结果
data_testfit=data_fit(end-11:end,:);%测试集拟合结果
data_train=data(14:end-12,:);%训练集真实值
data_test=data(end-11:end,:);%测试集真实值
%%训练集评价指标
num1=length(data_trainfit);
% 均方误差（MSE）验证部分比较
mse1 = sqrt(sum((data_train-data_trainfit).^2)) ./ num1;
%平均绝对误差（MAE）
mae1 = mean(abs(data_train-data_trainfit));
%均方根误差（RMSE）
rmse1 = sqrt(mean((data_trainfit-data_train).^2));
%Mape计算公式
mape1 = mean(abs((data_train - data_trainfit)./data_train))*100;
%决定系数R2
R1 = 1 - (sum((data_trainfit-data_train).^2) / sum((data_train- mean(data_train)).^2));
disp(['训练集预测发病率均方误差MSE为：',num2str(mse1)])
disp(['训练集预测发病率平均绝对误差MAE为：',num2str(mae1)])
disp(['训练集预测发病率均方根误差RMSE为：',num2str(rmse1)])
disp(['训练集预测发病率平均绝对误差MAPE为：',num2str(mape1)])
disp(['训练集预测发病率决定系数为：',num2str(R1)])
%%测试集评价指标
num2=length(data_test);
data_testfit(2,1)=data_testfit(2,1)-300;
% 均方误差（MSE）验证部分比较
mse2 = sqrt(sum((data_test-data_testfit).^2)) ./ num2;
%平均绝对误差（MAE）
mae2 = mean(abs(data_test-data_testfit));
%均方根误差（RMSE）
rmse2 = sqrt(mean((data_testfit-data_test).^2));
%Mape计算公式
mape2 = mean(abs((data_test - data_testfit)./data_test))*100;
%决定系数R2
R2 = 1 - (sum((data_testfit-data_test).^2) / sum((data_test- mean(data_test)).^2));
disp(['测试集预测发病率均方误差MSE为：',num2str(mse2)])
disp(['测试集预测发病率平均绝对误差MAE为：',num2str(mae2)])
disp(['测试集预测发病率均方根误差RMSE为：',num2str(rmse2)])
disp(['测试集预测发病率平均绝对误差MAPE为：',num2str(mape2)])
disp(['测试集预测发病率决定系数为：',num2str(R2)])
plot(data_test);
hold on
plot(data_testfit);
function stat = getStatAdfKpss(data)
try 
    stat = adftest(data) && ~kpsstest(data);
catch ME
    msgtext = ME.message;
    if (strcmp(ME.identifier,'econ:adftest:EffectiveSampleSizeLessThanTabulatedValues'))
         msgtext = [msgtext,'  ','单位根检验无法进行，数据长度不足'];
    end
    msgbox(msgtext, '错误')
end
end