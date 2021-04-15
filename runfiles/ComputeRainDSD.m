function [] = ComputeRainDSD(path)
% plot DSD of rainfall
% LNw,Dm,Mu

% Modified in April 1,2021, by yaru zhang

DisDSDPath=strcat(path,'UK\DSD\');

DisDSD_Count=load(strcat(DisDSDPath,'Chi_DSDP_2003_2017_1m.dat'));
DisDSD_Size=load(strcat(DisDSDPath,'Chi_DSD_Size'));
GaugeRain=load(strcat(DisDSDPath,'Chi_RainP_2003_2017_1m.dat'));

[OutputData] = Root_ComputeRainDSD(DisDSD_Count,DisDSD_Size,GaugeRain);

WriteGridData(DisDSDPath,'Chi_DSDInfos_2003_2017_1m.dat',OutputData,3);

end

function [OutputData] = Root_ComputeRainDSD(DisDSD_Count,DSD_Size,GaugeRain)

Sensor_A = 50;  % unit: cm^2
TimeCols = 6;
GaugeIndex = 4;
dt = 60;  % unit: second

[TotalTimes, BinsNum] = size(DisDSD_Count);
BinsNum = BinsNum - TimeCols;
Date_Given = GaugeRain(:,1:TimeCols);
Rain_Given = GaugeRain(:,TimeCols+GaugeIndex);
R_gag = Rain_Given* 3600/dt;
[dD] = ComputedD(DSD_Size);

R_Dis=zeros(TotalTimes,1);
DSD_Params=zeros(TotalTimes,3);
Nm=zeros(TotalTimes,BinsNum);
for t=1:TotalTimes
    curDisDSD_Count=DisDSD_Count(t,TimeCols+1:end)';
    vol = Rain_Given(t);
    if vol>0
        [Nm(t,:)] = ComputeNm(Sensor_A, dt, dD, curDisDSD_Count,DSD_Size);
        [DSD_Params_O] = ComputeDistParams(Nm(t,:), dD, DSD_Size);
        if min(DSD_Params_O)>-999
            DSD_Params(t,:) = DSD_Params_O;
            [N_D] = ComputeDisND(DSD_Params_O,DSD_Size);
            [R_Dis(t)] = ComputeDisRain(dD, N_D ,DSD_Size);
        end
    end
end

OutputData=zeros(TotalTimes,TimeCols+5);
OutputData(:,1:TimeCols) = Date_Given;
OutputData(:,TimeCols+1) = R_gag;
OutputData(:,TimeCols+2) = R_Dis;
OutputData(:,TimeCols+3:TimeCols+5) = DSD_Params;
OutputData = GetPositiveParams(OutputData, 0);

% [R_gag_N,R_dis_N] = ComputeHourlyRainfall(R_gag,R_Dis,Date_Given);

Dm = GetPositiveData(DSD_Params(:,1), 0);
Nw = GetPositiveData(DSD_Params(:,2), 0);
Mu = GetPositiveData(DSD_Params(:,3), -10);
LNw = log10(Nw);

[Nm_binM] = GetPositiveMean(Nm);

ylim_v=120;
plot(R_gag, R_Dis, 'k.');
hold on;
plot([0,ylim_v],[0,ylim_v],'k-');
xlabel('Measured rainfall (mm/h)');
ylabel('Model rainfall (mm/h)');
xlim([0,ylim_v]);
ylim([0,ylim_v]);

figure;
DSD_Size(1)=[];
Nm_binM(1)=[];
plot(DSD_Size, Nm_binM, 'k-');
xlabel('Drop diameter (mm)');
ylabel('Concentration (mm^-^1 m^-^3)');

figure;
ksdensity(LNw);
xlabel('Log_1_0 N_w, N_w (mm^-^1 m^-^3)');
ylabel('PDF');

figure;
ksdensity(Dm);
xlabel('D_m (mm)');
ylabel('PDF');

figure;
ksdensity(Mu);
xlabel('\mu');
ylabel('PDF');

end

function [R_gag_N,R_dis_N] = ComputeHourlyRainfall(R_gag,R_Dis,Dates)

[TotalTimes,~]=size(Dates);

lastDate = Dates(1,:);
acc_gauge = 0;
acc_dis = 0;
acc_n = 0;
R_gag_N = [];
R_dis_N = [];
N_num = 0;
for t=1:TotalTimes
    curDate = Dates(t,:);
    
    date_dif = curDate(1:4)-lastDate(1:4);
    if sum(date_dif.^2)>0
        N_num=N_num+1;
        R_gag_N(N_num)=acc_gauge;
        if acc_n>0
            R_dis_N(N_num)=acc_dis/acc_n;
        else
            R_dis_N(N_num)=0;
        end
        
        acc_gauge = 0;
        acc_dis = 0;
        acc_n = 0;
        lastDate = curDate;
    end
    
    acc_gauge=acc_gauge+R_gag(t);
    acc_dis=acc_dis+R_Dis(t);
    if R_Dis(t)>0
        acc_n=acc_n+1;
    end
    
end

end

function [DParams] = ComputeDistParams(Nm, dD,DSD_Size)

[dsn] = length(DSD_Size);

m4 = 0;
m3 = 0;
m2 = 0;
m6 = 0;
for i=1:dsn
    Di= DSD_Size(i);
    Nmi = Nm(i);
    if Nmi>0     
        fun4 = Di^4 * Nmi * dD(i);
        m4 = m4 +fun4;

        fun3 = Di^3 * Nmi * dD(i);
        m3 = m3 +fun3;

        fun2 = Di^2 * Nmi * dD(i);
        m2 = m2 +fun2;

        fun6 = Di^6 * Nmi * dD(i);
        m6 = m6 +fun6;
    end
end

Dm = m4/m3;
Nw = 256/6*(m3^5/m4^4);
yt = m4^2/(m2*m6);
Mu = ((7-11*yt)-sqrt((7-11*yt)^2-4*(yt-1)*(30*yt-12)))/(2*(yt-1));

DParams=-999*ones(1,3);
if Dm>=0 && Nw>=0 && Mu>=-10 && Mu<20
    DParams(1)=Dm;
    DParams(2)=Nw;
    DParams(3)=Mu;
end

end

function [N_D] = ComputeDisND(DParams,DSD_Size)

[dsn] = length(DSD_Size);

Dm = DParams(1);
Nw = DParams(2);
u = DParams(3);

N_D = zeros(dsn,1);
for i=1:dsn
    Di = DSD_Size(i);
    fu = (6 * (4+u)^(u+4))/(4^4 * gamma(u+4));
    N_D(i) = Nw * fu * (Di/Dm)^u * exp(-1 * (u+4) * (Di/Dm));
end

end

function [DX] = GetPositiveData(DO, tol)

dn = length(DO);

DX=[];
acc_n=0;
for i=1:dn
    if DO(i)>tol
        acc_n=acc_n+1;
        DX(acc_n)=DO(i);
    end
end

end

function [Params_N] = GetPositiveParams(Params, tol)

[dn,~] = size(Params);

Params_N=[];
acc_n=0;
for i=1:dn
    if Params(i,end-1)>tol
        acc_n=acc_n+1;
        Params_N(acc_n,:)=Params(i,:);
    end
end

end

function [Data_M] = GetPositiveMean(Data)

[nrows, ncols]=size(Data);

Data_M=zeros(1,ncols);
Data_New=[];
for i=1:ncols
    data_one = Data(:,i);
    data_one_positive = data_one(data_one>0);
    Data_M(i)=mean(data_one_positive);
end

end

function [DropSizes] = GetTimeSeriesData(DSD_Count,DSD_Size)

classnum = length(DSD_Count);

DropSizes=[];
acc_num=0;
for i=1:classnum
    curN=DSD_Count(i);
    if curN>0
        curS=DSD_Size(i);
        for j=1:curN
            acc_num=acc_num+1;
            DropSizes(acc_num)=curS;
        end
    end
end
end

function [R] = ComputeDisRain(dD, N_D, DSD_Size)

[dsn] = length(DSD_Size);

R_coef = pi * 6 / 10^4;

acc_rain = 0;
for i=1:dsn
    Di = DSD_Size(i);
    Ni = N_D(i);
    Vi = ComputeDropVel_AtlV(Di);
    acc_rain = acc_rain + Di^3 * Vi * Ni * dD(i);
end

R = R_coef * acc_rain;

end

function [Nm] = ComputeNm(Sensor_A, dt, dD, DSD_Count,DSD_Size)

[dsn] = length(DSD_Size);

Nm = zeros(1, dsn);
for i=1:dsn
    Di = DSD_Size(i);
    Ni = DSD_Count(i);
    Vi = ComputeDropVel_AtlV(Di);
    Nm(i) = Ni / (Sensor_A * 10^(-4) * dt * Vi * dD(i));
end

end

function [dD] = ComputedD(DSD_Size)

[dsn] = length(DSD_Size);
dD = zeros(1,dsn);

dD(1)=DSD_Size(1);
for i=2:dsn
    dD(i) = DSD_Size(i)-DSD_Size(i-1);
end

end

function [V_Atl] = ComputeDropVel_AtlV(Di)

V_Atl = 3.78 * Di^0.67;

end

function [V_Atl] = ComputeDropVel_Atl(Di)

V_Atl = 17.67 * (Di / 10)^0.67;

end

function [V_Upl] = ComputeDropVel_Upl(Di)

V_Upl = 4.874 * Di * exp(-0.195 * Di);

end

function [V_VD] = ComputeDropVel_VD(Di)

V_VD = 0.0561 * Di^3 - 0.912 * Di^2 + 5.03 * Di - 0.254;

end

function [] = WriteGridData(path,f_name,resultData,type)

new_path=strcat(path,f_name);
% dlmwrite(new_path,resultData,'-append','delimiter',' ');
if type==0 % integer
    dlmwrite(new_path,resultData,'delimiter','\t','precision','%.6g');
elseif type==1  % float with 2 precision
    dlmwrite(new_path,resultData,'delimiter','\t','precision','%.2f');
elseif type==2  % float with 4 precision
    dlmwrite(new_path,resultData,'delimiter','\t','precision','%.4f');
elseif type==3  % float with 4 precision
    dlmwrite(new_path,resultData,'delimiter','\t');
end

end