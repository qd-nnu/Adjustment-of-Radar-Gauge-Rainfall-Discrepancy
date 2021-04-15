function [] = RME_RunModel(path,events)
% STEP 3 of RME model
% Simulate the raindrop microphysical process for each event

mode=2;   % 1-consider drift, 2-consider evp, 3-consider both
RDs = 0.1:0.1:6;
RadarParams = [479500,330500,1000,1000,100,100];

path=strcat(path,'hameldon_hill\');
MOutputPath=strcat(path,'Model_Output\');

StormInfo=GetStormInfo(strcat(path,events));

Root_runRME(MOutputPath,StormInfo,2,RDs,RadarParams);
Root_runRME(MOutputPath,StormInfo,3,RDs,RadarParams);

end

function [] = Root_runRME(MOutputPath,StormInfo,mode,RDs,RadarParams)

AtmosPre='Atmos';
SavePathName = char('Drift','Evp','DriftEvp');

% read model information
LayersHeight = load(strcat(MOutputPath,'Layers_Height.txt'));
DomainNet = load(strcat(MOutputPath,'Domain_Net.txt'));
Dv_dp_Table = load(strcat(MOutputPath,'Dv_dp_Table.txt'));
GridNets = load(strcat(MOutputPath,'GridNets.txt'));

[TotalEvents,~]=size(StormInfo);

[DropsInit] = ComputeDropsizeInfo(GridNets,RDs);
[DropsInfo_Init] = DropsInitiation(DropsInit,LayersHeight,DomainNet);

for e=1:TotalEvents  
    eventfolder=sprintf('%.4d-%.2d-%.2d',StormInfo(e,1:3));  
    EventPath = strcat(MOutputPath,eventfolder,'/'); 
    AtmosPath = strcat(EventPath,AtmosPre,'/'); 
    
    SaveDataPath=strcat(MOutputPath,eventfolder,'/',SavePathName(mode,:),'/');      

    if ~exist(SaveDataPath)
        mkdir(SaveDataPath);
    end 
    
    [E_dur] = GetEventDurationByFile(AtmosPath,AtmosPre);
    
    for t=1:E_dur
        AtmosData = load(strcat(AtmosPath,AtmosPre,num2str(t)));
        fprintf('Read atmos data is computed...\n');

        [DropsSum] = SimRainMicroEvolution(mode,DropsInfo_Init,AtmosData,LayersHeight,DomainNet,Dv_dp_Table);
        
        switch mode
            case 1
                [GroundDrift] = GetDropGroundDrift(RDs,RadarParams,DropsSum);
                WriteGridData(SaveDataPath,strcat('Drift',num2str(t)),GroundDrift,3);
                
            case 2
                [GroundEvp] = GetDropGroundEvp(RDs,DropsSum);
                WriteGridData(SaveDataPath,strcat('Evp',num2str(t)),GroundEvp,4);
                
            case 3
                [GroundDrift] = GetDropGroundDrift(RDs,RadarParams,DropsSum);
                [GroundEvp] = GetDropGroundEvp(RDs,DropsSum);
                
                WriteGridData(SaveDataPath,strcat('Drift',num2str(t)),GroundDrift,3);
                WriteGridData(SaveDataPath,strcat('Evp',num2str(t)),GroundEvp,4);
        end
        fprintf('Event %s Time %d is computed...\n',eventfolder,t);
    end
end

end

function [DropsSum] = SimRainMicroEvolution(mode,DropsInfo_Init,AtmosData,LayersHeight,DomainNet,Dv_dp_Table)

global IDrop_S IDrop_DI IDrop_LI IDrop_DX IDrop_DY IDrop_DZ IDrop_AX IDrop_AY;
global IAtmo_VX IAtmo_VY IAtmo_VZ IAtmo_T IAtmo_P IAtmo_RH;

% constant variable
%------------------------
dom_size=5000;  % WRF domain size
dt = 60;  % unit: seconds
t_start = 0;
t_end = 1800;
dom_x_num=35;
dom_y_num=35;
%------------------------

iteration = (t_end-t_start)/dt;
dom_xy_num = dom_x_num*dom_y_num;

% drop size, domain_id, layer_id, dx, dy, dz
IDrop_S=1;IDrop_DI=2;IDrop_LI=3;
IDrop_DX=4;IDrop_DY=5;IDrop_DZ=6;
IDrop_AX=7;IDrop_AY=8;
% domain_id: 1-(x1,y1),2-(x2,y1),3-(x3,y1),...n-(x1,y2),...
% dz means the distance to the bottom of the layer
% dx, dy mean the distance to the center of the domain

IAtmo_VX=1;IAtmo_VY=2;IAtmo_VZ=3;IAtmo_T=4;IAtmo_P=5;IAtmo_RH=6;

[Drops_num, ~] = size(DropsInfo_Init);
DropsInfo = DropsInfo_Init;

switch mode
    case 1
        for i=1 : iteration  
            for d=1:Drops_num
                dInfo = DropsInfo(d,:);
                if dInfo(IDrop_S)>0 && dInfo(IDrop_LI)>0
                    [dAtmos] = GetDropAtmos(dInfo,AtmosData,dom_xy_num);
                    [dInfo] = DropMoveTick(dInfo,dAtmos, dt);
                    [dInfo] = SetDropDomain(dInfo,LayersHeight,dom_size,dom_x_num);

                    DropsInfo(d,:)=dInfo;
                end
            end
            fprintf('...Iteration %d is computed...\n',i);
        end
    case 2
        for i=1 : iteration  
            for d=1:Drops_num
                dInfo = DropsInfo(d,:);
                if dInfo(IDrop_S)>0 && dInfo(IDrop_LI)>0
                    [dAtmos] = GetDropAtmos(dInfo,AtmosData,dom_xy_num);
                    [dInfo] = DropEvaporationTick(dInfo,dAtmos, dt, Dv_dp_Table);

                    if dInfo(IDrop_S)>0   % check as the drop may evaporate
                        [dInfo] = DropMoveTick_NoDrift(dInfo, dt);
                        [dInfo] = SetDropDomain(dInfo,LayersHeight,dom_size,dom_x_num);
                    end
                    DropsInfo(d,:)=dInfo;
                end
            end
            fprintf('...Iteration %d is computed...\n',i);
        end
    case 3
        for i=1 : iteration  
            for d=1:Drops_num
                dInfo = DropsInfo(d,:);
                if dInfo(IDrop_S)>0 && dInfo(IDrop_LI)>0
                    [dAtmos] = GetDropAtmos(dInfo,AtmosData,dom_xy_num);
                    [dInfo] = DropEvaporationTick(dInfo,dAtmos, dt, Dv_dp_Table);

                    if dInfo(IDrop_S)>0   % check as the drop may evaporate
                        [dInfo] = DropMoveTick(dInfo,dAtmos, dt);
                        [dInfo] = SetDropDomain(dInfo,LayersHeight,dom_size,dom_x_num);
                    end
                    DropsInfo(d,:)=dInfo;
                end
            end
            fprintf('...Iteration %d is computed...\n',i);
        end
end

DropsSum=[DropsInfo(:,1),DropsInfo(:,7),DropsInfo(:,8)];

end

function [GroundDrift] = GetDropGroundDrift(RDs,RadarParams,DropsSum)

% RadarParams = [479500,330500,1000,1000,100,100];
DY = RadarParams(3);
DX = RadarParams(4);
Nrows = RadarParams(5);
Ncols = RadarParams(6);

StartX=RadarParams(2);
StartY=RadarParams(1)-DY*(Nrows-1);
EndX = RadarParams(2)+DX*(Ncols-1);
EndY = RadarParams(1);

RD_Num=length(RDs);
[TotalDrops,~]=size(DropsSum);

Dom_Num = TotalDrops/RD_Num;
GroundDrift = zeros(Dom_Num,RD_Num);

for i=1:RD_Num
    for j=1:Dom_Num
        current_drop = DropsSum((i-1)*Dom_Num+j,:);
        if current_drop(1)>0
            e_x= current_drop(2);
            e_y= current_drop(3);

            curCol=round((e_x-StartX)/DX);
            curRow=round((e_y-StartY)/DY);

            if curRow<1 || curRow>Nrows || curCol<1 || curCol>Ncols
                continue;
            end
            
            Domain_id = (curCol - 1) * Nrows + curRow;

            GroundDrift(j,i) = Domain_id;
        end
    end
end


end

function [GroundEvp] = GetDropGroundEvp(RDs,DropsSum)

RD_Num=length(RDs);
[TotalDrops,~]=size(DropsSum);

Dom_Num = TotalDrops/RD_Num;
GroundEvp = zeros(Dom_Num,RD_Num);

for i=1:RD_Num
    for j=1:Dom_Num
        GroundEvp(j,i) = DropsSum((i-1)*Dom_Num+j,1);
    end
end

end

function [InitDrops_All] = ComputeDropsizeInfo(GridNets,RDs)

RD_Num=length(RDs);

InitDrops_All = [];
for i=1:RD_Num
    InitDrops(:,2:4)=GridNets;
    InitDrops(:,1)=RDs(i);
    InitDrops_All=[InitDrops_All;InitDrops];
end

end

function [dInfo] = DropEvaporationTick(dInfo,dAtmos, dt,Dv_dp_Table)

global IDrop_S;
global IAtmo_T IAtmo_P IAtmo_RH;

Pa=1.2;
Pwv=0.554;
Pw = 1000;

dp = Pa - Pwv;

D = dInfo(IDrop_S);
% unit: mm
[V] = ComputeDropVelocity(D);
% unit: m/s

dz = V * dt;

u = ComputeKinViscosity(dAtmos(IAtmo_T));  
% unit: m2/s
% refer to https://www.engineeringtoolbox.com/air-absolute-kinematic-viscosity-d_601.html

Dv = ComputeWaterVaporDiff(dAtmos(IAtmo_T), dAtmos(IAtmo_P)); 
% unit: m2/s
% refer to https://www.researchgate.net/post/Binary_diffusion_coefficients_for_water_vapour_in_air_at_normal_pressure
fv = 0.78 + 0.308 * (u/Dv)^0.333 * ((V * D * 10^(-3) / u)^0.5);

[Dv_dp] = ComputeWaterVaporTerm(dAtmos(IAtmo_T), dAtmos(IAtmo_RH), Dv_dp_Table);
% unit: g/cm s

dm = (4*pi * (D/2) * 10^(-1) * fv * Dv_dp) * dt * 10^(-3);


% unit: kg
% refer to paper: The Evaporation, Temperature and Thermal Relaxation

% igore dm2 
% dm2 = (4*pi * (D/2) * 10^(-3)* fv * Dv * dp) * dt;
% unit: kg

% both dD and D_n work, has been checked
dD = 2 * dm * dz / (pi * Pw * (D * 10^(-3))^2 * dt * V) * 10^3;
% refer to Page 73 of A study of raindrop size distribution and their variation with height
dD=(1/3) * dD; % by QQY on 20190814
D_n = D - dD;
% D_n2 = 2 * ( (D/2 * 10^(-3))^3 - dm / (1.333 * pi * Pw ))^0.333 * 10^3;

if D_n>0
    dInfo(IDrop_S) = D_n;
else
    dInfo(IDrop_S) = 0;
end

end

function [dInfo] = DropMoveTick(dInfo,dAtmos, dt)

global IDrop_S IDrop_DX IDrop_DY IDrop_DZ IDrop_AX IDrop_AY;
global IAtmo_VX IAtmo_VY;
D = dInfo(IDrop_S);
[V_drop] = ComputeDropVelocity(D);

dis_x=dAtmos(IAtmo_VX)*dt;
dis_y=dAtmos(IAtmo_VY)*dt;
dis_z=V_drop*dt;

dInfo(IDrop_DX)=dInfo(IDrop_DX)+dis_x;
dInfo(IDrop_DY)=dInfo(IDrop_DY)+dis_y;
dInfo(IDrop_DZ)=dInfo(IDrop_DZ)-dis_z;
dInfo(IDrop_AX)=dInfo(IDrop_AX)+dis_x;
dInfo(IDrop_AY)=dInfo(IDrop_AY)+dis_y;

end

function [dInfo] = DropMoveTick_NoDrift(dInfo, dt)

global IDrop_S IDrop_DZ;
D = dInfo(IDrop_S);
[V_drop] = ComputeDropVelocity(D);

dis_z=V_drop*dt;

dInfo(IDrop_DZ)=dInfo(IDrop_DZ)-dis_z;

end

function [dInfo] = SetDropDomain(dInfo,LayersHeight, dom_size, dom_x_num)

global IDrop_DI IDrop_LI IDrop_DX IDrop_DY IDrop_DZ;

domain_id = dInfo(IDrop_DI);
layer_id = dInfo(IDrop_LI);
dx=dInfo(IDrop_DX);
dy=dInfo(IDrop_DY);
dz=dInfo(IDrop_DZ);
semi_size=dom_size/2;

while dz<0 && layer_id>0
    layer_id=layer_id-1;
    if layer_id>0
        dz=dz+LayersHeight(layer_id);
    else
        dz=0;
    end
end

if dx>semi_size
    domain_id=domain_id+1;
    dx=(dx-semi_size)-semi_size;
elseif dx< -1*semi_size
    domain_id=domain_id-1;
    dx=(dx+semi_size)+semi_size;
end

if dy>semi_size
    domain_id=domain_id+dom_x_num;
    dy=(dy-semi_size)-semi_size;
elseif dy< -1*semi_size
    domain_id=domain_id-dom_x_num;
    dy=(dy+semi_size)+semi_size;
end

dInfo(IDrop_DI)=domain_id;
dInfo(IDrop_LI)=layer_id;
dInfo(IDrop_DX)=dx;
dInfo(IDrop_DY)=dy;
dInfo(IDrop_DZ)=dz;

end

function [dAtmos] = GetDropAtmos(dInfo,AtmosData,dom_xy_num)

global IDrop_DI IDrop_LI 
domain_id = dInfo(IDrop_DI);
layer_id = dInfo(IDrop_LI);
atmos_id = (layer_id-1)*dom_xy_num+domain_id;
dAtmos=AtmosData(atmos_id,:);

end

function [Dv_dp] = ComputeWaterVaporTerm(T, RH, Dv_dp_Table)

TT = [5,15,25,35,100];
RHT = [15,25,35,45,55,65,75,85,95,100];

T_ind=length(TT);
for ti=1:length(TT)
    if T<TT(ti)
        T_ind=ti;
        break;
    end
end

RH_ind=length(RHT);
for ti=1:length(RHT)
    if RH<RHT(ti)
        RH_ind=ti;
        break;
    end
end

Dv_dp = Dv_dp_Table(RH_ind, T_ind);

Dv_dp = Dv_dp*10^(-6);
% unit g/ cm*s

Dv_dp = Dv_dp*10^(-1);
% unit kg/ m*s

end

function [Dv] = ComputeWaterVaporDiff(temp, P)

T = temp + 273;
P0 = 1;
T0 = 256;

Dv = 1.97 * 10^(-5) * (P0 / P) * (T / T0)^1.685;

end

function [u] = ComputeKinViscosity(temp)

if temp<0
    u=12.85;
elseif temp<5
    u=13.28;
elseif temp<10
    u=13.72;
elseif temp<15
    u=14.16;
elseif temp<20
    u=14.61;
elseif temp<25
    u=15.06;
else
    u=15.52;
end

u=u*10^(-6);

end

function [DropsInfo_N] = RemoveEvapDrops(DropsInfo)

global IDrop_S 
[Drops_num, ~] = size(DropsInfo);
DropsInfo_N = [];

acc_n=0;
for i=1:Drops_num
    if DropsInfo(i,IDrop_S)>0 
        acc_n=acc_n+1;
        DropsInfo_N(acc_n,:) = DropsInfo(i,:);
    end
end

end

function [DropsInfo_Init] = DropsInitiation(DropsInitData,LayersHeight,DomainNet)

global IDrop_S 
[Drops_num, ~] = size(DropsInitData);

DropsInfo_Init=zeros(Drops_num,8);

for i=1:Drops_num
    if DropsInitData(i,IDrop_S)>0
        [dz,layer_id] = GetGivenPointLayer(LayersHeight,DropsInitData(i,4));
        [dx,dy,domain_id] = GetGivenPointDomain(DomainNet,DropsInitData(i,2),DropsInitData(i,3));
        DropsInfo_Init(i,:)=[DropsInitData(i,1),domain_id,layer_id,dx,dy,dz,DropsInitData(i,2),DropsInitData(i,3)];
    end
    % with fix order
end

fprintf('Drops Initiation is computed...\n');

end

function [V_drop] = ComputeDropVelocity(D)

V_drop=3.78 * D^0.67;

end

function [dx,dy,domain_id] = GetGivenPointDomain(DomainNet,givenX,givenY)

[dom_num, ~]=size(DomainNet);

min_dis=10^8;
for i=1:dom_num
    dis = sqrt((givenX-DomainNet(i,1))^2+(givenY-DomainNet(i,2))^2);
    if dis<= min_dis
        min_dis=dis;
        domain_id=i;
    end
end

dx=givenX-DomainNet(domain_id,1);
dy=givenY-DomainNet(domain_id,2);

end

function [acc_n] = GetEventDurationByFile(inputPath,AtmosPre)

NetList = dir(inputPath);
[m, ~] = size(NetList);

acc_n=0;
for i= 3 : m
    file_name=NetList(i, 1).name;
    data_type_loc=strfind(file_name,AtmosPre);
    if data_type_loc > 0
        acc_n=acc_n+1;
    end
end

end

function [dH,LayerIndex] = GetGivenPointLayer(LayersHeight,givenZ)

LayerNumber=length(LayersHeight);

dH=0;
LayerIndex=1;
acc_height=0;
for i=1:LayerNumber
    acc_height=acc_height+LayersHeight(i);
    if givenZ<= acc_height
        topH=acc_height-givenZ;
        dH=LayersHeight(i)-topH;
        LayerIndex=i;
        break;
    end
end

end

function [StormInfo] = GetStormInfo(StormInfoPath)

fid=fopen(StormInfoPath,'r');
str=fgets(fid);

StormInfo=[];
while str>0
    [DateCell]=textscan(str,'%d-%d-%d:%d:%d:%d %d-%d-%d:%d:%d:%d');
    StormInfo=[StormInfo;cell2mat(DateCell)];
    str=fgets(fid);
end

StormInfo=double(StormInfo);
fclose(fid);

end

function [] = WriteGridData(path,f_name,resultData,type)

new_path=strcat(path,f_name);
% dlmwrite(new_path,resultData,'-append','delimiter',' ');
if type==0 % integer
    dlmwrite(new_path,resultData,'delimiter','\t','precision','%.6g');
elseif type==1  % float with 2 precision
    dlmwrite(new_path,resultData,'delimiter','\t','precision','%.2f');
elseif type==2  % float with 4 precision
    dlmwrite(new_path,resultData,'delimiter','\t','precision','%.1f');
elseif type==3  % float with 4 precision
    dlmwrite(new_path,resultData,'delimiter','\t','precision','%g');
elseif type==4  % float with 4 precision
    dlmwrite(new_path,resultData,'delimiter','\t');
end

end
