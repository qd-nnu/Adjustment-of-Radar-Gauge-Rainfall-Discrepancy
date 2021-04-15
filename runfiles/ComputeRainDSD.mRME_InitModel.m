  function [] = RME_InitModel(path,events)
% STEP 2 of RME model
% Get the event basic info to run the model
% including (1) Compute domain horizontal coordinates and raindrop size for each domain
%           (2) Atmosphere info (wind, temperature, pressure and RH)
%           (3) Compute raindrop size distribution
% All the above info are obtained for each time step

% Modified in April,1,2021, by yaru zhang

RadarLoc = [381060,428740,398];  
RadarParams = [479500,330500,1000,1000,100,10  0];

path=strcat(path,'hameldon_hill\');
% WRFOutPath=strcat('E:\WRF_Out\hameldon_hill\');
PolarPath = strcat(path,'Radar Polar\');
OutPath= strcat(path,'Model_Output\');

StormInfo=GetStormInfo(strcat(path,events));

BuildGrid_Net(RadarLoc,RadarParams,OutPath);

InitAtomsInfo(WRFOutPath,OutPath);

InitDualRadarDSD(StormInfo,PolarPath,OutPath);

DSDFilePost = char('D0','dBNw','u');
GetHourlyDualDSD(StormInfo,OutPath,'D0',-1);
GetHourlyDualDSD(StormInfo,OutPath,'dBNw',-1);
GetHourlyDualDSD(StormInfo,OutPath,'u',-999);

end

function [] = InitAtomsInfo(WRFOutPath,OutPath)

NetList = dir(WRFOutPath);
[m, ~] = size(NetList);
    
for i= 3 : m
       
    file_name=NetList(i, 1).name;
    data_type_loc=strfind(file_name,'.nc');
    if data_type_loc > 0
        date_loc=strfind(file_name,'_20');
        date_str=file_name(date_loc+1:date_loc+10);
        eventfolder=sprintf('%.4s-%.2s-%.2s',date_str(1:4),date_str(6:7),date_str(9:10));

        SaveAtmosPath=strcat(OutPath,eventfolder,'/Atmos/');      
        NCPath=strcat(WRFOutPath,file_name);

        if ~exist(SaveAtmosPath)
            mkdir(SaveAtmosPath);
        end

        GetAtmosData(NCPath,SaveAtmosPath);

        fprintf('%s is computed.\n',date_str);
    end
end

fclose all;

end

function [] = BuildGrid_Net(RadarLoc,RadarParams,OutPath)

DY = RadarParams(3);
DX = RadarParams(4);
NY = RadarParams(5);
NX = RadarParams(6);

XXs=RadarParams(2):DX:RadarParams(2)+DX*(NX-1);
YYs=RadarParams(1)-DY*(NY-1):DY:RadarParams(1);
[MXs,MYs]=meshgrid(XXs,YYs);

MXs_1D=MXs(:);
MYs_1D=MYs(:);

GridNets_1D=[MXs_1D,MYs_1D];
[GridNets_All] = GetBeamSpatialCoord(GridNets_1D,RadarLoc);

fprintf('Netork is generated.\n');
WriteGridData(OutPath,strcat('GridNets.txt'),GridNets_All,3);

end

function [] = InitDualRadarDSD(StormInfo,RadarPath,outputpath)

[TotalEvents,~]=size(StormInfo);

for e=1:TotalEvents
    event_date=sprintf('%.4d%.2d%.2d%.2d',StormInfo(e,1:4)); 
    eventfolder=sprintf('%.4d-%.2d-%.2d',StormInfo(e,1:3));  
    
    eventInputPath = strcat(RadarPath,event_date,'/Radar/'); 
    eventOutputPath = strcat(outputpath,eventfolder,'/DSD/'); 
    if exist(eventOutputPath, 'dir')
        rmdir(eventOutputPath, 's');
    end
    
    mkdir(eventOutputPath);

    NetList = dir(eventInputPath);
    [m, ~] = size(NetList);
    
    ZH=[];
    ZDR=[];
    for i= 3 : m
        file_name=NetList(i, 1).name;
        data_type1=strfind(file_name,'ZH');
        data_type2=strfind(file_name,'ZDR');
        data_type3=strfind(file_name,'KDP');
        date_loc=strfind(file_name,'_20');
        date_str=file_name(date_loc+1:date_loc+12);
        
        if data_type1 > 0
            [ZH] = load(strcat(eventInputPath,file_name));
        end
        
        if data_type2 > 0
            [ZDR] = load(strcat(eventInputPath,file_name));
        end
        
        if data_type3 > 0
            [KDP] = load(strcat(eventInputPath,file_name));
        end
        
        if ~isempty(ZH) && ~isempty(ZDR)
            [ZH_G] = CovertPolarToGrid(ZH);
            [ZDR_G] = CovertPolarToGrid(ZDR);
            [KDP_G] = CovertPolarToGrid(KDP);
            
            [R] = ComputeRainfall(ZH_G,ZDR_G);
            [Nw, D0, u] = DSDParamRetrieval(ZDR_G, ZH_G, KDP_G);

            WriteGridData(eventOutputPath,strcat('DSD_',date_str,'_dBNw'),Nw,1);
            WriteGridData(eventOutputPath,strcat('DSD_',date_str,'_D0'),D0,1);
            WriteGridData(eventOutputPath,strcat('DSD_',date_str,'_u'),u,1);
%             WriteGridData(eventOutputPath,strcat('HamRain_',date_str,'_R'),R,1);

            fprintf('Event%s-Date%s is computed.\n',event_date,date_str);

            ZH=[];
            ZDR=[];
            KDP=[];
        end    
    end
end
end

function [] = GetHourlyDualDSD(StormInfo,OutPath,DSDFilePost,tol)

[TotalEvents,~]=size(StormInfo);
GN=100;

for e=1:TotalEvents     
    outputfolder=sprintf('%.4d-%.2d-%.2d',StormInfo(e,1:3));  
    
    inputPath = strcat(OutPath,outputfolder,'/DSD/'); 
    outputPath = strcat(OutPath,outputfolder,'/HDSD/');
    if ~exist(outputPath, 'dir')
        mkdir(outputPath);
    end
    
    NetList = dir(inputPath);
    [m, ~] = size(NetList);
    
    file_name=NetList(3, 1).name;
    date_loc=strfind(file_name,'_20');
    last_str=file_name(date_loc+1:date_loc+10);
    
    DSD_All=zeros(GN,GN,1);
    acc_n=0;
    acc_h=0;
    for i= 3 : m
        file_name=NetList(i, 1).name;
        data_type_loc=strfind(file_name,DSDFilePost);
        if data_type_loc > 0
            date_loc=strfind(file_name,'_20');
            date_str=file_name(date_loc+1:date_loc+10);
            [DSD] = load(strcat(inputPath,file_name));
            
            if ~strcmp(last_str,date_str)
                [DSD_All] = NDsTemporalInterp(DSD_All,tol);
                acc_n=acc_n+1;
                WriteGridData(outputPath,strcat(DSDFilePost,num2str(acc_n)),DSD_All,1);
                
                DSD_All=zeros(GN,GN,1);
                last_str=date_str;
                acc_h=0;
            end
            
            acc_h=acc_h+1;
            DSD_All(:,:,acc_h)=DSD;
      
        end
    end
    
    acc_n=acc_n+1;
    [DSD_All] = NDsTemporalInterp(DSD_All,tol);
    WriteGridData(outputPath,strcat(DSDFilePost,num2str(acc_n)),DSD_All,1);
    fprintf('Event%s is computed.\n',outputfolder);
                
    fclose all;
end

end

function [DSD_Interp] = NDsTemporalInterp(DSD_All,tol)

[nrows,ncols,~]=size(DSD_All);

NDs_2D=zeros(nrows,ncols);
for i=1:nrows
    for j=1:ncols
        cur_DSD = DSD_All(i,j,:);
        cur_DSD_invalid = cur_DSD(cur_DSD>tol);
        if ~isempty(cur_DSD_invalid)
            NDs_2D(i,j)= median(cur_DSD_invalid);
        else
            NDs_2D(i,j)= -999;
        end
    end
end

DSD_Interp=NDs_2D;
% for i=1:nrows
%     for j=1:ncols
%         if NDs_2D(i,j)<=0
%             DSD_Interp(i,j)=FindNearestPositive(NDs_2D,i,j);
%         end
%     end
% end

end

function [pvalue] = FindNearestPositive(Data,i,j)

[nrows,ncols]=size(Data);
max_radius = nrows/2;

pvalue=0;
for r=1:max_radius
    sub_left=j-r;
    sub_right=j+r;
    sub_top=i-r;
    sub_bottom=i+r;
    
    if sub_left<=0
        sub_left=1;
    end
    
    if sub_right>ncols
        sub_right=ncols;
    end
    
    if sub_top<=0
        sub_top=1;
    end
    
    if sub_bottom>nrows
        sub_bottom=nrows;
    end
    data_sub=Data(sub_top:sub_bottom,sub_left:sub_right);
    pvalue = max(max(data_sub));
    if pvalue>0
        break;
    end
end
end

function [Grid] = CovertPolarToGrid(Polar)

Tol=-5;
dom_s=50000;
res=1000;
bin_length=600;
bin_number=140;   % only the inner part is considered
theta=0:359;

Polar_inner = Polar(:,1:bin_number);

rr=1:bin_length:bin_length*bin_number;
rr=rr';

% for i = 1:length(theta)
%     theta(i) = theta(i) - 90;
%     if (theta(i)<0)
%         theta(i) = theta(i) + 359;
%     end
% end
theta=2*pi*theta/360;

P_XX = rr*cos(theta);
P_YY = rr*sin(theta);
P_ZZ = Polar_inner';

G_X=-1*dom_s:res:dom_s-res;
G_Y=-1*dom_s:res:dom_s-res;
[G_XX,G_YY] = meshgrid(G_X,G_Y);

[g_rows,g_cols]=size(G_XX);
[p_rows,p_cols]=size(P_XX);
R_I = zeros(g_rows,g_cols);
R_I_num = zeros(g_rows,g_cols);

for i=1:p_rows
    for j=1:p_cols
        p_x=P_XX(i,j);
        p_y=P_YY(i,j);
        
        g_row = round((p_y + dom_s)/res);
        g_col = round((p_x + dom_s)/res);
        
        if g_row<1 || g_col<1 || g_row>g_rows || g_col>g_cols
            continue;
        end
        
        if P_ZZ(i,j)>R_I(g_row,g_col)
            R_I(g_row,g_col)=P_ZZ(i,j);
            R_I_num(g_row,g_col)=R_I_num(g_row,g_col)+1;
        end

    end
end

Grid = R_I;
% [Grid] = ComputePositiveMean(R_I,R_I_num);

end

function [R] = ComputeRainfall(ZH,ZDR)

a=0.0142;   % 0.0142 0.0159 0.0144  0.0067  0.00746
b=0.77;     % 0.77   0.737  0.761   0.927  0.945
c=-1.67;    % -1.67 -1.03   -1.51   -3.43   -4.76

[nrows,ncols]=size(ZH);
R=zeros(nrows,ncols);

for i=1:nrows
    for j=1:ncols
        if ZH(i,j)>0 && ZDR(i,j)>0
            R(i,j)=a * ZH(i,j)^b * ZDR(i,j)^c;
        end
    end
end

end

function [Nw, D0, u] = DSDParamRetrieval(ZDR, ZH, KDP)

[nrows,ncols]=size(ZDR);

Nw=zeros(nrows,ncols);
D0=zeros(nrows,ncols);
u=zeros(nrows,ncols);

for i=1:nrows
    for j=1:ncols
        zdr = ZDR(i,j);
        zh = ZH(i,j);
        kdp = KDP(i,j);
        [paramDSD] = DSDParamRetrieval_3Para(zdr, zh, kdp);
        Nw(i,j)=paramDSD(1);
        D0(i,j)=paramDSD(2);
        u(i,j)=paramDSD(3);
    end
end

end

function [R_I] = ComputePositiveMean(R_I,R_I_num)

[nrows,ncols]=size(R_I);

for i=1:nrows
    for j=1:ncols
        if R_I_num(i,j)>1
            R_I(i,j)=R_I(i,j)/R_I_num(i,j);
        end
    end
end

end

function [paramDSD] = DSDParamRetrieval_3Para(zdr, zh, kdp)
% refer to the paper: Comparison of Polarimetric Radar Drop Size
% Distribution Retrieval Algorithms by Edward A. Branges

% the algorithm is named as 'constrained-gamma', which is considered to be
% better than the algorithm 'beta-method' proposed by Bringi et al.
% see: "Comparison of Two Raindrop Size Distribution Retrieval Algorithms" 

% the equation work best for 0.3<=zdr<=3 dB
if zh>0 && zdr>0 && zdr<2
    zh_r = 0.1 * 10^zh;
   % by QQY on 20190813

   D0 = 0.171* zdr^3 - 0.725 * zdr^2 + 1.479 * zdr + 0.717;
  u = 6.084 * D0^2 - 29.85 * D0 + 34.64;


%     D0=2.575*zdr^3-5.406 *zdr^2+4.826 *zdr+ 0.434;
%     u = 20.92 * D0^2 - 58.02* D0 + 41.85;
    if D0>2
        D0
    end

    % The Nw calculaed by Brandes is too large, instead, the method to
    % compute water content by Kim is used
    % see: "Retrieval of Three-Dimensional Raindrop Size Distribution Using
    % X-Band Polarimetric Radar Data"
    if kdp <= 0.3 || zh <= 35
        W = 0.00383 * zh^0.55;
    else
        W = 0.991 * kdp^0.713;% by QQY on 20190813
        %W=0.01583*kdp^0.6652;  % by QQY on 20190813
    end
    
%     W = 5.589 * 10^(-4) * zh_r * 10^(0.223 * zdr^2 - 1.124 * zdr);    
%     W = 10^(-2.9) * zh_r * 10^(1.72 * zdr^2 - 2.48 * zdr - 0.5 * zdr^3 + 0.06 * zdr^4);
    Nw = W / (1.73 * 10^(-5) * D0^4);
    
    Nw_dB = log10(Nw);
    
    paramDSD=[Nw_dB, D0, u];
else
    paramDSD=[-1, -1, -999];
    
end

end

function [paramDSD] = DSDParamRetrieval_NormGamma(zdr, zh)
% refer to the paper: Comparison of Polarimetric Radar Drop Size
% Distribution Retrieval Algorithms by Edward A. Branges

% the algorithm is named as 'constrained-gamma', which is considered to be
% better than the algorithm 'beta-method' proposed by Bringi et al.
% see: Comparison of Two Raindrop Size Distribution Retrieval Algorithms 
% for X-Band Dual Polarization Observations

% the equation work best for 0.3<=zdr<=3 dB
if zh>0 && zdr>=0.3 && zdr<=3
    zh_r = 0.1 * 10^zh;
    
    D0 = 0.171* zdr^3 - 0.725 * zdr^2 + 1.479 * zdr + 0.717;
    u = 6.084 * D0^2 - 29.85 * D0 + 34.64;

%     W = 5.589 * 10^(-4) * zh_r * 10^(0.223 * zdr^2 - 1.124 * zdr);
    W = 10^(-2.9) * zh_r * 10^(1.72 * zdr^2 - 2.48 * zdr - 0.5 * zdr^3 + 0.06 * zdr^4);
    Nw = W / (1.73 * 10^(-5) * D0^4);
    
    if Nw>10^10
        Nw=Nw;
    end
    
    Nw_dB = log10(Nw);
    paramDSD=[Nw_dB, D0, u];
else
    paramDSD=[-1, -1, -999];
    
end

end

function [BeamCoor] = GetBeamSpatialCoord(PixelData,RadInfo)

RadX=RadInfo(1);
RadY=RadInfo(2);
RadH=RadInfo(3);

[pixel_number,~]=size(PixelData);

BeamCoor=zeros(pixel_number,3);
pixel_DD=zeros(pixel_number,1);
dH=zeros(pixel_number,1);
for i=1 : pixel_number
   pixel_D=sqrt((PixelData(i,1)-RadX)*(PixelData(i,1)-RadX)+...
    (PixelData(i,2)-RadY)*(PixelData(i,2)-RadY)); 
   pixel_DD(i)=pixel_D/1000;
   
   dH(i)=ComputeBeamDHeight(pixel_D)/2;
   beamHeight = ComputeBeamHeight(pixel_D)+RadH;

   BeamCoor(i,1)=PixelData(i,1);
   BeamCoor(i,2)=PixelData(i,2);
   BeamCoor(i,3)=int32(beamHeight);
%    BeamCoor(i,4)=pixel_D/1000;
end

end

function [beamHeight] = ComputeBeamHeight(distance)

REA=0.5; % radar elevation angle (in degrees),
IR=1.21;  % refractive index (1.21), 
Re=6371;  %  radius of the earth (6371 km) 

if distance<1000
    D=distance; % slant range observed on radar in kilometers (km)
else
    D=distance/1000;
end

beamHeight=D*sind(REA)+D^2/(2*IR*Re);

beamHeight=beamHeight*1000;  % in m

end

function [DHeight] = ComputeBeamDHeight(distance)

REA=0.5; % radar elevation angle (in degrees),
AB=0.95; % the angular beamwidth
Pi=3.14;

if distance<1000
    D=distance; % slant range observed on radar in kilometers (km)
else
    D=distance/1000;
end

bw=(AB*Pi/180)*D;

DHeight=bw*cos(REA)*1000;  % in m

end


function [NDs] = ComputeDisND(dbNws, D0s, us, RDs)

[nrows,ncols]=size(dbNws);
d_num=length(RDs);
NDs=zeros(nrows,ncols,d_num);

for d=1:d_num
    D = RDs(d);
    for i=1:nrows
        for j=1:ncols
            if dbNws(i,j)>0 && D0s(i,j)>0 && us(i,j)>-5
                Nw = 10^(dbNws(i,j));
                [NDs(i,j,d)] = ComputeND(Nw, D0s(i,j), us(i,j), D);
            end
        end
    end
end

end

function [N_D] = ComputeND(Nw, D0, u, D)

fu = (6 * (3.67+u)^(u+4))/(3.67^4 * gamma(u+4));
N_D = Nw * fu * (D/D0)^u * exp(-1 * (u+3.67) * (D/D0));

end

function [PM] = ComputePositiveMean_1D(Data)

[nl,nc]=size(Data);

acc_n=zeros(1,nc);
acc_sum=zeros(1,nc);
for i=1:nl
    for j=1:nc
        if Data(i,j)>0
            acc_n(j)=acc_n(j)+1;
            acc_sum(j)=acc_sum(j)+Data(i,j);
        end
    end
end

PM=zeros(1,nc);
for j=1:nc
    if acc_n(j)>0
        PM(j)=acc_sum(j)/acc_n(j);
    end
end

end

function [] = GetAtmosData(NCPath,SaveAtmosPath)

ncid = netcdf.open(NCPath,'NC_WRITE');
varid = netcdf.inqVarID(ncid,'U');
U_Data=netcdf.getVar(ncid,varid);
varid = netcdf.inqVarID(ncid,'V');
V_Data=netcdf.getVar(ncid,varid);
varid = netcdf.inqVarID(ncid,'W');
W_Data=netcdf.getVar(ncid,varid);
varid = netcdf.inqVarID(ncid,'TK');
TK_Data=netcdf.getVar(ncid,varid);
varid = netcdf.inqVarID(ncid,'P');
P_Data=netcdf.getVar(ncid,varid);
varid = netcdf.inqVarID(ncid,'RH');
RH_Data=netcdf.getVar(ncid,varid);

[~,SN_Number,Height_Number,Time_Number]=size(U_Data);
[WE_Number,~,~,~]=size(V_Data);

for t=1:Time_Number
    Save_name=strcat('Atmos',num2str(t));
    FullSavePath=strcat(SaveAtmosPath,Save_name);

    if exist(FullSavePath)
        delete(FullSavePath);
    end
    
    Wind_Number=0;
    AtmosData=[];
    for h=1:Height_Number
       U_2D_Data(:,:)=U_Data(:,:,h,t);
       V_2D_Data(:,:)=V_Data(:,:,h,t);
       W_2D_Data(:,:)=W_Data(:,:,h,t);
       TK_2D_Data(:,:)=TK_Data(:,:,h,t);
       P_2D_Data(:,:)=P_Data(:,:,h,t);
       RH_2D_Data(:,:)=RH_Data(:,:,h,t);
       for i=1:SN_Number
          for j=1: WE_Number
              Wind_Number=Wind_Number+1;
              AtmosData(Wind_Number,1)=U_2D_Data(j,i);
              AtmosData(Wind_Number,2)=V_2D_Data(j,i);
              AtmosData(Wind_Number,3)=W_2D_Data(j,i);
              AtmosData(Wind_Number,4)=TK_2D_Data(j,i);
              AtmosData(Wind_Number,5)=P_2D_Data(j,i);
              AtmosData(Wind_Number,6)=RH_2D_Data(j,i);
          end
       end
       
       fprintf('%d/%d is computed.\n',t,h);
    end
    
    WriteGridData(SaveAtmosPath,Save_name,AtmosData,1);
end

end

function f = fm(a,x)
    
f=a(1)*x.^a(2);

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
    dlmwrite(new_path,resultData,'delimiter','\t','precision','%.4f');
elseif type==3  % float with 4 precision
    dlmwrite(new_path,resultData,'delimiter','\t','precision','%d');
end

end
