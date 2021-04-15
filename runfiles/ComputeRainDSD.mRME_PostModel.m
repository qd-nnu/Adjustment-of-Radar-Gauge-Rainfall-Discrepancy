function [] = RME_PostModel(path,events)
% STEP 4 of RME model
% Process the result of the raindrop microphysical simulation

mode=3;   % 1-consider drift, 2-consider evp, 3-consider both
RDs = 0.1:0.1:6;
RadarParams = [479500,330500,1000,1000,100,100];

path=strcat(path,'hameldon_hill\');
RadarPath = strcat(path,'Radar Grid\');
MOutputPath=strcat(path,'Model_Output\');

StormInfo=GetStormInfo(strcat(path,events));
% 
Root_postRME(RadarPath,MOutputPath,StormInfo,1,RDs,RadarParams);
Root_postRME(RadarPath,MOutputPath,StormInfo,2,RDs,RadarParams);
Root_postRME(RadarPath,MOutputPath,StormInfo,3,RDs,RadarParams);
% Root_postRME_Hourly(MOutputPath,StormInfo,1,RDs,RadarParams);

end

function [] = Root_postRME_Hourly(MOutputPath,StormInfo,mode,RDs,RadarParams)

RunPName = char('Drift','Evp','DriftEvp');
SavePName = char('HRadar_D','HRadar_E','HRadar_DE');
HRadarPre='HRadar_O';  
DSDPre='HDSD';
AtmosPre='Atmos';
RadarFilePre='Radar';
DSDFilePre1='D0';
DSDFilePre2='dBNw';
DSDFilePre3='u';

[TotalEvents,~]=size(StormInfo);

for e=1:TotalEvents 
    eventfolder=sprintf('%.4d-%.2d-%.2d',StormInfo(e,1:3));  
    EPath = strcat(MOutputPath,eventfolder,'/'); 
    RadarEPath = strcat(EPath,HRadarPre,'/'); 
    ModelEPath = strcat(EPath,RunPName(mode,:),'/'); 
    AtmosEPath = strcat(EPath,AtmosPre,'/');
    DSDEPath = strcat(EPath,DSDPre,'/');
    OutputEPath = strcat(EPath,SavePName(mode,:),'/'); 
    
    if exist(OutputEPath, 'dir')
        rmdir(OutputEPath, 's');
    end
    mkdir(OutputEPath);
    
    [E_dur] = GetEventDurationByFile(RadarEPath,RadarFilePre);
    
    for t=1:E_dur
        if mode == 3
            ModelData_D = load(strcat(ModelEPath,RunPName(1,:),num2str(t)));
            ModelData_E = load(strcat(ModelEPath,RunPName(2,:),num2str(t)));
        else
        	ModelData = load(strcat(ModelEPath,RunPName(mode,:),num2str(t)));
        end
        
        RadarData = load(strcat(RadarEPath,RadarFilePre,num2str(t)));
        DSDData_D0 = load(strcat(DSDEPath,DSDFilePre1,num2str(t)));
        DSDData_dBNw = load(strcat(DSDEPath,DSDFilePre2,num2str(t)));
        DSDData_u = load(strcat(DSDEPath,DSDFilePre3,num2str(t)));

        [DSD] = ComputeDisND(DSDData_dBNw, DSDData_D0, DSDData_u, RDs);
        [DSD] = NDsInterp(DSD);

        switch mode
            case 1
                [N_RadarData] = GetDriftRadar(ModelData,DSD,RadarData,RDs);
            case 2
                [N_RadarData] = GetEvpRadar(ModelData,DSD,RadarData,RDs);
            case 3
                [EvpRadarData] = GetEvpRadar(ModelData_E,DSD,RadarData,RDs);
                [N_RadarData] = GetDriftRadar(ModelData_D,DSD,EvpRadarData,RDs);
        end
        
        WriteGridData(OutputEPath,strcat(SavePName(mode,:),num2str(t)),N_RadarData,1);
        fprintf('Event %s Time %d is computed...\n',eventfolder,t);

    end
    
end

fclose all;

end

function [] = Root_postRME(RadarPath,MOutputPath,StormInfo,mode,RDs,RadarParams)

SubLim =[0.5 699.5 0.5 1049.5]*1000;  % X_start,X_end,Y_start,Y_end
GridSize=1000;
maxFilePerHour=60;
tstepmin=1;
tstephour=1;

RunPName = char('Drift','Evp','DriftEvp');
SavePName = char('HRadar_D','HRadar_E','HRadar_DE');
DSDPre='DSD';
AtmosPre='Atmos';
RadarFilePre='Ham_';
RadarFilePost='_1km.dat';
DSDFilePre = 'DSD_';
DSDFilePost1='_D0';
DSDFilePost2='_dBNw';
DSDFilePost3='_u';

[TotalEvents,~]=size(StormInfo);

for e=1:TotalEvents 
    ta=datenum([StormInfo(e,1:4),0,0]); % starting date
    tb=datenum([StormInfo(e,7:10),0,0]); % ending date
    tb = tb + 1/24;
    
    eventfolder=sprintf('%.4d-%.2d-%.2d',StormInfo(e,1:3));  
    EPath = strcat(MOutputPath,eventfolder,'/'); 
    RadarEPath = strcat(RadarPath,eventfolder,'/'); 
    ModelEPath = strcat(EPath,RunPName(mode,:),'/'); 
    AtmosEPath = strcat(EPath,AtmosPre,'/');
    DSDEPath = strcat(EPath,DSDPre,'/');
    OutputEPath = strcat(EPath,SavePName(mode,:),'/'); 
    
    if exist(OutputEPath, 'dir')
        rmdir(OutputEPath, 's');
    end
    mkdir(OutputEPath);
    
    [E_dur] = GetEventDurationByFile(AtmosEPath,AtmosPre);
    
    t_0 = ta;
    for t=1:E_dur
        tm_0=t_0;
        if mode == 3
            ModelData_D = load(strcat(ModelEPath,RunPName(1,:),num2str(t)));
            ModelData_E = load(strcat(ModelEPath,RunPName(2,:),num2str(t)));
        else
        	ModelData = load(strcat(ModelEPath,RunPName(mode,:),num2str(t)));
        end

        RadarData=[];
        DSDData_D0=[];
        for i=1:maxFilePerHour
            if tm_0>=ta && tm_0<=tb
                cm_0 = datevec(tm_0);
                str_date=sprintf('%.4d%.2d%.2d%.2d%.2d',cm_0(1:5));  
                RaPath=strcat(RadarEPath,'Radar/',RadarFilePre,str_date,RadarFilePost);
                DSDPath=strcat(DSDEPath,DSDFilePre,str_date,DSDFilePost1);
                
                if exist(RaPath, 'file')
                    RadarData = load(RaPath);
                end
                
                if exist(DSDPath, 'file')
                    DSDPath2=strcat(DSDEPath,DSDFilePre,str_date,DSDFilePost2);
                    DSDPath3=strcat(DSDEPath,DSDFilePre,str_date,DSDFilePost3);
                    
                    DSDData_D0 = load(DSDPath);
                    DSDData_dBNw = load(DSDPath2);
                    DSDData_u = load(DSDPath3);
                end
                
                if ~isempty(RadarData) && ~isempty(DSDData_D0)
                    [DSD] = ComputeDisND(DSDData_dBNw, DSDData_D0, DSDData_u, RDs);
                    [DSD] = NDsInterp(DSD);
%                     [DSD] = NDsInterp(DSD); % second interp 
%                     [DSD] = NDsInterp(DSD); % third interp
                    switch mode
                        case 1
                            [N_RadarData] = GetDriftRadar(ModelData,DSD,RadarData,RDs);
                        case 2
                            [N_RadarData] = GetEvpRadar(ModelData,DSD,RadarData,RDs);
                        case 3
                            [EvpRadarData] = GetEvpRadar(ModelData_E,DSD,RadarData,RDs);
                            [N_RadarData] = GetDriftRadar(ModelData_D,DSD,EvpRadarData,RDs);
                    end
        
%                     ComparisonRadarImage(RadarData,N_RadarData,RadarParams);

                    WriteGridData(OutputEPath,strcat(SavePName(mode,:),'_',str_date),N_RadarData,1);
                    fprintf('Event %s Time %s is computed...\n',eventfolder,str_date);
                    
                    RadarData=[];
                    DSDData_D0=[];
                end
            end
            tm_0 = tm_0 + tstepmin/60/24;
            
        end 
        
        t_0 = t_0 + tstephour/24;
    end
    
end

fclose all;

end

function [EvpRadarData_2D] = GetDriftEvpRadar(ModelData_D,ModelData_E,DSD,RadarData,RDs)

[nrows,ncols] = size(RadarData);
[dom_num,RD_num]=size(ModelData_D);

Radar_One=RadarData(:);

ori_sum=zeros(dom_num,1);
evp_sum=zeros(dom_num,1);
for d=2:RD_num
    ND_2D=DSD(:,:,d);
    ND_1D = ND_2D(:);
    
    D_ori = RDs(d);
    vd = 3.78 * D_ori^0.67;
    dD = RDs(d)-RDs(d-1);

    for i=1:dom_num
        D_evp = ModelData(i,d);
        
        ori_sum(i)=ori_sum(i) + D_ori^3 * vd * ND_1D(i) * dD;
        evp_sum(i)=evp_sum(i) + D_evp^3 * vd * ND_1D(i) * dD;
    end
end

retain_ratio=zeros(dom_num,1);
for i=1:dom_num
    if ori_sum(i)>0
        retain_ratio(i) = evp_sum(i)/ ori_sum(i);
    end
end

EvpRadarData = Radar_One .* retain_ratio;

EvpRadarData_2D=reshape(EvpRadarData,[nrows,ncols]);

end

function [DriftRadarData_2D] = GetDriftRadar(ModelData,DSD,RadarData,RDs)

[nrows,ncols] = size(RadarData);
[dom_num,RD_num]=size(ModelData);

Radar_One=RadarData(:);
N_RadarData=zeros(dom_num,1);
Weig_sum = zeros(dom_num,1);

for d=1:RD_num
    D = RDs(d);
    DSD_2D=DSD(:,:,d);
    DSD_One = DSD_2D(:);

    for i=1:dom_num
        s_radar=Radar_One(i);
        s_dom= ModelData(i,d);
        if s_dom>0
            weig = DSD_One(s_dom)*D^3;
            Weig_sum(s_dom)=Weig_sum(s_dom)+weig;
            N_RadarData(s_dom)=N_RadarData(s_dom)+weig*s_radar;
        end
    end
end

[DriftRadarData] = ComputePositiveMean(Weig_sum, N_RadarData);
DriftRadarData_2D=reshape(DriftRadarData,[nrows,ncols]);

end

function [EvpRadarData_2D] = GetEvpRadar(ModelData,DSD,RadarData,RDs)

[nrows,ncols] = size(RadarData);
[dom_num,RD_num]=size(ModelData);

Radar_One=RadarData(:);

ori_sum=zeros(dom_num,1);
evp_sum=zeros(dom_num,1);
for d=2:RD_num
    ND_2D=DSD(:,:,d);
    ND_1D = ND_2D(:);
    
    D_ori = RDs(d);
    vd = 3.78 * D_ori^0.67;
    dD = RDs(d)-RDs(d-1);

    for i=1:dom_num
        D_evp = ModelData(i,d);
        
        ori_sum(i)=ori_sum(i) + D_ori^3 * vd * ND_1D(i) * dD;
        evp_sum(i)=evp_sum(i) + D_evp^3 * vd * ND_1D(i) * dD;
    end
end

retain_ratio=zeros(dom_num,1);
for i=1:dom_num
    if ori_sum(i)>0
        retain_ratio(i) = evp_sum(i)/ ori_sum(i);
    end
end

EvpRadarData = Radar_One .* retain_ratio;

EvpRadarData_2D=reshape(EvpRadarData,[nrows,ncols]);

end

function [NDs_N] = NDsInterp(NDs)

[nrows,ncols,d_num]=size(NDs);

NDs_N=NDs;
for d=1:d_num
    NDs_2D=NDs(:,:,d);
    for i=2:nrows-1
        for j=2:ncols-1
            if NDs_2D(i,j)<=0
%                 NDs_N(i,j,d)=max([NDs_2D(i-1,j),NDs_2D(i+1,j),NDs_2D(i,j-1),NDs_2D(i,j+1)]);
                NDs_N(i,j,d)=FindNearestPositive(NDs_2D,i,j);
            end
        end
    end
end

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

function [NDs] = ComputeDisND(dbNws, D0s, us, RDs)

[nrows,ncols]=size(dbNws);
d_num=length(RDs);
NDs=zeros(nrows,ncols,d_num);

for d=1:d_num
    D = RDs(d);
    for i=1:nrows
        for j=1:ncols
            if dbNws(i,j)>0 && D0s(i,j)>0 && us(i,j)>-3 && us(i,j)<40
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

function [Radar_M] = ComputePositiveMean(Weig_sum, RadarData)

[nl]=length(Weig_sum);
Radar_M=zeros(nl,1);

for i=1:nl
    if Weig_sum(i)>0
        Radar_M(i)=RadarData(i)/Weig_sum(i);
    end
end

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






