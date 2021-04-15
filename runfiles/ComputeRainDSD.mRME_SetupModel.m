function [] = RME_SetupModel(path)
% STEP 1 of RME model
% Get the perminent info to run the model, including domain_net and layers
% Input parameters:
% path: data storage path
% Output parameters:
% height 

% Modified in April 10,2021, by yaru zhang

path=strcat(path,'hameldon_hill\');
DataPath=strcat(path,'Model_Input\');
OutPath=strcat(path,'Model_Output\');
WRFOutPath=strcat(DataPath,'wrfout_2016-03_01-00_01-18.nc');

[NX, LayerHeight] =BuildDomain_Net(WRFOutPath);

WriteGridData(OutPath,'Domain_Net.txt',NX,3);
WriteGridData(OutPath,'Layers_Height.txt',LayerHeight,1);

end

function [NX, LayerHeight] = BuildDomain_Net(WRFOutPath)
% 


CEN_WE=381060;
CEN_SN=428740;
Grid_Size=5000;
Grid_Number_WE=17; % half per aspect
Grid_Number_SN=17;
        
ncid = netcdf.open(WRFOutPath,'NC_WRITE');

varid = netcdf.inqVarID(ncid,'PH');
PH_Data=netcdf.getVar(ncid,varid);
varid = netcdf.inqVarID(ncid,'PHB');
PHB_Data=netcdf.getVar(ncid,varid);
varid = netcdf.inqVarID(ncid,'HGT');
HGT_Data=netcdf.getVar(ncid,varid);

dx=Grid_Size*Grid_Number_WE;
dy=Grid_Size*Grid_Number_SN;
startX=CEN_WE-dx;
endX=CEN_WE+dx;
startY=CEN_SN-dy;
endY=CEN_SN+dy;
[Grid_X,Grid_Y]=meshgrid(startX:Grid_Size:endX,startY:Grid_Size:endY); 

CX=Grid_X';
CY=Grid_Y';
NewX=CX(:);
NewY=CY(:);
NX(:,1)=NewX(:);
NX(:,2)=NewY(:);

[WE_Number,SN_Number,Height_Number,Time_Number]=size(PH_Data);

MeanHeight=zeros(Height_Number,1);
for h=1:Height_Number
    for t=1:Time_Number
       for i=1:WE_Number
          for j=1: SN_Number
              PH=PH_Data(i,j,h,t);
              PHB=PHB_Data(i,j,h,t);
              HGT=HGT_Data(i,j,t);
              curHeight(t,i,j)=(PH+PHB)/9.81-HGT;
          end
       end
    end
    
    MeanHeight(h)=mean(mean(mean(curHeight)));
end

LayerHeight=zeros(Height_Number-1,1);
for h=1:Height_Number-1
    LayerHeight(h)=MeanHeight(h+1)-MeanHeight(h);
end


end

function [] = WriteGridData(path,f_name,resultData,type)

new_path=strcat(path,f_name);
% dlmwrite(new_path,resultData,'-append','delimiter',' ');
if type==0 % output type is reasonably selected by the size and type of data
    dlmwrite(new_path,resultData,'delimiter','\t','precision','%.6g');
elseif type==1  % float with 2 precision
    dlmwrite(new_path,resultData,'delimiter','\t','precision','%.2f');
elseif type==2  % float with 4 precision
    dlmwrite(new_path,resultData,'delimiter','\t','precision','%.4f');
elseif type==3  % integer
    dlmwrite(new_path,resultData,'delimiter','\t','precision','%d');
end

end
