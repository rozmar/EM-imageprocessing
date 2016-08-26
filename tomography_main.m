%% locating the files, ordering
close all
clear all
locations=marcicucca_locations;
cd([locations.tgtardir,'EMdata/Tomogram/'])
alltherawdirs=uipickfiles;
for dirnum=1:length(alltherawdirs)
    
%     dirs.raw=[locations.tgtardir,'EMdata/Tomogram/Human 0410122 AZ1 scale 100nm jpgs/'];
%     dirs.filtered=[locations.tgtardir,'EMdata/FilteredTomogram/Human 0410122 AZ1 scale 100nm jpgs/'];
    dirs.raw=[alltherawdirs{dirnum},'/'];
    disp(dirs.raw)
    temphely=strfind(dirs.raw,'/Tomogram/');
    dirs.filtered=[dirs.raw(1:temphely),'Filtered',dirs.raw(temphely+1:end)];
    files=dir(dirs.raw);
    files([files.isdir])=[];
    clear filedata
    for fnum=1:length(files)
        doth=strfind(files(fnum).name,'.');
        doth=doth(end);
        filedata(fnum).name=files(fnum).name;
        filedata(fnum).ext=files(fnum).name(doth+1:end);
        filedata(fnum).num=str2num(files(fnum).name(doth-3:doth-1));
    end
    filedata(~strcmp({filedata.ext},'jpg'))=[];
    [~,ix]=sort([filedata.num]);
    filedata=filedata(ix);
    %% loadin the files
    A=imread([dirs.raw,filedata(1).name]);
    dims=size(A);
    dims(3)=length(filedata);
    rawdata=zeros(dims,'uint8');
    tic
    for fnum=1:length(filedata)
        A=imread([dirs.raw,filedata(fnum).name]);
        rawdata(:,:,fnum)=A(:,:,1);
        
    end
    toc
    
    %% filtering
    
    pixelxsize=.000379;
    pixelysize=.000379;
    pixelzsize=.006216;
    
    filtersize=.001;
    filterxsize=filtersize;
    filterysize=filtersize;
    filterzsize=filtersize;
    
    filterstep=([filterxsize/pixelxsize,filterysize/pixelysize,filterzsize/pixelzsize]);
    
    h=fspecial3('gaussian',round(filterstep*2*2.354));
    % h=ones([2,2,round(filterstep(3))]);
    % h=h/round(filterstep(3)*4);
    filtereddata=zeros(size(rawdata),'uint8');
    tic
    filtereddata=imfilter(rawdata,h);
    toc
    
    %% saving
    a=dir(dirs.filtered);
    if isempty(a)
        mkdir(dirs.filtered)
    end
    for fnum=1:length(filedata)
        imwrite(filtereddata(:,:,fnum),[dirs.filtered,filedata(fnum).name],'jpg')
    end
end
