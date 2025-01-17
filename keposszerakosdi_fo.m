%% kepek beolvasasa
clear all
close all
plotthestuff=0;
locations=marcicucca_locations;
dirs.alapkonyvtar=[locations.EMdir,'EMdata/'];
blacklist={'SA-MAG_X10k_5A4','SA-MAG_X8000_7A9'};
% dirs.tifkonyvtar=[dirs.alapkonyvtar,'rawtif/Elmi_patosz/juci/0509164cs/bu2/'];
% dirs.montage2ddir=[dirs.alapkonyvtar,'2dmontage/Elmi_patosz/juci/0509164cs/bu2/'];
% dirs.tifkonyvtar=[dirs.alapkonyvtar,'rawtif/Elmi_patosz/juci/0403183j_AIS/'];
% dirs.montage2ddir=[dirs.alapkonyvtar,'2dmontage/Elmi_patosz/juci/0403183j_AIS/'];
% dirs.montage3ddir=[dirs.alapkonyvtar,'3dmontage/Elmi_patosz/juci/0403183j_AIS/'];
% dirs.maskdir=[dirs.alapkonyvtar,'2dmasks/Elmi_patosz/juci/0403183j_AIS/'];
% dirs.reconstructSERdir=[dirs.alapkonyvtar,'3dmontage/Elmi_patosz/juci/0403183j_AIS/ReconstructSeries/'];
dirs.tifkonyvtar=[dirs.alapkonyvtar,'rawtif/Elmi_patosz/juci/Rosehip szinapszishelyek/0506102mg_bu9/'];
dirs.montage2ddir=[dirs.alapkonyvtar,'2dmontage/Elmi_patosz/juci/Rosehip szinapszishelyek/0506102mg_bu9/'];
dirs.tifkonyvtar=[dirs.alapkonyvtar,'rawtif/Elmi_patosz/juci/Rosehip szinapszishelyek/1404163rm_bu8/'];
dirs.montage2ddir=[dirs.alapkonyvtar,'2dmontage/Elmi_patosz/juci/Rosehip szinapszishelyek/1404163rm_bu8/'];

dirstodo=uipickfiles('FilterSpec',dirs.alapkonyvtar);
for diri=1:length(dirstodo)
    dirs.tifkonyvtar=[dirstodo{diri},'/'];
    temp=strfind(dirs.tifkonyvtar,'rawtif');
    dirs.montage2ddir=[dirs.tifkonyvtar(1:temp-1),'2dmontage',dirs.tifkonyvtar(temp+6:end)];
    a=dir(dirs.montage2ddir);
    if isempty(a)
        mkdir(dirs.montage2ddir)
    end
%%
% alapkonyvtar='/home/rozmar/Mount/TGTAR_1/ANALYSISdata/marci/_elmi/';
filek=dir(dirs.tifkonyvtar);
filek([filek.isdir])=[];
todel=zeros(size(filek));
for fnum=1:length(filek) %only tif files are to remain
    if ~any(strfind(filek(fnum).name,'.tif')) | ~(strfind(filek(fnum).name,'.tif')+3==length(filek(fnum).name)) | any(strfind(filek(fnum).name,'hely'))
        todel(fnum)=1;
    end
end
filek(find(todel))=[];
if any(strfind([filek.name],'('))
    for fnum=1:length(filek)
        fname=filek(fnum).name;
        %     filek(fnum).nyerskep=imread([tifkonyvtar,fname]);
        if any(strfind(filek(fnum).name,'('))
            hyp1=strfind(fname,'(');
            hyp2=strfind(fname,')');
            filek(fnum).z=(fname(1:hyp1(end)-1));
            filek(fnum).idx=str2num(fname(hyp1(end)+1:hyp2(end)-1))+1;
        else
            dot=strfind(fname,'.');
            filek(fnum).z=(fname(1:dot(end)-1));
            filek(fnum).idx=1;
            
        end
        
        
    end
else
    for fnum=1:length(filek)
        fname=filek(fnum).name;
        %     filek(fnum).nyerskep=imread([tifkonyvtar,fname]);
        hyp=strfind(fname,'_');
        filek(fnum).z=(fname(1:hyp(end)-1));%hyp(end-1)+
        filek(fnum).idx=str2num(fname(hyp(end)+1:end-4));
        
    end
end

%clear blacklsited files
todel=[];
for i=1:length(blacklist)
    todel=[todel,find(strcmp({filek.z},blacklist{i}))];
end
filek(todel)=[];
keposszerakosdi_2D_main(filek,dirs,plotthestuff)
end
return

% %% játszós
% mthresh=500;
% matchthreshold=1;
% close all
% pic=imread('/home/rozmar/Mount/HDDS/Backup/Users/marci/juci2.jpg');
% pic = rgb2gray(pic);
% pic2=imrotate(pic(200:500,100:300),100);
%
%     ptsbase  = detectSURFFeatures(pic,'MetricThreshold',mthresh,'NumOctaves',2,'NumScaleLevels',15);
%     ptsnew = detectSURFFeatures(pic2,'MetricThreshold',mthresh,'NumOctaves',2,'NumScaleLevels',15);
%     [featuresBase,  validPtsBase]  = extractFeatures(pic,  ptsbase);
%     [featuresNew, validPtsNew] = extractFeatures(pic2, ptsnew);
%     indexPairs = matchFeatures(featuresBase, featuresNew,'MatchThreshold',matchthreshold);
%     matchedbase=validPtsBase(indexPairs(:,1));
%     matchednew=validPtsNew(indexPairs(:,2));
%     [tform, inliernew, inlierbase] = estimateGeometricTransform(matchednew, matchedbase, 'similarity');
%
% figure(1)
%     clf
%     subplot(2,1,1)
%     imagesc(pic)
%     hold on
%     plot(validPtsBase.selectStrongest(100),'showOrientation',true);
%     colormap gray
%     subplot(2,1,2)
%     imagesc(pic2)
%     hold on
%     plot(validPtsNew.selectStrongest(100),'showOrientation',true);
%     colormap gray
%     figure(2)
%     clf
%     showMatchedFeatures(pic,pic2,inlierbase,inliernew,'montage');
% %     tform.T(3,1)=tform.T(3,1);
% %     tform.T(3,2)=tform.T(3,2);
%   tform.T(3,1)=0;
%     tform.T(3,2)=0;
%     outpic=imwarp(pic2,tform);
%     figure(3)
%     clf
%     subplot(2,1,1)
%     imagesc(pic);
%     subplot(2,1,2)
%     imagesc(outpic)
%% imtransform based on reconstruct
%read files
jpgquality=75;
files=dir(dirs.reconstructSERdir);
files([files.isdir])=[];
for filei=1:length(files)
    disp(['filenum: ',num2str(filei)])
    fid=fopen([dirs.reconstructSERdir,files(filei).name]);
    A=textscan(fid,'%s');
    fclose(fid);
    idx=find(strcmp(A{1},'xcoef="'));
    if ~isempty(idx)
        xcoef=[str2num(A{1}{idx+2});str2num(A{1}{idx+3});str2num(A{1}{idx+1})];
        idx=find(strcmp(A{1},'ycoef="'));
        ycoef=[str2num(A{1}{idx+2});str2num(A{1}{idx+3});str2num(A{1}{idx+1})];
        T=[xcoef,ycoef,[0;0;1]];
        tform=affine2d(T);
        for i=1:length(A{1})
            if any(strfind(A{1}{i},'src="'))
                filename=A{1}{i}(6:end-1);
                disp(['filename: ',filename]);
            end
            
        end
        pic=imread([dirs.montage3ddir,filename(1:end-3),'tif']);
        newpic=imwarp(pic,tform);
        imwrite(newpic,[dirs.reconstructSERdir,filename],'jpg','Bitdepth',8,'Quality',jpgquality);
    end
end

%%
% T=Torig(:,[2,3,1])
% 
% 

%% mintavételezés változtatása és 3D rekonstrukcio Reconstruct programhoz
jpgquality=75;
overwrite=1;
szorzo=5;
sdvalup=3;
sdvaldown=3;
gauci=round([3,3]*szorzo);
pixdif=0; %200
downsize=round(5*szorzo);
pixdif=ceil(pixdif/downsize);
mthresh= 50;
matchthreshold=10;
maxdist=10;
switches.pixelarea=.923^2;
switches.thresholdstart=0;
switches.plotthestuff=1;
h=fspecial('Gaussian',gauci*6,gauci(1));
filek=dir(dirs.montage2ddir);
filek([filek.isdir])=[];
todel=[];
for fnum=1:length(filek)
    fname=filek(fnum).name;
    %     filek(fnum).nyerskep=imread([tifkonyvtar,fname]);
    hyp=strfind(fname,'_');
    filek(fnum).ID=(fname(1:hyp(end)-1));%hyp(end-1)+
    filek(fnum).idx=str2num(fname(hyp(end)+1:end-4));
    if ~any(strfind(filek(fnum).ID,'X15k'))
        todel=[todel,fnum];
    end
end
filek(todel)=[];
filek(1)=[]; %az elso file tul kicsi...

tforms=struct;
for roundi =1:1 %first round -rotation - second round - revision - 

    fnum=1;
    lastacceptedrotation=0;
    while fnum<length(filek)
        fnum=fnum+1;
        fname=filek(fnum).name;
        a=dir([dirs.montage3ddir,filek(fnum).name]);%(1:end-4),'.jpg']);
            if fnum==2
                basepic=imread([dirs.montage2ddir,filek(1).name]);
                basepicdims=size(basepic);
                basepic=uint8(basepic);
%                 
                imwrite(basepic,[dirs.montage3ddir,filek(1).name],'tif');
                imwrite(basepic,[dirs.montage3ddir,filek(1).name(1:end-4),'.jpg'],'jpg','Bitdepth',8,'Quality',jpgquality);
            else
%                 basepic=imread([dirs.montage3ddir,filek(fnum-1).name(1:end-4),'.jpg']);
                basepic=imread([dirs.montage3ddir,filek(fnum-1).name]);
                basepic=uint8(basepic);
            end
            basepic(basepic==0)=NaN;
            tempic=basepic;
            tempiclows=NaN([ceil(size(tempic)/downsize),downsize]);
            for it=1:downsize
                tempiclow=downsample(tempic,downsize,it-1);
                ttt=downsample(tempiclow',downsize,it-1)';
                tempiclows(1:size(ttt,1),1:size(ttt,2),it)=ttt;
            end
            tempiclow=nanmin(double(tempiclows),[],3);
            basepiclow=tempiclow;
            if isempty(a)
                newpic=imread([dirs.montage2ddir,fname]);
            else
                newpic=imread([dirs.montage3ddir,fname]);
%                 newpic=imread([dirs.montage3ddir,fname(1:end-4),'.jpg']);
            end
            newpicdims=size(newpic);
            newpic=uint8(newpic);
            newpic(newpic==0)=NaN;
            tempic=newpic;
            tempiclows=NaN([ceil(size(tempic)/downsize),downsize]);
            for it=1:downsize
                tempiclow=downsample(tempic,downsize,it-1);
                ttt=downsample(tempiclow',downsize,it-1)';
                tempiclows(1:size(ttt,1),1:size(ttt,2),it)=ttt;
            end
            tempiclow=nanmin(double(tempiclows),[],3);
            newpiclow=tempiclow;
            
            basepiclow(basepiclow==0)=NaN;
            newpiclow(newpiclow==0)=NaN;
            basepicloworig=basepiclow;
            newpicloworig=newpiclow;
            %     tempiclow=downsample(tempic,downsize);
            %     newpiclow=downsample(tempiclow',downsize)';
            %

            % equalizing the image
            nanIDX = find(isnan(newpiclow));
            iter=0;
            while iter<220
                
                iter=iter+1;
                if mod(iter,2)==0
                    nanIDX(end)=[];
                    newpiclow(nanIDX) = newpiclow(nanIDX+1);
                    nanIDX      = find(isnan(newpiclow));
                else
                    nanIDX(1)=[];
                    newpiclow(nanIDX) = newpiclow(nanIDX-1);
                    nanIDX      = find(isnan(newpiclow));
                end
                if mod(iter,10)==0
                    newpiclow=newpiclow';
                end
            end
            newpiclow=newpicloworig-imfilter(newpiclow,h,'replicate');
            
            nanIDX = find(isnan(basepiclow));
            iter=0;
            while  iter<220
                
                iter=iter+1;
                if mod(iter,2)==0
                    nanIDX(end)=[];
                    basepiclow(nanIDX) = basepiclow(nanIDX+1);
                    nanIDX      = find(isnan(basepiclow));
                else
                    nanIDX(1)=[];
                    basepiclow(nanIDX) = basepiclow(nanIDX-1);
                    nanIDX      = find(isnan(basepiclow));
                end
                if mod(iter,10)==0
                    basepiclow=basepiclow';
                end
            end
            basepiclow=basepicloworig-imfilter(basepiclow,h,'replicate');
            % equalizing the image
            temp1=double(basepiclow)-double((nanmedian(basepiclow(:))-sdvaldown*nanstd(basepiclow(:))));
            temp1=temp1/(nanmedian(temp1(:))+sdvalup*nanstd(temp1(:)));%max(temp1(:));
            temp1(temp1<0)=0;
            temp1(temp1>1)=1;
            temp2=double(newpiclow)-double((nanmedian(newpiclow(:))-sdvaldown*nanstd(newpiclow(:))));
            temp2=temp2/(nanmedian(temp2(:))+sdvalup*nanstd(temp2(:)));%max(temp2(:));
            temp2(temp2<0)=0;
            temp2(temp2>1)=1;
            
            
            ptsbase  = detectSURFFeatures(temp1,'MetricThreshold',mthresh,'NumOctaves',1,'NumScaleLevels',8);
            ptsnew = detectSURFFeatures(temp2,'MetricThreshold',mthresh,'NumOctaves',1,'NumScaleLevels',8);
            
            %A szélhez közeli talált cuccokat törli
            
            todel=zeros(size(ptsbase));  %
            for i=1:length(ptsbase)
                coord=ptsbase(i).Location;
                coordrund=round(coord);
                if any(coord<pixdif) | coord(2)>size(temp1,1)-pixdif | coord(1)>size(temp1,2)-pixdif
                    todel(i)=1;
                else
                    nearestxnan=min([find(isnan(temp1(coordrund(2):size(temp1,1),coordrund(1))),1,'first'),find(isnan(temp1(coordrund(2):-1:1,coordrund(1))),1,'first')]);
                    nearestynan=min([find(isnan(temp1(coordrund(2),coordrund(1):size(temp1,2))),1,'first'),find(isnan(temp1(coordrund(2),coordrund(1):-1:1)),1,'first')]);
                    addnum=(ptsbase(i).Scale*15+12)/3;
                    if nanmin([nearestxnan-addnum,nearestynan-addnum])<pixdif
                        todel(i)=1;
                    end
                end
            end
            ptsbase(find(todel))=[];
            deleted=length(find(todel));
            todel=zeros(size(ptsnew));
            for i=1:length(ptsnew)
                coord=ptsnew(i).Location;
                coordrund=round(coord);
                if any(coord<pixdif) | coord(2)>size(temp2,1)-pixdif | coord(1)>size(temp2,2)-pixdif% | ptsnew(i).Scale<15
                    todel(i)=1;
                else
                    nearestxnan=min([find(isnan(temp2(coordrund(2):size(temp2,1),coordrund(1))),1,'first'),find(isnan(temp2(coordrund(2):-1:1,coordrund(1))),1,'first')]);
                    nearestynan=min([find(isnan(temp2(coordrund(2),coordrund(1):size(temp2,2))),1,'first'),find(isnan(temp2(coordrund(2),coordrund(1):-1:1)),1,'first')]);
                    addnum=(ptsnew(i).Scale*15+12)/3;
                    if nanmin([nearestxnan-addnum,nearestynan-addnum])<pixdif
                        todel(i)=1;
                    end
                end
            end
            ptsnew(find(todel))=[];
            deleted=deleted+length(find(todel));
            disp(deleted)
            [featuresBase,  validPtsBase]  = extractFeatures(temp1,  ptsbase);
            [featuresNew, validPtsNew] = extractFeatures(temp2, ptsnew);
            inliernew=[];
            matchthresholdnow=matchthreshold;
            maxdistnow=maxdist;
            while length(inliernew)<9 & maxdistnow<100
                disp(['matchthreshold=',num2str(matchthresholdnow),'  maxdist=',num2str(maxdistnow)]);
                indexPairs = matchFeatures(featuresBase, featuresNew,'MatchThreshold',matchthresholdnow);
                if matchthresholdnow<=90
                    matchthresholdnow=matchthresholdnow+10;
                else
                    maxdistnow=maxdistnow+3;
                end
                matchedbase=validPtsBase(indexPairs(:,1));
                matchednew=validPtsNew(indexPairs(:,2));
                if length(matchedbase)>3
                    [tform, inliernew, inlierbase] = estimateGeometricTransform(matchednew, matchedbase, 'similarity','MaxDistance',maxdistnow);
                else
                    tform=[];
                    inliernew=[];
                    inlierbase=[];
                end
            end
            if ~exist('inlierbase','var') | isempty(inlierbase)
                return
            end
            tform.T(3,1)=downsize*tform.T(3,1);
            tform.T(3,2)=downsize*tform.T(3,2);
            rotationdegree=-1*double(atan2d( tform.T(1,2),tform.T(1,1)));
            
            if abs(rotationdegree)>1
                dotherotate=1;
            else
                dotherotate=0;
            end
                
            
%             if roundi==1
%                 dotherotate=1;
%             elseif roundi==2
%                 
%                 if abs(rotationdegree)>3
%                     if lastacceptedrotation~=0 & abs(rotationdegree-lastacceptedrotation)<3
%                         dotherotate=1;
%                     else
%                         figure(333)
%                         subplot(2,2,1)
%                         imagesc(temp1)
%                         title('base image image')
%                         subplot(2,2,3)
%                         imagesc(temp2)
%                         title('original image')
%                         subplot(2,2,4)
%                         imagesc(imrotate(temp2,rotationdegree))
%                         title('rotated image')
%                         button = questdlg(['Rotate the lower image with ',num2str(rotationdegree),' degrees?'],'ANSWER!','Yes','No','Try Again','Yes');
%                         if strcmp(button,'Yes')
%                             dotherotate=1;
%                             lastacceptedrotation=rotationdegree;
%                         elseif strcmp(button,'Try Again')
%                             fnum=fnum-1;
%                             dotherotate=0;
%                         else
%                             lastacceptedrotation=0;
%                             dotherotate=0;
%                         end
%                     end
%                 else
%                     lastacceptedrotation=0;
%                     dotherotate=0;
%                 end
%             elseif roundi==3
%                 tforms2(fnum).T=tform.T;
%                 tforms2(fnum).rotationdegree=rotationdegree;
%                 disp(['for ',filek(fnum).name, ' the residual rotation is ',num2str(rotationdegree)]);
%                 dotherotate=0;
%                  if abs(rotationdegree)>3
%                         fnum=fnum-1;
%                     end
%             end
            
            
            
            if dotherotate==1
                outpic=imrotate(newpic,rotationdegree);
                disp(['rotating: ',filek(fnum).name, '    ',num2str(rotationdegree),' degrees'])
            else
                outpic=newpic;
                tforms(fnum).T=tform.T;
                tforms(fnum).rotationdegree=rotationdegree;
                tforms(fnum).X=tforms(fnum).T(3,2);
                tforms(fnum).Y=tforms(fnum).T(3,1);
                disp(['for ',filek(fnum).name, ' the residual rotation is ',num2str(rotationdegree)]);
            end
            
             outbw=outpic>0;
                stats = regionprops(outbw,'BoundingBox','Area');
                [~,idx]=sort([stats.Area],'descend');
                stats=stats(idx(1));
                stats.BoundingBox=round(stats.BoundingBox);
                outpic=outpic(stats.BoundingBox(2)+1:stats.BoundingBox(2)+stats.BoundingBox(4)-1,stats.BoundingBox(1)+1:stats.BoundingBox(1)+stats.BoundingBox(3)-1);
                lastacceptedrotation=rotationdegree;
                %     outpic=imwarp(newpic,tform);
                
                if plotthestuff==1
                    figure(1)
                    clf
                    subplot(2,1,1)
                    imagesc(temp1)
                    hold on
                    plot(validPtsBase.selectStrongest(100),'showOrientation',true);
                    colormap gray
                    subplot(2,1,2)
                    imagesc(temp2)
                    hold on
                    plot(validPtsNew.selectStrongest(100),'showOrientation',true);
                    colormap gray
                    
                    figure(2)
                    clf
                    showMatchedFeatures(temp1,temp2,inlierbase,inliernew,'montage');
                    figure(3)
                    clf
                    subplot(2,1,1)
                    imagesc(basepic);
                    subplot(2,1,2)
                    imagesc(outpic)
                end
                
                outpic(isnan(outpic))=0;
                imwrite(outpic,[dirs.montage3ddir,filek(fnum).name],'tif');
                 imwrite(outpic,[dirs.montage3ddir,filek(fnum).name(1:end-4),'.jpg'],'jpg','Bitdepth',8,'Quality',jpgquality);
                if dotherotate==1
                    fnum=fnum-1;
                end
        %     return
        %         pause(3)
    end
end
return
%% generate and apply offset values then rotate and save image
tforms(1).Xoffset=0;
tforms(1).Yoffset=0;
for i=2:length(tforms)
    if abs(tforms(i).X)>10000
        tforms(i).Xoffset=round(tforms(i-1).Xoffset);
    else
        tforms(i).Xoffset=round(tforms(i-1).Xoffset-tforms(i).X);
    end
    if abs(tforms(i).Y)>10000
        tforms(i).Yoffset=round(tforms(i-1).Yoffset);
    else
        tforms(i).Yoffset=round(tforms(i-1).Yoffset-tforms(i).Y);
    end
end
for i=1:length(tforms)
    tforms(i).Xoffset=tforms(i).Xoffset-max([tforms.Xoffset]);
    tforms(i).Yoffset=tforms(i).Yoffset-max([tforms.Yoffset]);
end
figure(1234)
plot([tforms.Xoffset])
hold on
plot([tforms.Yoffset],'r-')
for fnum=1:length(filek)
     fname=filek(fnum).name;
      pic=imread([dirs.montage3ddir,filek(fnum).name]);
      pic=[uint8(zeros(abs(tforms(fnum).Xoffset),size(pic,2)));pic];
      pic=[uint8(zeros(size(pic,1),abs(tforms(fnum).Yoffset))),pic];
      pic=imrotate(pic,-90);
      imwrite(pic,[dirs.montage3ddir,filek(fnum).name],'tif')
      disp(['final file: ',filek(fnum).name])
end
return
%% üres maszkok legyártása
switches.pixelarea=.923^2;
switches.thresholdstart=0.5;
switches.minarea=1000;
switches.maxarea=3e+05;
switches.plotthestuff=1;
gauci0=[8,8];
gauci=[3,3];
h0=fspecial('Gaussian',gauci0*6,gauci(1));
h=fspecial('Gaussian',gauci*6,gauci(1));
filek=dir(dirs.montage2ddir);
filek([filek.isdir])=[];
for fnum=1:length(filek)
    fname=filek(fnum).name;
    %     filek(fnum).nyerskep=imread([tifkonyvtar,fname]);
    hyp=strfind(fname,'_');
    filek(fnum).ID=(fname(1:hyp(end)-1));%hyp(end-1)+
    filek(fnum).idx=str2num(fname(hyp(end)+1:end-4));
end
for fnum=1:length(filek)
    fname=filek(fnum).name;
    a=dir([dirs.maskdir,fname]);
    if isempty(a)
        basepic=imread([dirs.montage2ddir,fname]);
        basepic=double(basepic);
        basepic(basepic==0)=NaN;
        basepic=abs(basepic-nanmax(basepic(:)));
        tempic0=imfilter(basepic,h0,'replicate');
        tempic=imfilter(basepic,h,'replicate');
        dfperf0=(tempic-tempic0)./tempic0;
        tempic=dfperf0;%>nanmedian(dfperf0(:))+nanstd(dfperf0(:))*.5;
        threshold=switches.thresholdstart;
        SD=nanstd(tempic(:));
        MED=nanmedian(tempic(:));
        BW=tempic>MED+SD*threshold;
        BW=bwareaopen(BW,round(switches.minarea/switches.pixelarea));
        figure(1)
        clf
        subplot(2,1,1)
        imagesc(tempic)
        caxis([nanstd(dfperf0(:))*-2, nanstd(dfperf0(:))*2])
        subplot(2,1,2)
        imagesc(BW)
        return
        %     imwrite(BW,[dirs.maskdir,fname],'tif')
        
        
    end
end

%% teli maszkok legyártása
switches.pixelarea=.923^2;
switches.thresholdstart=0;
switches.minarea=100;
switches.maxarea=3e+05;
switches.plotthestuff=1;
gauci0=[15,15];
gauci=[3,3];
h0=fspecial('Gaussian',gauci0*6,gauci(1));
h=fspecial('Gaussian',gauci*6,gauci(1));
filek=dir(dirs.montage2ddir);
filek([filek.isdir])=[];
for fnum=1:length(filek)
    fname=filek(fnum).name;
    %     filek(fnum).nyerskep=imread([tifkonyvtar,fname]);
    hyp=strfind(fname,'_');
    filek(fnum).ID=(fname(1:hyp(end)-1));%hyp(end-1)+
    filek(fnum).idx=str2num(fname(hyp(end)+1:end-4));
end
for fnum=1:length(filek)
    fname=filek(fnum).name;
    a=dir([dirs.maskdir,fname]);
    if isempty(a)
        basepic=imread([dirs.montage2ddir,fname]);
        basepic=double(basepic);
        basepic(basepic==0)=NaN;
        basepic=abs(basepic-nanmax(basepic(:)));
        tempic0=imfilter(basepic,h0,'replicate');
        tempic=imfilter(basepic,h,'replicate');
        dfperf0=(tempic-tempic0)./tempic0;
        tempic=basepic;%>nanmedian(dfperf0(:))+nanstd(dfperf0(:))*.5;
        % %     figure(1)
        % %     clf
        % %     imagesc(tempic)
        % %     caxis([0, nanmax(dfperf0(:))])
        % %     return
        %         tempic=abs(tempic-nanmax(tempic(:)));
        threshold=switches.thresholdstart;
        SD=nanstd(tempic(:));
        MED=nanmedian(tempic(:));
        BW=tempic>MED+SD*threshold;
        BW=bwareaopen(BW,round(switches.minarea/switches.pixelarea));
        [mask,num]=bwlabel(BW,4);
        stats=regionprops(mask,'Area','BoundingBox','Centroid','ConvexArea','EquivDiameter','Extent','FilledArea','MajorAxisLength','Perimeter','PixelList','PixelIdxList','Solidity');
        areas=[stats.ConvexArea]*switches.pixelarea;
        roundnesses=[stats.EquivDiameter]./[stats.MajorAxisLength];
        soliditys=[stats.Solidity];
        extents=[stats.Extent];
        cellnum=length(stats);
        
        while max(areas)>switches.maxarea
            threshidxes=find(areas>switches.maxarea);
            for i=1:length(threshidxes)
                
                threshnow=threshold;
                
                newpic=(zeros(size(tempic)));
                newpic(stats(threshidxes(i)).PixelIdxList)=tempic(stats(threshidxes(i)).PixelIdxList);
                nBW=newpic>MED+SD*threshnow;
                
                newstats=regionprops(nBW,'ConvexArea','EquivDiameter','MajorAxisLength','PixelIdxList','Solidity','Extent');
                
                nnn=1;
                while sum([newstats.ConvexArea])*switches.pixelarea>switches.minarea & nnn==1
                    threshnow=threshnow+.25;
                    nBW=newpic>MED+SD*threshnow;
                    newstats=regionprops(nBW,'ConvexArea','EquivDiameter','MajorAxisLength','PixelIdxList','Solidity','Extent');
                    [L, nnn] = bwlabel(nBW, 4);
                    figure(2)
                    cla
                    imagesc(nBW)
                    title(num2str(threshnow))
                end
                allpixindlist=[];
                for iii=1:length(newstats)
                    allpixindlist=[allpixindlist;newstats(iii).PixelIdxList];
                end
                BW(stats(threshidxes(i)).PixelIdxList)=false;
                if sum([newstats.ConvexArea])*switches.pixelarea>switches.minarea
                    BW(allpixindlist)=true;
                end
                
                if switches.plotthestuff==1
                    figure(1)
                    clf
                    subplot(1,3,1)
                    imagesc(tempic);
                    colormap(linspecer)
                    title(['original'])
                    subplot(1,3,2)
                    imagesc(mask);
                    colormap(linspecer)
                    title(['mask - threshold: ',num2str(threshold)])
                    subplot(1,3,3)
                    IMAGEnew=tempic;
                    IMAGEnew(bwconvhull(BW,'objects',4)>0)=max(IMAGEnew(:));
                    imagesc(IMAGEnew);
                    colormap(linspecer)
                    title(['masked original'])
                    %             pause
                end
                
            end
            BW=bwareaopen(BW,round(switches.minarea/switches.pixelarea));
            [mask,num]=bwlabel(BW,4);
            stats=regionprops(mask,'Area','BoundingBox','Centroid','ConvexArea','EquivDiameter','Extent','FilledArea','MajorAxisLength','Perimeter','PixelList','PixelIdxList','Solidity');
            areas=[stats.ConvexArea]*switches.pixelarea;
            roundnesses=[stats.EquivDiameter]./[stats.MajorAxisLength];
            soliditys=[stats.Solidity];
            extents=[stats.Extent];
        end
        imwrite(BW,[dirs.maskdir,fname],'tif')
        
        %         temp1=bwconvhull(BW,'objects',4);
        %         ptsbase  = detectSURFFeatures(temp1,'MetricThreshold',mthresh,'NumOctaves',4,'NumScaleLevels',7);
        %         [featuresBase,  validPtsBase]  = extractFeatures(temp1,  ptsbase);
        %         figure(1231);
        %         clf
        %         imshow(temp1);
        %         hold on;
        %         plot(validPtsBase.selectStrongest(500),'showOrientation',true);
    end
end
%% z montázs
mthresh=1000;

filek=dir(dirs.maskdir);
filek([filek.isdir])=[];
for fnum=1:length(filek)
    fname=filek(fnum).name;
    %     filek(fnum).nyerskep=imread([tifkonyvtar,fname]);
    hyp=strfind(fname,'_');
    filek(fnum).ID=(fname(1:hyp(end)-1));%hyp(end-1)+
    filek(fnum).idx=str2num(fname(hyp(end)+1:end-4));
end

for fnum=2:length(filek)
    
    fname=filek(fnum).name;
    if fnum==2
        basepic=imread([dirs.maskdir,filek(1).name]);
    end
    newpic=imread([dirs.maskdir,fname]);
    %         temp1=bwconvhull(basepic,'objects',8);
    %         temp2=bwconvhull(newpic,'objects',8);
    temp1=abs(1-imfill(basepic,'holes'));
    temp2=abs(1-imfill(newpic,'holes'));
    ptsbase  = detectSURFFeatures(temp1,'MetricThreshold',mthresh,'NumOctaves',15,'NumScaleLevels',4);
    ptsnew = detectSURFFeatures(temp2,'MetricThreshold',mthresh,'NumOctaves',15,'NumScaleLevels',4);
    [featuresBase,  validPtsBase]  = extractFeatures(temp1,  ptsbase);
    [featuresNew, validPtsNew] = extractFeatures(temp2, ptsnew);
    figure(1231);
    clf
    subplot(2,1,1)
    imshow(temp1);
    hold on;
    plot(validPtsBase.selectStrongest(500),'showOrientation',true);
    subplot(2,1,2)
    imshow(temp2);
    hold on;
    plot(validPtsNew.selectStrongest(500),'showOrientation',true);
    %         pause
    
    indexPairs = matchFeatures(featuresBase, featuresNew);%
    matchedbase=validPtsBase(indexPairs(:,1));
    matchednew=validPtsNew(indexPairs(:,2));
    [tform, inliernew, inlierbase] = estimateGeometricTransform(matchednew, matchedbase, 'similarity');
    figure(1)
    clf
    showMatchedFeatures(temp1,temp2,inlierbase,inliernew,'montage');
    hold on
    
    pause
end


