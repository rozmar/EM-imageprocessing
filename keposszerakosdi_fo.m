%% kepek beolvasasa
clear all
close all
plotthestuff=0;
locations=marcicucca_locations;
dirs.alapkonyvtar=[locations.EMdir,'EMdata/'];
dirs.tifkonyvtar=[dirs.alapkonyvtar,'rawtif/Elmi_patosz/juci/0410181mg/a dendrit/'];
dirs.montage2ddir=[dirs.alapkonyvtar,'2dmontage/Elmi_patosz/juci/0410181mg/a dendrit/'];
dirs.maskdir=[dirs.alapkonyvtar,'2dmasks/Elmi_patosz/juci/0410181mg/a dendrit/'];
% alapkonyvtar='/home/rozmar/Mount/TGTAR_1/ANALYSISdata/marci/_elmi/';
filek=dir(dirs.tifkonyvtar);
filek([filek.isdir])=[];
todel=zeros(size(filek));
for fnum=1:length(filek)
    if ~(strfind(filek(fnum).name,'.tif')+3==length(filek(fnum).name)) | any(strfind(filek(fnum).name,'hely'))
        todel(fnum)=1;
    end
end
filek(find(todel))=[];
for fnum=1:length(filek)
    fname=filek(fnum).name;
    %     filek(fnum).nyerskep=imread([tifkonyvtar,fname]);
    hyp=strfind(fname,'_');
    filek(fnum).z=(fname(1:hyp(end)-1));%hyp(end-1)+
    filek(fnum).idx=str2num(fname(hyp(end)+1:end-4));
    
end
%% xy montázs
zvals=unique({filek.z});
zpix=struct;
gauci=[1 1];
gaucinagy=[5 5];
        h=fspecial('gaussian',gauci*6,gauci(1));
        hnagy=fspecial('gaussian',gaucinagy*6,gaucinagy(1));
 stdszorzo=4;    
for znum=1:length(zvals)
    idxes=find(strcmp({filek.z},zvals(znum)));
    fname=filek(idxes(1)).name;
    hyps=strfind(fname,'_');
    a=dir([dirs.montage2ddir,fname(1:hyps(end)-1),'.tif']);
    if isempty(a)
    for picnum=2:length(idxes)
        disp([fname,'  picture ',num2str(picnum),'/',num2str(length(idxes))]);
        fname=filek(idxes(picnum-1)).name;
        %         basepic=filek(idxes(picnum-1)).nyerskep;
        %         newpic=filek(idxes(picnum)).nyerskep;
        if picnum==2
            basepic=imread([dirs.tifkonyvtar,fname]);
            basepic=double(basepic);
        end
        fname=filek(idxes(picnum)).name;
        newpic=imread([dirs.tifkonyvtar,fname]);
        newpic=double(newpic);
        % %         basepic0=zeros(size(basepic)+size(newpic)*2);
        % %         basepic0(size(newpic,1)+1:size(newpic,1)+size(basepic,1),size(newpic,2)+1:size(newpic,2)+size(basepic,2))=basepic;
        % %         basepic=basepic0;
        % %         c = normxcorr2(newpmaxdist=1;ic,basepic);
        
        
        STDval=nanstd(basepic(:));
        temp1=double(basepic)-double((nanmedian(basepic(:))-stdszorzo*STDval));%double(nanmin(basepic(:)));%;
        STDval=nanstd(temp1(:));
        temp1=temp1/(nanmedian(temp1(:))+stdszorzo*STDval);%max(temp2(:));nanmax(temp1(:));%
        temp1=imfilter(temp1,h,'replicate');
        STDval=nanstd(newpic(:));
        temp2=double(newpic)-double((nanmedian(newpic(:))-stdszorzo*STDval));%double(nanmin(newpic(:)));%
        STDval=nanstd(temp2(:));
        temp2=temp2/(nanmedian(temp2(:))+stdszorzo*STDval);%nanmax(temp2(:));%%max(temp2(:));
        temp2=imfilter(temp2,h,'replicate');
%         
%     temp1=double(basepic)-double(nanmin(basepic(:)));
%         temp1=temp1/(nanmedian(temp1(:))+2*nanstd(temp1(:)));%max(temp1(:));
%         temp2=double(newpic)-double(nanmin(newpic(:)));
%         temp2=temp2/(nanmedian(temp2(:))+2*nanstd(temp2(:)));%max(temp2(:));
    
        temp1(temp1<0)=0;
    temp1(temp1>1)=1;
    
    temp2(temp2<0)=0;
    temp2(temp2>1)=1;
    
        
        mthresh=1000;
        inlierbase=[];
        matchedbase=1;
        inlinerbases=1;
        matchthreshold=5;
         maxdist=1;
        while ((length(inlierbase)<13 | length(inlierbase)>100) & (400<=mthresh)) & length(inlierbase)<1000 & length(inlierbase)/length(matchedbase)<.9
            if length(inlierbase)<13
                if maxdist<=5
                    maxdist=maxdist+1;
                elseif matchthreshold+10>100
                    mthresh=mthresh-100;
                else  
                    matchthreshold=matchthreshold+10;
                end
            else
                matchthreshold=matchthreshold/1.5;
%                 mthresh=mthresh+50;
            end
            ptsbase  = detectSURFFeatures(temp1,'MetricThreshold',mthresh,'NumOctaves',3,'NumScaleLevels',9);
            ptsnew = detectSURFFeatures(temp2,'MetricThreshold',mthresh,'NumOctaves',3,'NumScaleLevels',9);
            [featuresBase,  validPtsBase]  = extractFeatures(temp1,  ptsbase);
            [featuresNew, validPtsNew] = extractFeatures(temp2, ptsnew);
            indexPairs = matchFeatures(featuresBase, featuresNew,'MatchThreshold',matchthreshold);
            matchedbase=validPtsBase(indexPairs(:,1));
            matchednew=validPtsNew(indexPairs(:,2));
            if length(matchedbase)>1
                [tform, inliernew, inlierbase] = estimateGeometricTransform(matchednew, matchedbase, 'similarity','MaxDistance',maxdist);
            else
                inlinerbase=[];
            end
            disp([mthresh,matchthreshold,maxdist,length(matchedbase),length(inlierbase)])
            
        end
        Tinv  = tform.invert.T;
        xchange=floor(Tinv(3,2));
        ychange=floor(Tinv(3,1));
            
        base0=NaN(max(size(basepic,1),size(newpic,1))+abs(xchange),max(size(basepic,2),size(newpic,2))+abs(ychange),2);
        if ychange<0
            basey=1;
            newy=abs(ychange)+1;
        else
            basey=abs(ychange)+1;
            newy=1;
        end
        if xchange<0
            basex=1;
            newx=abs(xchange)+1;
        else
            basex=abs(xchange)+1;
            newx=1;
        end
        base0(basex:basex+size(basepic,1)-1,basey:basey+size(basepic,2)-1,1)=basepic;
        base0(newx:newx+size(newpic,1)-1,newy:newy+size(newpic,2)-1,2)=newpic;
        xchangenow=5000;
        ychangenow=5000;
        while abs(xchangenow)+abs(ychangenow)>2
            % második körben pontosítás
            overlapbaseixes=~isnan(base0(:,:,1))&~isnan(base0(:,:,2));
            [row, col] = ind2sub(size(base0(:,:,1)), find(overlapbaseixes));
            overlappingbase=(base0(min(row):max(row),min(col):max(col),1));
            overlappingnew=(base0(min(row):max(row),min(col):max(col),2));
            
            overlappingbase=double(overlappingbase)-double(nanmin(overlappingbase(:)));
            overlappingbase=overlappingbase/(nanmedian(overlappingbase(:))+2*nanstd(overlappingbase(:)));%max(overlappingbase(:));
            
            overlappingnew=double(overlappingnew)-double(nanmin(overlappingnew(:)));
            overlappingnew=overlappingnew/(nanmedian(overlappingnew(:))+2*nanstd(overlappingnew(:)));%max(overlappingnew(:));
            
            ptsbase2  = detectSURFFeatures(overlappingbase,'MetricThreshold',mthresh,'NumOctaves',3,'NumScaleLevels',12);
            ptsnew2 = detectSURFFeatures(overlappingnew,'MetricThreshold',mthresh,'NumOctaves',3,'NumScaleLevels',12);
            [featuresBase2,  validPtsBase2]  = extractFeatures(overlappingbase,  ptsbase2);
            [featuresNew2, validPtsNew2] = extractFeatures(overlappingnew, ptsnew2);
            indexPairs2 = matchFeatures(featuresBase2, featuresNew2,'MatchThreshold',matchthreshold);
            matchedbase2=validPtsBase2(indexPairs2(:,1));
            matchednew2=validPtsNew2(indexPairs2(:,2));
            [tform2, inliernew2, inlierbase2] = estimateGeometricTransform(matchednew2, matchedbase2, 'similarity');
            Tinv2  = tform2.invert.T;
            xchangenow=floor(Tinv2(3,2));
            ychangenow=floor(Tinv2(3,1));
            xchange=xchange+xchangenow;
            ychange=ychange+ychangenow;
            disp([floor(Tinv2(3,2)), floor(Tinv2(3,1))])
%             base0=NaN(size(basepic,1)+abs(xchange),size(basepic,2)+abs(ychange),2);
            base0=NaN(max(size(basepic,1),size(newpic,1))+abs(xchange),max(size(basepic,2),size(newpic,2))+abs(ychange),2);
            if ychange<0
                basey=1;
                newy=abs(ychange)+1;
            else
                basey=abs(ychange)+1;
                newy=1;
            end
            if xchange<0
                basex=1;
                newx=abs(xchange)+1;
            else
                basex=abs(xchange)+1;
                newx=1;
            end
            base0(basex:basex+size(basepic,1)-1,basey:basey+size(basepic,2)-1,1)=basepic;
            base0(newx:newx+size(newpic,1)-1,newy:newy+size(newpic,2)-1,2)=newpic;
        end
        
        
        % két kép intenzitásának összefésülése
        overlapbaseixes=~isnan(base0(:,:,1))&~isnan(base0(:,:,2));
        [row, col] = ind2sub(size(base0(:,:,1)), find(overlapbaseixes));
        overlappingbase=(base0(min(row):max(row),min(col):max(col),1));
        overlappingnew=(base0(min(row):max(row),min(col):max(col),2));
        
        if plotthestuff==1
            figure(11)
            clf
            plot(overlappingnew,overlappingbase,'kx')
            hold on
        end
        
        overlapdims=size(overlappingnew);
%         eddigx=round(overlapdims(1)/1);
%         eddigy=round(overlapdims(2)/1);
        overlappingnew=imfilter(overlappingnew,hnagy);
        overlappingbase=imfilter(overlappingbase,hnagy);
%         overlappingnew=overlappingnew(1:eddigx,1:eddigy);
%         overlappingbase=overlappingbase(1:eddigx,1:eddigy);
        overlappingnew(1)=[];
        overlappingbase(1)=[];
        idxtodel=find(isnan(overlappingnew)|isnan(overlappingbase));
        overlappingnew(idxtodel)=[];
        overlappingbase(idxtodel)=[];
%         overlappingnew= deleteoutliers(overlappingnew, .05, 1);
%         overlappingbase= deleteoutliers(overlappingbase, .05, 1);
%         idxtodel=find(isnan(overlappingnew)|isnan(overlappingbase));
%         overlappingnew(idxtodel)=[];
%         overlappingbase(idxtodel)=[];
%         sdc=2;
%         idxtodel=find((overlappingnew)>median(overlappingnew)+std(overlappingnew)*sdc|(overlappingnew)<median(overlappingnew)-std(overlappingnew)*sdc|(overlappingbase)>median(overlappingbase)+std(overlappingbase)*sdc|(overlappingbase)<median(overlappingbase)-std(overlappingbase)*sdc);
%         overlappingnew(idxtodel)=[];
%         overlappingbase(idxtodel)=[];
        if plotthestuff==1
            plot(overlappingnew,overlappingbase,'rx')
        end
        p=polyfit(overlappingnew,overlappingbase,1);
        preverz=polyfit(overlappingbase,overlappingnew,1);
        if plotthestuff==1
            xek=[round(min(overlappingnew)):round(max(overlappingnew))];
            plot(xek,polyval(p,xek),'k-','LineWidth',3)
            ylabel('intensity of old file')
            xlabel('intensity of new file')
            figure(12312322)
            clf
            subplot(2,1,1)
            hist(overlappingbase(:),1000)
            title('base histogram')
            subplot(2,1,2)
            hist(overlappingnew(:),1000)
            title('new histogram')
        end
        difi=overlappingbase-overlappingnew;
        diff=mean(overlappingbase(:))-mean(overlappingnew(:));
        if plotthestuff==1
            axx=[min(overlappingbase(:)):(max(overlappingbase(:))-min(overlappingbase(:)))/1000:max(overlappingbase(:))];
            figure(22)
            subplot(2,2,1)
            hist(overlappingbase(:),axx)
%             caxis(cmap)
            subplot(2,2,2)
            hist(overlappingnew(:),axx)
%             caxis(cmap)
            subplot(2,2,3)
            hist(difi(:))
            %                 caxis(cmap)
            subplot(2,2,4)
            hist(polyval(p,overlappingnew),axx)
%             caxis(cmap)
        end               %
         base0(:,:,2)=polyval(p,base0(:,:,2));
%          base0(:,:,1)=polyval(preverz,base0(:,:,1));
        %képek szép összefésülése
%         kozeprow=round(mean([min(col),max(col)]));
%         if ychange>0
%             base0(:,kozeprow:max(col),2)=NaN;
%             base0(:,min(col):kozeprow-1,1)=NaN;
%         else
%             base0(:,kozeprow:max(col),1)=NaN;
%             base0(:,min(col):kozeprow-1,2)=NaN;
%         end
        % új képet rátesszük a régi tetejére
        idxtonan=find(~isnan(base0(:,:,2)));
        tt=base0(:,:,1);
        tt(idxtonan)=NaN;
        base0(:,:,1)=tt;
        base=nanmean(base0,3);
%         base2=base0(:,:,2);
%         base(~isnan(base0(:,:,2)))=base2(~isnan(base0(:,:,2)));
if plotthestuff==1
    figure(1)
    clf
    subplot(2,2,1)
    imagesc(basepic)
    title('basepic')
    subplot(2,2,2)
    imagesc(newpic)
    title('newpic')
    subplot(2,2,3)
    showMatchedFeatures(temp1,temp2,inlierbase,inliernew,'montage');
    
    
    %         showMatchedFeatures(temp1,temp2,matchedbase,matchednew);
    %         imagesc(c)
    subplot(2,2,4)
    imagesc(base)
    title('merge')
    %                 pause
end
        basepic=base;
        [x,y]=ind2sub(size(basepic),find(~isnan(basepic)));
        basepic=basepic(min(x):max(x),min(y):max(y));
    end
    
    outpic=basepic-min(basepic(:));
    outline=outpic;
    outline(find(isnan(outline)))=[];
    outline=deleteoutliers(outline,.05);
    outpic=outpic/max(outline(:));
    hyps=strfind(fname,'_');
    imwrite(uint16(floor(outpic*2^16)),[dirs.montage2ddir,fname(1:hyps(end)-1),'.tif'],'tif')
%     zpix(znum).image=basepic;
end
end
return
%% mintavételezés változtatása
sdvalup=2;
sdvaldown=2;
gauci=[3,3];
pixdif=100;
downsize=30;
pixdif=ceil(pixdif/downsize);
mthresh=200;
matchthreshold=.1;
switches.pixelarea=.923^2;
switches.thresholdstart=0;
switches.plotthestuff=1;
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
for fnum=2:length(filek)
    fname=filek(fnum).name;
    a=dir([dirs.maskdir,fname]);
    if fnum==2
        clear tempiclows
        basepic=imread([dirs.montage2ddir,filek(2).name]);
        basepic(basepic==0)=NaN;
        tempic=imfilter(basepic,h,'replicate');
        tempiclows=NaN([ceil(size(tempic)/downsize),downsize]);
        for it=1:downsize
            tempiclow=downsample(tempic,downsize,it-1);
            ttt=downsample(tempiclow',downsize,it-1)';
            tempiclows(1:size(ttt,1),1:size(ttt,2),it)=ttt;
        end
        tempiclow=nanmin(double(tempiclows),[],3);
        basepiclow=tempiclow;
    end
    newpic=imread([dirs.montage2ddir,fname]);
    newpic(newpic==0)=NaN;
    tempic=imfilter(newpic,h,'replicate');
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
    
%     tempiclow=downsample(tempic,downsize);
%     newpiclow=downsample(tempiclow',downsize)';
%      
    temp1=double(basepiclow)-double((nanmedian(basepiclow(:))-sdvaldown*nanstd(basepiclow(:))));
    temp1=temp1/(nanmedian(temp1(:))+sdvalup*nanstd(temp1(:)));%max(temp1(:));
    temp1(temp1<0)=0;
    temp1(temp1>1)=1;
    temp2=double(newpiclow)-double((nanmedian(newpiclow(:))-sdvaldown*nanstd(newpiclow(:))));
    temp2=temp2/(nanmedian(temp2(:))+sdvalup*nanstd(temp2(:)));%max(temp2(:));
    temp2(temp2<0)=0;
    temp2(temp2>1)=1;
    
    
    ptsbase  = detectSURFFeatures(temp1,'MetricThreshold',mthresh,'NumOctaves',2,'NumScaleLevels',30);
    ptsnew = detectSURFFeatures(temp2,'MetricThreshold',mthresh,'NumOctaves',2,'NumScaleLevels',30);

    %A szélhez közeli talált cuccokat törli

        todel=zeros(size(ptsbase));  %
        for i=1:length(ptsbase)
            coord=ptsbase(i).Location;
            coordrund=round(coord);
            if any(coord<pixdif) | coord(2)>size(temp1,1)-pixdif | coord(1)>size(temp1,2)-pixdif
                todel(i)=1;
            else
                nearestxnan=min([find(temp1(coordrund(2):size(temp1,1),coordrund(1))==0,1,'first'),find(temp1(coordrund(2):-1:1,coordrund(1))==0,1,'first')]);
                nearestynan=min([find(temp1(coordrund(2),coordrund(1):size(temp1,2))==0,1,'first'),find(temp1(coordrund(2),coordrund(1):-1:1)==0,1,'first')]);
                addnum=ptsbase(i).Scale*6+9;
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
                nearestxnan=min([find(temp2(coordrund(2):size(temp2,1),coordrund(1))==0,1,'first'),find(temp2(coordrund(2):-1:1,coordrund(1))==0,1,'first')]);
                nearestynan=min([find(temp2(coordrund(2),coordrund(1):size(temp2,2))==0,1,'first'),find(temp2(coordrund(2),coordrund(1):-1:1)==0,1,'first')]);
                addnum=ptsnew(i).Scale*6+9;
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
    indexPairs = matchFeatures(featuresBase, featuresNew,'MatchThreshold',matchthreshold);
    matchedbase=validPtsBase(indexPairs(:,1));
    matchednew=validPtsNew(indexPairs(:,2));
    if length(matchedbase)>1
        [tform, inliernew, inlierbase] = estimateGeometricTransform(matchednew, matchedbase, 'similarity');
    else
        inlinerbase=[];
    end
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

return
end
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


