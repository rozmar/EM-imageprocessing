function basepic=keposszerakosdi_mergeoverlappingpix(basepic,newpic,stdszorzo,h,hnagy,plotthestuff)

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
                
                ptsbase2  = detectSURFFeatures(overlappingbase,'MetricThreshold',mthresh,'NumOctaves',1,'NumScaleLevels',12);
                ptsnew2 = detectSURFFeatures(overlappingnew,'MetricThreshold',mthresh,'NumOctaves',1,'NumScaleLevels',12);
                [featuresBase2,  validPtsBase2]  = extractFeatures(overlappingbase,  ptsbase2);
                [featuresNew2, validPtsNew2] = extractFeatures(overlappingnew, ptsnew2);
                indexPairs2 = matchFeatures(featuresBase2, featuresNew2,'MatchThreshold',matchthreshold);
                matchedbase2=validPtsBase2(indexPairs2(:,1));
                matchednew2=validPtsNew2(indexPairs2(:,2));
                %%
                baselocations=[matchedbase2.Location];
                xdistfromedge=min([abs(baselocations(:,1)-size(overlappingbase,2)),baselocations(:,1)],[],2);
                ydistfromedge=min([abs(baselocations(:,2)-size(overlappingbase,1)),baselocations(:,2)],[],2);
                basedistfomedge=min(xdistfromedge,ydistfromedge);
                newlocations=[matchednew2.Location];
                xdistfromedge=min([abs(newlocations(:,1)-size(overlappingbase,2)),newlocations(:,1)],[],2);
                ydistfromedge=min([abs(newlocations(:,2)-size(overlappingbase,1)),newlocations(:,2)],[],2);
                newdistfomedge=min(xdistfromedge,ydistfromedge);
                distfomedge=min(basedistfomedge,newdistfomedge);
                matchednew2=matchednew2(distfomedge>median(distfomedge));
                matchedbase2=matchedbase2(distfomedge>median(distfomedge));
                %%
                %             pause
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
        