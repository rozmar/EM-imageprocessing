function montage=keposszerakosdi_mergepicturesbyrelativelocations(dirs,filek,idxes,locations,hnagy,plotthestuff)
fname=filek(idxes(1)).name;
    image=imread([dirs.tifkonyvtar,fname]);     
montage=NaN([ceil(max(locations(:,2))-min(locations(:,2))+size(image,1)),ceil(max(locations(:,1))-min(locations(:,1))+size(image,2))]+1);
        for picnum=1:length(idxes)
            fname=filek(idxes(picnum)).name;
            image=imread([dirs.tifkonyvtar,fname]);
            if picnum==1
                montage([1:size(image,1)]+locations(picnum,2),[1:size(image,2)]+locations(picnum,1))=image;
            else
                montagenew=NaN(size(montage));
                montagenew([1:size(image,1)]+locations(picnum,2),[1:size(image,2)]+locations(picnum,1))=image;
                overlapbaseixes=~isnan(montage)&~isnan(montagenew);
                [row, col] = ind2sub(size(montage), find(overlapbaseixes));
                overlappingbase=(montage(min(row):max(row),min(col):max(col)));
                overlappingnew=(montagenew(min(row):max(row),min(col):max(col)));
                
                
                % két kép intenzitásának összefésülése
                
                if plotthestuff==1
                    figure(11)
                    clf
                    plot(overlappingnew,overlappingbase,'kx')
                    hold on
                end
                
                overlapdims=size(overlappingnew);
                overlappingnew=imfilter(overlappingnew,hnagy);
                overlappingbase=imfilter(overlappingbase,hnagy);
                overlappingnew(1)=[];
                overlappingbase(1)=[];
                idxtodel=find(isnan(overlappingnew)|isnan(overlappingbase));
                overlappingnew(idxtodel)=[];
                overlappingbase(idxtodel)=[];
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
                end
                
                
                
                
                %
                image=polyval(p,double(image));
                montage([1:size(image,1)]+locations(picnum,2),[1:size(image,2)]+locations(picnum,1))=image;
                progressbar(picnum/length(idxes))
            end
        end