function keposszerakosdi_2D_main(filek,dirs,plotthestuff)
%% xy montázs - új megközelítés
zvals=unique({filek.z});
zpix=struct;
gauci=[1 1];
gaucinagy=[5 5];
h=fspecial('gaussian',gauci*6,gauci(1));
hnagy=fspecial('gaussian',gaucinagy*6,gaucinagy(1));
stdszorzo=4;
% percentileboundaries=[.01 99.99];
for znum=1:length(zvals)
    idxes=find(strcmp({filek.z},zvals(znum)));
    [~,ix]=sort([filek(idxes).idx]);
    idxes=idxes(ix);
    fname=filek(idxes(1)).name;
    hyps=strfind(fname,'_');
    a=dir([dirs.montage2ddir,zvals{znum},'.tif']);
    tformsmatrix=cell(length(idxes));
    matchedpairs=zeros(length(idxes));
    if isempty(a) & length(idxes)>=2
        for picnum=1:length(idxes)-1
            for picnum2=picnum+1:length(idxes)
                disp([fname,'   -  picture ',num2str(picnum),'/',num2str(length(idxes)), 'VS ',  'picture ',num2str(picnum2),'/',num2str(length(idxes))]);
                fname=filek(idxes(picnum)).name;
                basepic=imread([dirs.tifkonyvtar,fname]);
                basepic=double(basepic);
                fname=filek(idxes(picnum2)).name;
                newpic=imread([dirs.tifkonyvtar,fname]);
                newpic=double(newpic);
                dosecondround=true;
                [tformsmatrix{picnum,picnum2},matchedpairs(picnum,picnum2)]=keposszerakosdi_givetransformationforpicpair(basepic,newpic,stdszorzo,h,dosecondround,hnagy,plotthestuff);
            end
        end
        %% relativ helyek kiszamitasa
        for picnum1=1:length(idxes)
            for picnum2=1:length(idxes)
                    if isstruct(tformsmatrix{picnum2,picnum1})
                        tformsmatrix{picnum1,picnum2}.tform=tformsmatrix{picnum2,picnum1}.tform.invert;
                        tformsmatrix{picnum1,picnum2}.matchedpairs=matchedpairs(picnum2,picnum1);
                        tformsmatrix{picnum2,picnum1}.matchedpairs=matchedpairs(picnum2,picnum1);
                        if isobject( tformsmatrix{picnum2,picnum1}.tform2)
                            tformsmatrix{picnum1,picnum2}.tform2=tformsmatrix{picnum2,picnum1}.tform2.invert;
                        else
                            tformsmatrix{picnum1,picnum2}.tform2=tformsmatrix{picnum2,picnum1}.tform2;
                        end
                    end
            end
        end
        %% új megközelítés, először a leginkább összepasszolókat teszem össze
        matchedpairss=matchedpairs;
        doneidxes=[];
        
        symstext='syms';
        for i=1:size(matchedpairs,1)
%             symstext=[symstext,' x',num2str(i)];
            symstext=[symstext,' x',num2str(i),' y',num2str(i)];
        end
        eval(symstext);
        eqnx1 = x1==0;
        eqny1 = y1==0;
        NEXT=1;
        while length(doneidxes)<size(matchedpairs,1)
            prevNEXT=NEXT;
            matchedpairsss=matchedpairss;
            while prevNEXT==NEXT
                [aa,idx]=nanmax(matchedpairsss(:));
                
                [a,b]=ind2sub(size(matchedpairsss),idx);
                
                if isempty(doneidxes) | (~any(a==doneidxes) && any(b==doneidxes)) | (any(a==doneidxes) && ~any(b==doneidxes))
                    xdifi=-(tformsmatrix{a,b}.tform.T(3,1)-1*tformsmatrix{a,b}.tform2.T(3,1));
                    ydifi=-(tformsmatrix{a,b}.tform.T(3,2)-1*tformsmatrix{a,b}.tform2.T(3,2));
                    NEXT=NEXT+1;
                    eval(['eqnx',num2str(NEXT),'= x',num2str(a),'-x',num2str(b),'==',num2str(xdifi),';']);
                    eval(['eqny',num2str(NEXT),'= y',num2str(a),'==y',num2str(b),'+',num2str(ydifi),';']);
                    doneidxes=unique([doneidxes,a,b]);
                end
                matchedpairsss(idx)=NaN;
            end
            matchedpairss(idx)=NaN;
        end
        
        solstrx='equationsToMatrix([';
        solstry='equationsToMatrix([';
        for ii=1:NEXT;
            solstrx=[solstrx,'eqnx',num2str(ii),', '];
            solstry=[solstry,'eqny',num2str(ii),', '];
        end
        solstrx=[solstrx(1:end-2),'], ['];
        solstry=[solstry(1:end-2),'], ['];
        for i=1:size(matchedpairs,1)
%             solstr=[solstr,'x',num2str(i),',  '];
            solstrx=[solstrx,'x',num2str(i),', '];
            solstry=[solstry,'y',num2str(i),', '];
        end
        solstrx=[solstrx(1:end-2),']);'];
        solstry=[solstry(1:end-2),']);'];
        [A,B]=eval(solstrx);
        X = linsolve(A,B);
        [A,B]=eval(solstry);
        Y = linsolve(A,B);
        loc=(double([X,Y]));
         loc=(bsxfun(@(x,y) x-y,loc,min(loc)))+1;
        loc=round(loc);
        montage=keposszerakosdi_mergepicturesbyrelativelocations(dirs,filek,idxes,loc,hnagy,plotthestuff);
% % % % % % % % % % %         %%
% % % % % % % % % % %          locations=struct;
% % % % % % % % % % %          for i=1:length(idxes)
% % % % % % % % % % %             locations(i).ncoords=0;
% % % % % % % % % % %             locations(i).sourcepic=[];
% % % % % % % % % % %             locations(i).hasacoord=false;
% % % % % % % % % % %             locations(i).order=NaN;
% % % % % % % % % % %          end
% % % % % % % % % % %           locations(1).coord(1,:)=[0,0];
% % % % % % % % % % %           locations(1).ncoords=1;
% % % % % % % % % % %            locations(1).hasacoord=true;
% % % % % % % % % % %             locations(1).sourcepic=NaN;
% % % % % % % % % % %             locations(1).order=1;
% % % % % % % % % % %             proba=0;
% % % % % % % % % % %             while proba<5 | (sum([locations.hasacoord])<length(locations) & proba<1000)
% % % % % % % % % % %                 proba=proba+1;
% % % % % % % % % % %                 idxestocheck=find([locations.hasacoord]);
% % % % % % % % % % %                 for picnum2i=1:length(idxestocheck)
% % % % % % % % % % %                     picnum2=idxestocheck(picnum2i);
% % % % % % % % % % %                     for picnum1=1:length(idxes)
% % % % % % % % % % %                          if isstruct(tformsmatrix{picnum2,picnum1})% & ~locations(picnum1).hasacoord
% % % % % % % % % % % %                              disp('yaaay')
% % % % % % % % % % %                               if isobject(tformsmatrix{picnum2,picnum1}.tform2)
% % % % % % % % % % %                                 locnow=tformsmatrix{picnum2,picnum1}.tform.T(3,1:2)+tformsmatrix{picnum2,picnum1}.tform2.T(3,1:2)+ mean(locations(picnum2).coord,1);
% % % % % % % % % % %                                 rellocnow=tformsmatrix{picnum2,picnum1}.tform.T(3,1:2)+tformsmatrix{picnum2,picnum1}.tform2.T(3,1:2);
% % % % % % % % % % %                               else
% % % % % % % % % % %                                   locnow=tformsmatrix{picnum2,picnum1}.tform.T(3,1:2) + mean(locations(picnum2).coord,1);
% % % % % % % % % % %                                   rellocnow=tformsmatrix{picnum2,picnum1}.tform.T(3,1:2);
% % % % % % % % % % %                               end
% % % % % % % % % % %                               if ~any( picnum2==[locations(picnum1).sourcepic])
% % % % % % % % % % %                                    locations(picnum1).ncoords=locations(picnum1).ncoords+1;
% % % % % % % % % % %                                    NEXT=locations(picnum1).ncoords;
% % % % % % % % % % %                                     locations(picnum1).hasacoord=true;
% % % % % % % % % % %                                     if isnan(locations(picnum1).order)
% % % % % % % % % % %                                         locations(picnum1).order=nanmax([locations.order])+1;
% % % % % % % % % % %                                     end
% % % % % % % % % % %                                     locations(picnum1).sourcepic(NEXT)=picnum2;
% % % % % % % % % % %                                     locations(picnum1).coord(NEXT,:)=locnow;
% % % % % % % % % % %                                     locations(picnum1).relcoord(NEXT,:)=rellocnow;
% % % % % % % % % % %                                     locations(picnum1).matchedpairnum(NEXT,:)=tformsmatrix{picnum2,picnum1}.matchedpairs';
% % % % % % % % % % %                                     
% % % % % % % % % % %                               end
% % % % % % % % % % %                          end
% % % % % % % % % % %                          progressbar(sum([locations.hasacoord])/length([locations.hasacoord]))
% % % % % % % % % % %                     end
% % % % % % % % % % %                 end
% % % % % % % % % % %             end
% % % % % % % % % % %             for i=1:length(locations)
% % % % % % % % % % %                 locations(i).numberofneighbours=sum(~isnan(locations(i).sourcepic));
% % % % % % % % % % %             end
% % % % % % % % % % %             if proba<1000
% % % % % % % % % % %             order=[locations.order];
% % % % % % % % % % %             %%
% % % % % % % % % % %             loc=[];
% % % % % % % % % % %             idxesnew=[];
% % % % % % % % % % %             clear locationsnew
% % % % % % % % % % % % % % % %             for i=1:length(order)
% % % % % % % % % % % % % % % %                 idxesnew(i)=idxes(find(order==i));
% % % % % % % % % % % % % % % %                 locationsnew(i)=locations(find(order==i));
% % % % % % % % % % % % % % % %             end
% % % % % % % % % % % % % % % %             idxes=idxesnew;
% % % % % % % % % % % % % % % %             locations=locationsnew;
% % % % % % % % % % % %             loc=zeros(length(locations),2);
% % % % % % % % % % %             for ii=1:length(locations)
% % % % % % % % % % %                 i=find(order==ii);
% % % % % % % % % % % %                 loc(i,:)=median(locations(i).coord,1); %régi cucc
% % % % % % % % % % %                 if ii==1
% % % % % % % % % % %                     loc(i,:)=[1,1];
% % % % % % % % % % %                 else
% % % % % % % % % % %                     kellhetapic=(locations(i).sourcepic);
% % % % % % % % % % %                     matchedpairok=(locations(i).matchedpairnum);
% % % % % % % % % % %                     picrellocs=locations(i).relcoord;
% % % % % % % % % % %                     neededpic=find(order(kellhetapic)<ii);
% % % % % % % % % % %                     picrellocs=picrellocs(neededpic,:);
% % % % % % % % % % %                     matchedpairok=matchedpairok(neededpic,:);
% % % % % % % % % % %                     kellhetapic=kellhetapic(neededpic);
% % % % % % % % % % %                     locationok=[];
% % % % % % % % % % %                     matchedpairss=[];
% % % % % % % % % % %                     for prevpici=1:length(kellhetapic)
% % % % % % % % % % %                         locationok(prevpici,:)=loc(prevpici,:)+picrellocs(prevpici,:);
% % % % % % % % % % %                         matchedpairss(prevpici)=matchedpairok(prevpici);
% % % % % % % % % % %                     end
% % % % % % % % % % %                     [~,neededpair]=max(matchedpairok);
% % % % % % % % % % %                     loc(i,:)=locationok(neededpair,:);
% % % % % % % % % % %                 end
% % % % % % % % % % %             end
% % % % % % % % % % % %             locations=loc;
% % % % % % % % % % %             fname=filek(idxes(1)).name;
% % % % % % % % % % %         image=imread([dirs.tifkonyvtar,fname]);
% % % % % % % % % % %         loc=(bsxfun(@(x,y) x-y,loc,min(loc)))+1;
% % % % % % % % % % %         %          locations=(bsxfun(@(x,y) x+y,locations,size(image)));
% % % % % % % % % % %         loc=round(loc);
% % % % % % % % % % %        
% % % % % % % % % % %         montage=keposszerakosdi_mergepicturesbyrelativelocations(dirs,filek,idxes,loc,hnagy,plotthestuff);
  
        outpic=montage-min(montage(:));
        outline=outpic;
        outline(find(isnan(outline)))=[];
%             outline=deleteoutliers(outline,.05);
                outline=sort(outline);
figure(122)
        clf
        subplot(2,1,1)
        hist(outline,1000)
%             outline=outline(1:round(length(outline)*.99));
        outline=outline(round(length(outline)*.001):round(length(outline)*.999));
        subplot(2,1,2)
        hist(outline,1000)
        outpic=(outpic-min(outline(:)))/(max(outline(:))-min(outline));

        imwrite(uint8(floor(outpic*2^8)),[dirs.montage2ddir,zvals{znum},'.tif'],'tif')
        %     zpix(znum).image=basepic;
        elseif isempty(a) & length(idxes)==1
         fname=filek(idxes(1)).name;
         basepic=imread([dirs.tifkonyvtar,fname]);
            basepic=double(basepic);
            outpic=basepic-min(basepic(:));
        outline=outpic;
        
        outline(find(isnan(outline)))=[];
        outline=sort(outline);
        outline=outline(round(length(outline)*.001):round(length(outline)*.999));
        
        outpic=(outpic-min(outline(:)))/(max(outline(:))-min(outline));
        
        imwrite(uint8(floor(outpic*2^8)),[dirs.montage2ddir,zvals{znum},'.tif'],'tif')
            else
                disp('Couldn''t merge pictures... or already done. :)')
                
            end
    end

end
% %%
% NEXT=0;
% for i=1:length(locations)
%     for j=1:length(locations(i).sourcepic)
%         if ~isnan(locations(i).sourcepic(j))
%             NEXT=NEXT+1;
%             eval(['eqn',num2str(NEXT),'= x',num2str(i),'==x',num2str(locations(i).sourcepic(j)),'+',num2str(locations(i).relcoord(j,1))]);
%         end
%     end
% end
% solstr='solve([';
% for ii=1:NEXT;
%     solstr=[solstr,'eqn',num2str(ii),', '];
% end
% solstr=[solstr(1:end-2),'], ['];
% for i=1:length(locations)
%     solstr=[solstr,'x',num2str(i),', '];
% end
% solstr=[solstr(1:end-2),']);'];
% sol=eval(solstr);
%%
% %% old version
% %% xy montázs - összes kép egymás után sorban
% zvals=unique({filek.z});
% zpix=struct;
% gauci=[1 1];
% gaucinagy=[5 5];
% h=fspecial('gaussian',gauci*6,gauci(1));
% hnagy=fspecial('gaussian',gaucinagy*6,gaucinagy(1));
% stdszorzo=4;
% for znum=1:length(zvals)
%     idxes=find(strcmp({filek.z},zvals(znum)));
%     [~,ix]=sort([filek(idxes).idx]);
%     idxes=idxes(ix);
%     fname=filek(idxes(1)).name;
%     hyps=strfind(fname,'_');
%     a=dir([dirs.montage2ddir,zvals{znum},'.tif']);
%     if isempty(a) & length(idxes)>=2
%         for picnum=2:length(idxes)
%             disp([fname,'  picture ',num2str(picnum),'/',num2str(length(idxes))]);
%             fname=filek(idxes(picnum-1)).name;
%             %         basepic=filek(idxes(picnum-1)).nyerskep;
%             %         newpic=filek(idxes(picnum)).nyerskep;
%             if picnum==2
%                 basepic=imread([dirs.tifkonyvtar,fname]);
%                 basepic=double(basepic);
%             end
%             fname=filek(idxes(picnum)).name;
%             newpic=imread([dirs.tifkonyvtar,fname]);
%             newpic=double(newpic);
%             
%          basepic=keposszerakosdi_mergeoverlappingpix(basepic,newpic,stdszorzo,h,hnagy,plotthestuff);
%         end
%          
%         outpic=basepic-min(basepic(:));
%         outline=outpic;
%         outline(find(isnan(outline)))=[];
%         %     outline=deleteoutliers(outline,.05);
%         outline=sort(outline);
%         %     outline=outline(1:round(length(outline)*.95));
%         outpic=outpic/max(outline(:));
%         
%         imwrite(uint8(floor(outpic*2^8)),[dirs.montage2ddir,zvals{znum},'.tif'],'tif')
%         %     zpix(znum).image=basepic;
%     elseif isempty(a) & length(idxes)==1
%          fname=filek(idxes(1)).name;
%          basepic=imread([dirs.tifkonyvtar,fname]);
%             basepic=double(basepic);
%             outpic=basepic-min(basepic(:));
%         outline=outpic;
%         outline(find(isnan(outline)))=[];
%         %     outline=deleteoutliers(outline,.05);
%         outline=sort(outline);
%         %     outline=outline(1:round(length(outline)*.95));
%         outpic=outpic/max(outline(:));
%         
%         imwrite(uint8(floor(outpic*2^8)),[dirs.montage2ddir,zvals{znum},'.tif'],'tif')
%     end
% end