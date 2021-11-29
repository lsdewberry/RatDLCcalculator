%sciread group averages - xyf
% only used for outputs formatted like Hdutyfactor (has a scalar and a std
% for each trial)
%% load in/initialize
close all;
%load('sciread_output.mat')
clearvars -except output
%select variable
var2examine = 'avgvelocity';
spotcheck=1;
fnames = fieldnames(output);
if ~isfield(output,var2examine);
    fprintf('field not found. Field names are:\n')
    disp(fnames)
    error('var2examine field not found')
end
    
%% filter all data to only needed data
%fnames = fieldnames(output);
Allfields = struct2cell(output);  %Creates a numfield x M x N cell array
% filter by err field 
filtindex =strcmp(Allfields(strcmp('Err', fieldnames(output)),:,:),'gud');
Allfields = Allfields(:,:,filtindex);
% filter by direction
filtindex =strcmp(Allfields(strcmp('direction', fieldnames(output)),:,:),'right');
Allfields = Allfields(:,:,filtindex);
% filter by number of cycles
filtindex =cell2mat(Allfields(strcmp('numcycles', fieldnames(output)),:,:))>=3;
Allfields = Allfields(:,:,filtindex);
% filter by velocity std field 
filtindex =cell2mat(Allfields(strcmp('velocitystd', fieldnames(output)),:,:))<=10;
Allfields = Allfields(:,:,filtindex);
% filter by variable selected
clear filtindex
fnames = {var2examine, 'rat','group','stim','daysPO','pixpercm'};
for f=1:length(fnames)
filtindex(f,:) =strcmp(fnames{f}, fieldnames(output));
end
filtindex=find(sum(filtindex));
Allfields = Allfields(filtindex,:,:);
% Allfields(feilds,1,trials)includes every rat,tp, etc
% rename rat name to ratname_s if it's a stimulated trial
Sindex = find(cell2mat(squeeze(Allfields(strcmp('stim',fnames),:,:))));
Allfields(strcmp('rat',fnames),:,Sindex)= strcat(Allfields(strcmp('rat',fnames),:,Sindex),'s');
Allfields(strcmp('group',fnames),:,Sindex)= strcat(Allfields(strcmp('group',fnames),:,Sindex),'s');
%clear output %commented out for easy testing

%% group trials for each rat/timepoint
colors=cbrewer('qual','Paired',12);
rats   = unique(Allfields(strcmp(   'rat', fnames),:,:));% find rat names
groups = unique(Allfields(strcmp( 'group', fnames),:,:));% find group names
tps    = unique(Allfields(strcmp('daysPO', fnames),:,:));% find daysPO timepoints
tps=tps([1,4,2,3]);
vars = {'x','y','f'};
for r=1:length(rats)
    Rindex  =strcmp(Allfields(strcmp(   'rat', fnames),:,:),rats{r});%find trials for selected rat_s
        for tp=1:length(tps)
            TPindex =strcmp(Allfields(strcmp('daysPO', fnames),:,:),tps{tp});%find trials where daysPO is tp
            filtindex = Rindex & TPindex;
            if sum(filtindex)==0
                continue
            end
            mrtp      = Allfields(strcmp(var2examine, fnames),:,filtindex); %trials of selected rat at selected timepoint and stim condition
            for trial = 1:size(mrtp,3)
                %mpixpercm = find(filtindex,trial);
                %mpixpercm = Allfields{strcmp( 'pixpercm', fnames),:,mpixpercm(trial)};
                mrtp1=mrtp{trial};
                rtp{r,tp,trial,1}=abs(mrtp1.avg);%/mpixpercm;
            end
            clear mrtp 
            for v=1
                 
                 [nrtp{r,tp,v}]=mean(cell2mat(rtp   (r,tp,:,v)));
                 %[nrtpstd{r,tp,v},~]        =avgtrial(rtp(r,tp,:,3),rtpstd(r,tp,:,v),0,colors(tp*2+v-2,:),[rats{r},' ',tps{tp},' ',vars{v}]);%get std
                 %[Y, F]=avgtrial(F,Y,graph?,color,name)
                 %nrtp(r,tp,percentofwhole,v)
                 %need to add std so that nrtp(r,tp,v,avg/std)
                 mystd=std(cell2mat(rtp   (r,tp,:,v)));
                 nrtpG{r}=Allfields{strcmp('group', fnames),:,find(filtindex,1,'first')};
                 if 0 %spotcheck
                     figure(v); hold on;
                     errorbar(tp,nrtp{r,tp,v},mystd,'o','DisplayName',[rats{r},' ',tps{tp}])
                     fig1=figure(v); hold on; xlabel('timepoint');legend('Location','eastoutside');
                     fig1.Position=[1925+(v-1)*850,80,850,800-40];set(gca,'LooseInset',get(gca,'TightInset'));
                     %errorbar(nrtp{r,tp,3},nrtp{r,tp,v},nrtpstd{r,tp,v},'Color',colors(tp*2+v-2,:),[rats{r},' ',tps{tp},' ',vars{v}]);
                 end
            end
        
        end
        if 0 %spotcheck
            xlim([0 5]);
            y1=ylim;  y1(1)=0; ylim(y1);
            %saveas(fig1,[rats{r},num2str(v),'velocity_avg.png']);
            pause(3);
            %spotcheck
        end
        close all
end

%% group into treatments and plot xvf,yvf,xvy
colors2(:,:,1)=cbrewer('seq','BuGn',6);
colors2(:,:,2)=cbrewer('seq','OrRd',6);
colors2(:,:,3)=cbrewer('seq','Purples',6);
colors2(:,:,4)=cbrewer('seq','Greys',6);
colors2=colors2(2:end-1,:,:);

for g=1:length(groups)
    Gindex = strcmp(nrtpG,groups{g});%find trials where group is g
    for tp=1:length(tps)
%         if sum(filtindex)==0
%             continue
%         end

        for v=1
            mgrptp{v} = cell2mat(nrtp(Gindex,tp,v));
        end
        
        for v=1
            grptp{g,tp, v} = nanmean(mgrptp{v});
            grptpstd{g,tp,v} = std(mgrptp{v});
        end
        if 1 %g==1||g==2 %0 %spotcheck
            figure(2);hold on;
            %errorbar(tp+g/8-.5,grptp{g,tp,1},grptpstd{g,tp,1},'o','Color',colors2(tp,:,g),'DisplayName',[groups{g},' ',tps{tp}],'LineWidth',1); hold on; legend('Location','eastoutside');
            errorbar(tp+g/8-.5,grptp{g,tp,1},grptpstd{g,tp,1},'o','Color',colors2(4,:,g),'DisplayName',[groups{g}],'LineWidth',1); hold on; legend('Location','eastoutside');
            xlim([0 5]);
            y1=ylim;  y1(1)=0;            ylim(y1);
            %figure(3);errorbar(tp,grptp{g,tp,2},grptpstd{g,tp,2},'o','Color',colors2(tp,:,g),'DisplayName',[groups{g},' ',tps{tp}],'LineWidth',1); hold on; legend('Location','eastoutside');
            %figure(4);errorbar(tp,grptp{g,tp,2},grptpstd{g,tp,1},grptpstd{g,tp,1},grptpstd{g,tp,2},grptpstd{g,tp,2},'Color',colors2(tp,:,g),'DisplayName',[groups{g},' ',tps{tp},' ',vars{2}],'LineWidth',1); hold on; legend('Location','eastoutside');
        end
    end
end
if 1 %spotcheck
figure(2);hold on; legend;xlabel('timepoint');ylabel(var2examine);
end

%% cool graphs
colors3(:,:,1)=cbrewer('seq','BuGn',20);
colors3(:,:,2)=cbrewer('seq','OrRd',20);
colors3(:,:,3)=cbrewer('seq','Purples',20);
colors3(:,:,4)=cbrewer('seq','Greys',20);
colors3(:,:,5)=cbrewer('seq','Blues',20);
colors3=colors3(2:end-1,:,:);
 %grptp{g,tp,v}
 %pre
gs=1;
coolgraphs(cell2mat(grptp(gs,2,1)),cell2mat(grptp(gs,2,2)),cell2mat(grptp(gs,2,3)),colors3(:,:,1))
title('pre');
%on
gs=2;
coolgraphs(cell2mat(grptp(gs,2,1)),cell2mat(grptp(gs,2,2)),cell2mat(grptp(gs,2,3)),colors3(:,:,2))
title('on');
%np sham
gs=3;
coolgraphs(cell2mat(grptp(gs,2,1)),cell2mat(grptp(gs,2,2)),cell2mat(grptp(gs,2,3)),colors3(:,:,3))
title('np sham');
%sham sham
gs=4;
coolgraphs(cell2mat(grptp(gs,2,1)),cell2mat(grptp(gs,2,2)),cell2mat(grptp(gs,2,3)),colors3(:,:,4))
title('sham sham');
%bln
gs=[1,3,4];
coolgraphs(nanmean(cell2mat(grptp(gs,1,1))),nanmean(cell2mat(grptp(gs,1,2))),nanmean(cell2mat(grptp(gs,1,3))),colors3(:,:,5))
title('baseline');

%%
function [] = coolgraphs(x,y,col,mycolormap)
%initialize
z = zeros(size(x));
figure
colormap(mycolormap)
%graph
surface([x;x], [y;y], [z;z], [col;col],...
    'facecol', 'no',...
    'edgecol', 'interp',...
    'linew',3);
ylim([0,.6]);xlim([-1,9])
%
end

function [datasnew, percentxsnew] = avgtrial(dataf,data, spotcheck,color,name)
    
    winsize = 5;%size in percentage of stride
    %% plot as percent of stride
    for j=1:size(dataf,3)%iterate over each trial=j
        percentf = dataf{j};
        mytrial  = data {j};
        if spotcheck %plot each trial
            hold on;
            plot(percentf, mytrial,':','Color',color/j,'DisplayName',[name,' trial ',num2str(j)])
            %plot(percentf, mytrial,':','Color',[color,1],'HandleVisibility','off') %plot individual strides
        end
        winnum=1;%index for which window bin we are on
        for k = 0:winsize:100-winsize %window moving over data in trial to get only points in winsize % of stride
%             winstridey = mytrialy(percentf>=k & percentf<k+winsize);
%             winstridex = mytrialx(percentf>=k & percentf<k+winsize);
%             winf       = percentf(percentf>=k & percentf<k+winsize);
            
            wintrials{j,winnum} =  mytrial(percentf>=k & percentf<k+winsize);%save strides
            winfs     {j,winnum}= percentf(percentf>=k & percentf<k+winsize);%save f values for strides
            winnum=winnum+1;
        end
    end
    %j
    if exist('wintrials', 'var')==0
        pause
    end
    msg=['testing sizes line 420, should be true: ', num2str(size(wintrials,2)==100/winsize)];
    %initialize some things
    nwintrials=cell(size(wintrials,1),1);
    nwinfs=cell(size(wintrials,1),1);
    weights=cell(size(wintrials,1),1);
    for m=1:size(wintrials,2) %loop thru the number of bins that exist. should = 100/winsize
        wins=wintrials(:,m);%select only bin m for all strides
        minps = min(cellfun('size',wins,1));
        if minps<=0
            %msg=['need to increase winsize? winsize=',num2str(winsize)]
            %or maybe exclude the stride missing the data?
            %could decide based on stride length outliers?
            minps=1;
        end
        for n=1:size(wintrials,1) %loop thru the strides (n) for each bin (m).
            winnm=wintrials{n,m}; %select the bin in question
            winfnm=winfs{n,m}; %select fs for the bin in question
            if length(winnm)>minps
                newis = round(linspace(1,length(winnm),minps));
                nwinnm = winnm (newis);
                nwinfnm = winfnm(newis);
                nwintrials{n} = [nwintrials{n}, nwinnm'];%
                nwinfs{n}      = [     nwinfs{n},nwinfnm'];%
                weights{n}     = [    weights{n},ones(size(nwinnm'))];
            elseif length(winnm)==minps
                nwintrials{n} = [nwintrials{n}, winnm'];%
                nwinfs{n}      = [     nwinfs{n},winfnm'];%
                weights{n}     = [    weights{n},ones(size(winnm'))];
            else %if i don't have enough data i'm interpolating.
                msg2{n,m} = ['no data for ',name,' trial ',num2str(n),', percentf ',num2str(m)];
                winmminus1 = nan; winmplus1  = nan; winnminus1=nan;winfallns=nan;
                try    winmminus1 = wintrials{n,m-1};  end %same stride bin previous
                try    winmplus1  = wintrials{n,m+1};  end %same stride bin next
                try    winnminus1 = nwintrials{n}(end);  end %the avg bin we just appended (this works, i just think it's better to nan these)
                for smalln=1:size(winfs,1) %all stride f for bin
                    winfallns = [winfallns,(winfs{smalln,m})];      
                end 
                nwinnm = nanmean([ winmminus1, winmplus1,winnminus1]);%what if mminus one and plus 1 are both empty??
                nwinfnm= nanmean([winfallns]);
                
                if isnan(mean([nwinnm,nwinfnm])) %catch error
                    msg=(['no data for ',name,' trial ',num2str(n),', percentf ',num2str(m)]);
                    weights{n}     = [    weights{n},nan];%don't care abt these points
                else
                    weights{n}     = [    weights{n},0.5*ones(size(nwinnm))];%don't care abt these points
                end
                nwintrials{n} = [nwintrials{n}, nwinnm];%
                nwinfs{n}      = [     nwinfs{n},nwinfnm];%
                
                %alternative - make winsize higher? or include 2 bins for this point?
            end
            if sum([size(nwintrials{n})==size(nwinfs{n}),size(nwintrials{n})==size(weights{n}),-4])
            error('yo our stuff is not lining up right')
            end
        end
    end
    nwintrials = cell2mat(nwintrials); nwinfs = cell2mat(nwinfs); weights = cell2mat(weights); 
    % nanmean with weights for data
    medianstrides = median(nwintrials,1);    
    test = abs(nwintrials-medianstrides);  
    test = sum(test,2);
    [test,TFtest] = rmoutliers(test,'median');
    weights = weights(~TFtest,:);
    avnwinstrides = nansum(nwintrials(~TFtest,:).*weights,1)./nansum(weights,1);
    avwinstridesSTDEV = sqrt(nansum(weights.*((nwintrials(~TFtest,:)-avnwinstrides).^2))./(nansum(weights,1)-(nansum(weights.^2,1)/nansum(weights))));

    % nanmean with weights for percentf
    avnwinfs = nansum(nwinfs(~TFtest,:).*weights,1)./nansum(weights,1);
    avwinfsSTDEV =sqrt(nansum(weights.*((nwinfs(~TFtest,:)-avnwinfs).^2))./(nansum(weights,1)-(nansum(weights.^2,1)/nansum(weights))));
    
    if spotcheck
        hold on;
        errorbar(avnwinfs', avnwinstrides',avwinstridesSTDEV','Color',color,'LineWidth',2,'DisplayName',[name,' avg'])
    end
    
    datasnew=avnwinstrides;
    percentxsnew=avnwinfs;
end