%% choose path

clear variables;
%clearvars -except output <- if you want to keep output and only edit a few
%trials in output you can use this.
%scipath = 'D:\Baseline_vagotomy_carlos\all DLC excel files';
DLCOutputpath = 'F:\Baseline_vagotomy_carlos\all DLC excel files';
%% 
% The output path is your Current Folder in matlab

 
%% variable setup

files = dir(DLCOutputpath);
files = files(3:end);%should be 3 - because the first two entries from 'dir' are '...' invalid 
%% 
% I ended up just getting pixels/frame in this code and converting it later, 
% except for the rough velocity calculation at the end. (this is not the velocity 
% I actually analyze for final graphs. It's just there to give me an idea which 
% trials have issues).

pixeltocm=15; %this is a rough value used to calculate a rough velocity. The real pixtocm varied for each day.
framestosecs = 500/1; %frames/sec
%% 
% What likelyhood value (output from DeepLabCut) you will accept. A higher value 
% is stricter. Between 0 and 1.

L=0.9; %tolerance cutoff for DLC's likelyhood metric. 0.9 works well, but it can be changed if needed.
%% 
% % define colors for graphing

colorsP = cbrewer('qual', 'Paired', 8+4);
colorsS = cbrewer('qual','Set2',4+1);
colors = [colorsP(1:8,:);colorsS(1:4,:);colorsP(9:12,:);colorsS(5,:)];
%% 
% % filters to smooth jitter from the line for finding FSTO.

d1 = designfilt('lowpassiir','FilterOrder',3,'HalfPowerFrequency',0.05,'DesignMethod','butter');
d2 = designfilt('lowpassiir','FilterOrder',3,'HalfPowerFrequency',0.04,'DesignMethod','butter'); %this is the stricter one lol
%% 
% % initialize output structure

output(size(files, 1)) = struct();%if you want to edit only some of the output trials you should comment out this line
%% The loop over each trial begins
% edit this if you want to only run some of the files. default is 1:size(files,1);

for filenum = 2:5 %edit this if you want to only run some of the files. default is 1:size(files,1);
    i=filenum; %for convienience in typing coding only
    close all %close windows that were opened on the previous loop
    spotcheck = true; %Loud or quiet - do you want it to graph stuff. can also change to %round(rand(1)*rand(1)); to only graph some trials
    pausevalue = false; %if I want to pause after graphing each few lines. Otherwise it graphs all quick.
% import excel files

    myfile = files(filenum).name;
    disp(['filenum ',num2str(filenum), ': ',myfile]);%prints the file we are on in the command window
    output(i).filname = myfile; %I will have a massive 'output' structure at the end of this. I really need to initialize output tbh.
    scifile = [DLCOutputpath, '\', myfile];
    [frame, nosex,nosey,nosel,RightForex,RightForey,RightForel,LeftForex,LeftForey,LeftForel,RightHindToex,RightHindToey,RightHindToel,RightHindMidx, RightHindMidy, RightHindMidl, RightHindHeelx,RightHindHeely,RightHindHeell,LeftHindToex,LeftHindToey,LeftHindToel,LeftHindMidx, LeftHindMidy, LeftHindMidl, LeftHindHeelx,LeftHindHeely,LeftHindHeell,MirrorRightHindx,MirrorRightHindy,MirrorRightHindl,MirrorLeftHindx,MirrorLeftHindy,MirrorLeftHindl,Backx,Backy,Backl,TailBasex,TailBasey,TailBasel,ABforcePlateCenterx,ABforcePlateCentery,ABforcePlateCenterl,CDforcePlateCenterx,CDforcePlateCentery,CDforcePlateCenterl] = sciimport(scifile);
    
%% 
% this is where we begin analysis - i put a 'try' statement here so that the 
% code will continue to the next trial if it runs into an error with a trial. 
% This is useful so that you can run it overnight and then troubleshoot/get rid 
% of trials with errors in the morning.

try
% find direction of travel using nose

    Nx= nosex(nosel>L); %filter by likelyhood
    if length(Nx)<1
        output(i).byEye='length(Nx)<1';
        continue %stop doing this trial if there are no nose points
    end
    if Nx(end)>Nx(1)
        direction = 'right'; %if the animal is going to the right, I trust its right paw
        directionfactor=-1; %used so I don't have a bimodal distribution because of direction of travel
    elseif Nx(end)<Nx(1)
        direction = 'left '; %if the animal is going to the left, I trust its left paw
        directionfactor=1; %used so I don't have a bimodal distribution because of direction of travel
    else
        error('direction?')
    end
    output(i).direction = direction;
% setup figure

    if spotcheck == 1
        fig1=figure(1); hold on; %figure for entire trial
        xlabel('frames in video');
        legend('Location','eastoutside');
        fig1.Position=[0,40,1000,1000-40];
        set(gca,'LooseInset',get(gca,'TightInset'));
        % the force plate points 
          plot(median(frame(ABforcePlateCenterl>(L))),median(ABforcePlateCenterx(ABforcePlateCenterl>(L)))/10,'o','DisplayName','ABforcePlateCenter/10');
          plot(median(frame(CDforcePlateCenterl>(L))),median(CDforcePlateCenterx(CDforcePlateCenterl>(L)))/10,'o','DisplayName','CDforcePlateCenter/10');
        fig2 = figure(2); hold on; % figure for average stride
        xlabel('% stride from frame with toe off to next frame with toe off');
        legend('Location','eastoutside');
        fig2.Position=[1000,40,1000-40,1000-40];
        set(gca,'LooseInset',get(gca,'TightInset'));
    end
% Rename variables 
% I make it so it only counts the one in front. (one of my main outcomes is 
% foot location/angle in swing, and I don't trust the toe height of the back foot 
% because the body is in the way.)

    if strcmp(direction,'left ')
        %im just going to rename variables
        ForeX = LeftForex;        ForeY = LeftForey;        ForeL = LeftForel;
        HindToeX = LeftHindToex;  HindToeY = LeftHindToey;  HindToeL = LeftHindToel;
        HindMidX = LeftHindMidx;  HindMidY = LeftHindMidy;  HindMidL = LeftHindMidl;
        HindHeelX = LeftHindHeelx;HindHeelY = LeftHindHeely;HindHeelL = LeftHindHeell;
    elseif strcmp(direction,'right')
        ForeX = RightForex;        ForeY = RightForey;        ForeL = RightForel;
        HindToeX = RightHindToex;  HindToeY = RightHindToey;  HindToeL = RightHindToel;
        HindMidX = RightHindMidx;  HindMidY = RightHindMidy;  HindMidL = RightHindMidl;
        HindHeelX = RightHindHeelx;HindHeelY = RightHindHeely;HindHeelL = RightHindHeell;
    else
        error('direction?')
    end
% smoothing the lowest value Y things 
% to try to find the floor from them

    HTx = HindToeX(HindToeL>L); HTf=frame(HindToeL>L); color=colors(1,:); %hind toe x
    HTy = (-HindToeY(HindToeL>L)); %y is inverted because origin is top left
    [sHTx,sHTy,sHTf] = basiccmooth(HTx, HTy, HTf, d2);
    %Hind midfoot smoothing
    HMx= HindMidX(HindMidL>L); HMf=frame(HindMidL>L);
    HMy= -HindMidY(HindMidL>L); 
    [sHMx,sHMy,sHMf] = basiccmooth(HMx, HMy, HMf, d1);
    %Hind Heel smoothing
    HHx= HindHeelX(HindHeelL>L); HHf=frame(HindHeelL>L);
    HHy = -HindHeelY(HindHeelL>L);
    [sHHx,sHHy,sHHf] = basiccmooth(HHx, HHy, HHf, d1);%using og filter
%% find floor

    if false%check here if your floor does not slope
            floorfit = median( [min(sHTy),min(sHMy),min(sHHy)]);
            floorfit = polyfit([1,2],[floorfit,floorfit],1); %just make a straight line at y=floor
            evals(i) = "floor";
    else
            [floorfit] = floorfind(sHTy,sHMy,sHHy);
    end
       sHTy=sHTy- (floorfit(1)*(1:length(sHTy))+floorfit(2))'; 
%% FSTO
% this is where we segment the trial
% 
% most tweaks will happen in the fsto  function
% 
% If you are using AGATHA FSTO segmenting, make sure it is in the right format 
% (cell array with each row as a trial, in each cell is a matrix with footstrikes 
% in column 1 and toeoffs in column 2. The last entry may often be 0 or NaN.

    segmentingMethod = 2
    if segmentingMethod == 1 %if you want to import FSTO values (in frames) from other sources (e.g. AGATHA)
        %define strides based off of FSTOs from AGATHA format here?
    elseif segmentingMethod == 2
        [output(i).Hdutyfactor.avg, frameStrideLength, output(i).numcycles, strides, output(i).Hdutyfactor.std,output(i).Hstridelength.std, output(i).byEye,FSTOs]=fstoY(sHTy,sHTf, directionfactor,spotcheck, 1,[0 0 0], 'HT FSTO'); 
        output(i).FSTOs=FSTOs; %last one may be 0 if you don't see a toeoff after the last footstrike
        FSTOsFullStrides=FSTOs(FSTOs(:,2)~=0,:);
        output(i).stanceDuration.avg = mean(FSTOsFullStrides(:,2)-FSTOsFullStrides(:,1));
        output(i).stanceDuration.std =  std(FSTOsFullStrides(:,2)-FSTOsFullStrides(:,1));
    elseif segmentingMethod == 3
        [output(i).Hdutyfactor.avg, frameStrideLength, output(i).numcycles, strides, output(i).Hdutyfactor.std,output(i).Hstridelength.std, output(i).byEye,FSTOs]=fsto(sHTx,sHTf, directionfactor,spotcheck, 1,[0 0 0], 'HT FSTO'); 
        output(i).FSTOs=FSTOs; %last one may be 0 if you don't see a toeoff after the last footstrike
        FSTOsFullStrides=FSTOs(FSTOs(:,2)~=0,:);
        output(i).stanceDuration.avg = mean(FSTOsFullStrides(:,2)-FSTOsFullStrides(:,1));
        output(i).stanceDuration.std =  std(FSTOsFullStrides(:,2)-FSTOsFullStrides(:,1));
    else 
        error('choose a valid segmenting method in the dropdown off the FSTO section');
    end
% fsto plotting and stridelength calculations

    if size(strides,1)<1 || sum(sum((isnan(strides))))
          output(i).numcycles     = output(i).numcycles;
          if spotcheck
            saveas(fig1,[myfile,'_fulltrial.png']);
          end
          i=i+1;
          continue %go to next iteration of for loop because there aren't enough strides in this trial
    end
%% 
% %we need to find stridelength in the x data based on strides which was found 
% in f data

    strideIndex=NaN(size(strides,1),2);%initializing array
    for j=1:size(strides,1)
    [~,strideIndex(j,1)] = min(abs(sHTf-strides(j,1)));
    [~,strideIndex(j,2)] = min(abs(sHTf-strides(j,2)));
    end
    output(i).Hstridelength.avg = mean(abs(sHTx(strideIndex(:,1))-sHTx(strideIndex(:,2))));
    output(i).Hstridelength.std =  std(abs(sHTx(strideIndex(:,1))-sHTx(strideIndex(:,2))));
    stridesXstart=sHTx(unique(strideIndex));
%% analysis of sided data
% hind toe

      [output(i).HT.x.avg, output(i).HT.f.avg, output(i).HT.x.std,output(i).HT.f.std] = avgforstride (sHTx,sHTf, strides, stridesXstart,spotcheck,color,'HindToeX'); 
        if spotcheck
            figure(1)
            plot(sHTf, sHTx/10,'Color',color,'DisplayName','HindToeX/10')
            if pausevalue
            pause
            end
        end
       color=colors(2,:);%rand(1,3);
      [~,output(i).ToeProminance.avg,~,output(i).ToeProminance.std] = pulloutYstuff(sHTy,sHTf,spotcheck,color,'HindToeY');%mean 1/2 peak width and mean prominance
      [output(i).HT.y.avg, ~,output(i).HT.y.std,~] =   avgforstride(sHTy,sHTf, strides, [],spotcheck,color,'HindToeY'); 
        if spotcheck && pausevalue
            pause
        end
% midfoot processing

    color=colors(3,:);%rand(1,3);
      [output(i).HM.x.avg, output(i).HM.f.avg,output(i).HM.x.std,output(i).HM.f.std] = avgforstride (sHMx,sHMf, strides, stridesXstart,spotcheck,color,'HindMidX');
        if spotcheck
            figure(1);hold on;
            plot(sHMf,sHMx/10,'Color',color,'DisplayName','HindMidX/10')
            if pausevalue
            pause
            end
        end
     color=colors(4,:);
     sHMy=sHMy-((floorfit(1)*(1:length(sHMy))+floorfit(2)))';
      [output(i).HM.y.avg, ~,output(i).HM.y.std] = avgforstride (sHMy,sHMf, strides, [],spotcheck,color,'HindMidY'); 
        if spotcheck
            figure(1);hold on;
            plot(sHMf,sHMy,'Color',color,'DisplayName','HindMidY')
            if pausevalue
            pause
            end
        end
% Hind Heel

    color=colors(5,:);%rand(1,3);
      [output(i).HH.x.avg, output(i).HH.f.avg,output(i).HH.x.std,output(i).HH.f.std] = avgforstride (sHHx,sHHf, strides, stridesXstart,spotcheck,color,'HindHeelX');
        if spotcheck
            figure(1);
            %plot(frame(HindHeelL>L),HHx,'DisplayName','HindHeelX')
            plot(sHHf,sHHx/10,'Color',color,'DisplayName','HindHeelX/10')
        end
        color=colors(6,:);%rand(1,3);%y is inverted because origin is top left
      sHHy=sHHy-((floorfit(1)*(1:length(sHHy))+floorfit(2)))';
      [~,output(i).HeelProminance.avg,~,output(i).HeelProminance.std] = pulloutYstuff(sHHy,sHHf,spotcheck,color,'HindHeelY');%mean 1/2 peak width and mean prominance
      [output(i).HH.y.avg, ~,output(i).HH.y.std] = avgforstride (sHHy,sHHf, strides, [],spotcheck,color,'HindHeelY');
        if spotcheck && pausevalue
            pause
        end
%% checking. 
% gives you the option to pause, view the graph, and put in an evaluation.

       if spotcheck && false%check here if you want to evaluate each graph as it pops up
           output(i).byEye = input('evaluation. examples include g for good, e for need to edit, b for bad\n','s');
       elseif output(i).byEye==""
           output(i).byEye = 'not evaluated'; 
       else
       end
% forefoot

    Fx = ForeX(ForeL>L);Ff=frame(ForeL>L);color=colors(7,:);%rand(1,3);
    Fy= -ForeY(ForeL>L);Fy=Fy-((floorfit(1)*(1:length(Fy))+floorfit(2)))';
    [sFx,sFy,sFf] = basiccmooth(Fx, Fy, Ff, d2);%using d2 filter    
      [output(i).F.y.avg, ~,output(i).F.y.std] = avgforstride (sFy,sFf, strides, [],spotcheck,color,'HindHeelX');
      [output(i).F.x.avg, output(i).F.f.avg,output(i).F.x.std,output(i).F.f.std] = avgforstride(sFx,sFf, strides, stridesXstart,spotcheck,color,'ForeX');
     color=colors(8,:);%rand(1,3);
      [output(i).F.y.avg, ~,output(i).F.y.std] = avgforstride(sFy,sFf, strides, [],spotcheck,color,'ForeY'); 
        if spotcheck
            figure(1);hold on; 
            plot(sFf,sFy,'Color',color,'DisplayName','ForeY');
            if pausevalue
            pause
            end
        end
% tailbase to midfoot distance

    TBMF = sqrt( ((HindMidY-TailBasey).^2)+((HindMidX-TailBasex).^2) ); color=colors(9,:);%rand(1,3);%Tailbase to midfoot length
        TBMF =  TBMF(HindMidL>L & TailBasel>L);        TBMFf =frame(HindMidL>L & TailBasel>L);
        [sTBMF,~,sTBMFf] = basiccmooth(TBMF, TBMF,TBMFf, d1);%using og filter   
      [output(i).TBMF.y.avg, output(i).TBMF.f.avg,output(i).TBMF.y.std,output(i).TBMF.f.std] = avgforstride(sTBMF,sTBMFf, strides, [], spotcheck,color,'TailBaseMidFootLength');
        if spotcheck
            figure(1);hold on;
            %plot(TBMFf,TBMF,'DisplayName','TailBasetoMidFootLength');
            plot(sTBMFf,sTBMF/10,'Color',color,'DisplayName','TailBasetoMidFootLength/10')
            if pausevalue
            pause
            end
        end
% HT-HM-HH
% the angle at the midfoot from the hind toe and hind heel. Shows how flexed 
% the foot is. Tends to be lower range in sciatica.

        color=colors(10,:);%rand(1,3);
      [angle, goodframes] = findangle(HindToeX,HindToeY,HindToeL,HindMidX,HindMidY,HindMidL, HindHeelX,HindHeelY,HindHeelL,frame,L,spotcheck,color,'toe-midfoot-heel angle'); %2nd is where the angle is
        [output(i).HT_HM_HH.y.avg, output(i).HT_HM_HH.f.avg,output(i).HT_HM_HH.y.std,output(i).HT_HM_HH.f.std] = avgforstride(angle',goodframes, strides, [],spotcheck,color,'toe-midfoot-heel angle');
        if spotcheck && pausevalue
            pause
        end
% foot angle at toeoff
% FSTOs is in frames. we should get an average, so let's do it over half a second, 
% at 500frames/s. This window size can be edited by editing the dropdowns.

        toeoffs1 = FSTOs(:,2); %this is in frames, need to find index that matches
        toeoffs  = NaN(size(toeoffs1(toeoffs1~=0)));
        for j=1:length(toeoffs1(toeoffs1~=0))
            [~,toeoffs(j)] = min(abs(goodframes-toeoffs1(j)));
        end
        FootAngleToeOff=NaN(size(toeoffs));
        for k=1:length(toeoffs)
        winstart = toeoffs(k)-125;winstart(winstart<1)=1; 
        winend   = toeoffs(k)+ 125;winend(winend>length(angle))=length(angle);           
        FootAngleToeOff(k) = median(angle(winstart:winend));%in degrees
        end
        output(i).FootAngleToeOff.avg = mean(FootAngleToeOff);
        output(i).FootAngleToeOff.std = std(FootAngleToeOff);
% HT-HM-TB
% angle between hind toe, hind midfoot, and tailbase (at the midfoot). Not sure 
% how useful this is tbh.

        color=colors(11,:);%rand(1,3);
      [angle, goodframes] = findangle(HindToeX,HindToeY,HindToeL,HindMidX,HindMidY,HindMidL, TailBasex,TailBasey,TailBasel,frame,L,spotcheck,color,'toe-midfoot-tailbase angle'); %2nd is where the angle is
        [output(i).HT_HM_TB.y.avg, output(i).HT_HM_TB.f.avg,output(i).HT_HM_TB.y.std,output(i).HT_HM_TB.f.std] = avgforstride(angle',goodframes, strides, [],spotcheck,color,'toe-midfoot-tailbase angle');
        if spotcheck && pausevalue
            pause
        end
% HT-TB-B
% angle between hind toe, tailbase, and back (at the tailbase). Not sure how 
% useful this is tbh.

        color=colors(12,:);%rand(1,3);
      [angle, goodframes] = findangle(HindMidX,HindMidY,HindMidL, TailBasex,TailBasey,TailBasel,Backx, Backy, Backl, frame,L,spotcheck,color,'toe-tailbase-back angle'); %2nd is where the angle is
        [output(i).HT_TB_B.y.avg, output(i).HT_TB_B.f.avg,output(i).HT_TB_B.y.std,output(i).HT_TB_B.f.std] = avgforstride(angle',goodframes, strides, [],spotcheck,color,'toe-tailbase-back angle');
        if spotcheck && pausevalue
            pause
        end
% angle of virtual limb (toe to tailbase)

        color=colors(12,:);%rand(1,3);
      [angle, goodframes] = findangle(HindToeX,HindToeY,HindToeL, TailBasex,TailBasey,TailBasel,TailBasex, HindToeY, HindToeL, frame,L,spotcheck,color,'toe-tailbase angle'); %2nd is where the angle is
        [output(i).HT_TB.y.avg, output(i).HT_TB.f.avg,output(i).HT_TB.y.std,output(i).HT_TB.f.std] = avgforstride(angle',goodframes, strides, [],spotcheck,color,'toe-tailbase angle');
        if spotcheck && pausevalue
            pause
        end
% angle of virtual limb (toe to back)

        color=colors(12,:);%rand(1,3);
      [angle, goodframes] = findangle(HindToeX,HindToeY,HindToeL, Backx,Backy,Backl,Backx, HindToeY, HindToeL, frame,L,spotcheck,color,'toe-back angle'); %2nd is where the angle is
        [output(i).HT_B.y.avg, output(i).HT_B.f.avg,output(i).HT_B.y.std,output(i).HT_B.f.std] = avgforstride(angle',goodframes, strides, [],spotcheck,color,'toe-back angle');
        if spotcheck && pausevalue
            pause
        end
%% nonsided
% measures that aren't the left or right side only
% nose x
% useful for judjing if the animal has stopped, but can also be used for velocity

    Nf = frame(nosel>L);color=colors(13,:);%rand(1,3);%Nx was defined earlier
    Ny= -nosey(nosel>L);Ny=Ny-((floorfit(1)*(1:length(Ny))+floorfit(2)))';
       [sNx,sNy, sNf] = basiccmooth(Nx, Ny,Nf, d1);%using og filter
     [output(i).N.x.avg, output(i).N.f.avg,output(i).N.x.std,output(i).N.f.std] = avgforstride(sNx,sNf, strides, stridesXstart,spotcheck,color,'NoseX');
       if spotcheck
            figure(1);hold on;
            plot(sNf,sNx/10,'Color',color,'DisplayName','NoseX/10')
       end
     color=colors(14,:);%rand(1,3);%height of nose relative to floor
      [output(i).N.y.avg, ~,output(i).N.y.std,~] = avgforstride(sNy,sNf, strides, [],spotcheck,color,'NoseY'); 
        if spotcheck
            figure(1);hold on;
            plot(sNf,sNy,'Color',color,'DisplayName','NoseY')
            if pausevalue
            pause
            end
        end
% back x
% useful for velocity

    Bx= Backx(Backl>L); Bf=frame(Backl>L); color=colors(15,:);%rand(1,3);
    By= -Backy(Backl>L); By=By-((floorfit(1)*(1:length(By))+floorfit(2)))';
        [sBx,sBy, sBf] = basiccmooth(Bx, By,Bf, d1);%using og filter
       [output(i).B.x.avg, output(i).B.f.avg,output(i).B.x.std,output(i).B.f.std] = avgforstride(sBx,sBf, strides, stridesXstart,spotcheck,color,'BackX');
        if spotcheck
            figure(1);hold on;
            plot(sBf,sBx/10, 'Color',color,'DisplayName','BackX/10')
        end
    color=colors(16,:);%rand(1,3);%height of back relative to floor already defined
       [output(i).B.y.avg, ~,output(i).B.y.std,~] = avgforstride(sBy,sBf, strides, [],spotcheck,color,'BackY'); 
        if spotcheck
            figure(1);hold on;
            plot(sBf,sBy/2, 'Color',color,'DisplayName','BackY/2')
            if pausevalue
            pause
            end
        end
% backslope
% we thought animals with sciatica might be more hunched over

    backslope = (-(Backy-TailBasey)./(Backx-TailBasex));  color=colors(17,:);%rand(1,3);%(y-y/x-x) %y is negative because the origin is in the top left
        backslope = backslope(Backl>L & TailBasel>L); bsf =frame(Backl>L & TailBasel>L);
        [sbackslope,~,sbsf] = basiccmooth(backslope,backslope, bsf, d1);
        backslope = -directionfactor*(sbackslope);% needed.
        backangle = atand(sbackslope);
      [nBackA, nBAf,nBAstd] = avgforstride(backangle,sbsf, strides, [],spotcheck,color,'BackAngle'); 
        if spotcheck == 1
            figure(1); hold on;
            %plot(bsf,-directionfactor*backslope*10,'DisplayName','backslope*10');
            plot(sbsf,backangle,'Color',color,'DisplayName','backangle');
            if pausevalue
            pause
            end
        end
% rough velocity
% this isn't really correct velocity because you haven't added the pixels per 
% cm for each trial. just use as a reference.

     velocity = (Nx(end)-Nx(1))/pixeltocm; %distance in pixels*cm/pix = cm
     velocity = velocity/((Nf(end)-Nf(1))/framestosecs); %cm/ frames*s/frames=s)->cm/s
     output(i).RoughVelocity = abs(velocity);
% stepwidth

     %smooth the mirror feet locations
     MRHy= -MirrorRightHindy(MirrorRightHindl>L); MLHy=-MirrorLeftHindy(MirrorLeftHindl>L);
     MRHf=  frame(MirrorRightHindl>L);MLHf=  frame(MirrorLeftHindl>L);
     [sMRHy,~, sMRHf] = basiccmooth(MRHy,MRHy,MRHf, d1);%using og filter
     [sMLHy,~, sMLHf] = basiccmooth(MLHy,MLHy,MLHf, d1);%using og filterdelength.avg = mean(abs(MirrorRightHindy(strideIndex(:,1))-MirrorRightHindy(strideIndex(:,2))));
     %so basically the mirrored feet are mostly in lines with random little
     %peaks. So we find the locations of those peaks, then find the value
     %of the y position halfway between them, and use that to find the
     %width.
     [~,locsR] = findpeaks(-sMRHy,sMRHf,'MinPeakProminence',4,'MinPeakDistance',25,'Annotate','extents');
     [~,locsL] = findpeaks(-sMLHy,sMLHf,'MinPeakProminence',4,'MinPeakDistance',25,'Annotate','extents');
     %locsR=[sMRHf(1);locsR];
     %locsL=[sMLHf(1);locsL];
     locs=min([length(locsR),length(locsL)]);
     selectedR=NaN(locs-1,1);selectedL=NaN(locs-1,1);
     for n=1:locs-1
        [~,interum] = min(abs(sMRHf-mean([locsR(n),locsR(n+1)])));%the index of the frame halfway between the peaks 
        selectedR(n) = sMRHy(interum);%the y value halfway between the peaks
        [~,interum] = min(abs(sMLHf-mean([locsL(n),locsL(n+1)])));
        selectedL(n) = sMLHy(interum);
     end   
     output(i).stepwidth.avg         = mean(abs(selectedR-selectedL));
     output(i).stepwidth.variability  = std(abs(selectedR-selectedL));
 
%% 
% %trial is done, saving and iterating

    %pause(3) 
    if spotcheck %saves the graph images to whatever your active directory is.
     saveas(fig1,[myfile,'_fulltrial.png']);
     saveas(fig2,[myfile,'_avgstride.png']);
    end
catch ME
    disp('this trial had an error');
    disp(['error: ',ME.message, ' ', ME.stack(1).name, ' line ', num2str(ME.stack(1).line)]);
    output(i).byEye=['error: ',ME.message, ' ', ME.stack(1).name, ' line ', num2str(ME.stack(1).line)];
end %ends the 'try' statement
end
outputFileName = ['output',num2str(now),'.mat'];
save(outputFileName, 'output')%saves the output structure without overwrighting because the 'now' function gets the current date and time as a serial date number
allDone = 'Yay!'
%% Functions
%% finding an angle
% given the x,y,and likelyhood values for 3 points over time
% 
% it finds the angle at point 2

function [angle3, goodframes] = findangle(x1,y1,l1, x2, y2, l2, x3, y3, l3,frame, L,spotcheck,color,name)
     d1 = designfilt('lowpassiir','FilterOrder',12, ...
     'HalfPowerFrequency',0.05,'DesignMethod','butter');

     angle = [x1,y1, x2, y2, x3, y3];
     angle = angle((l1>L & l2>L & l3>L),:);
     goodframes = frame((l1>L & l2>L & l3>L),:);
     P0 = angle(:,3:4);
     P1 = angle(:,1:2);
     P2 = angle(:,5:6);
     angle3 = NaN(length(goodframes),1);
     for i=1:length(goodframes)
        n1 = (P2(i,:) - P0(i,:)) / norm(P2(i,:) - P0(i,:));  % Normalized vectors
        n2 = (P1(i,:) - P0(i,:)) / norm(P1(i,:) - P0(i,:));
        angle3(i) = atan2(norm(det([n2; n1])), dot(n1, n2));  % Stable
     end
     try
        angle3=filtfilt(d1, angle3);
     catch
         disp('issue smoothing the angle vector');
     end
     angle3=angle3*180/pi; %convert from radians to degrees
     if spotcheck==1
         figure(1); hold on;
         plot(goodframes,angle3,'Color',color,'DisplayName',[name,' degrees']);
     end
end
%% finding a floor when it slopes

function [P] = floorfind(HTy,HMy,HHy)
datas = {HTy,HMy,HHy};
floory=NaN(3,3);floorx=NaN(3,3);
for i =1:3
    data = datas{i};
    [floory(1,i),floorx(1,i)]=min(data(1:end/3));
    [floory(2,i),floorx(2,i)]=min(data(end/3+1:2*end/3));
    [floory(3,i),floorx(3,i)]=min(data(2*end/3+1:end));
    floorx(2,i)=floorx(2,i)+length(data(1:end/3));
    floorx(3,i)=floorx(3,i)+length(data(1:2*end/3));
end
floorx=reshape(floorx,9,1);
floory=reshape(floory,9,1);
P = polyfit(floorx,floory,1);
%% 
% for reference, floorvec = P(1)*(1:length(data))+P(2);

end
%% finding the prominance of the height value
% finds locations and prominance of the max height for each stride/peak.

function [meanwidth,meanprominance,stdwidth,stdprominance] = pulloutYstuff(data,frames,spotcheck,color,name)
sdata=data;
[pks,locs,w,p] = findpeaks(sdata,frames, 'MinPeakProminence',4,'MinPeakDistance',25,'Annotate','extents');

%Getting rid of ouliers
[~,TFp] = rmoutliers(p,'mean', 'ThresholdFactor',2); %gets rid of anything more than 2 SDs from mean prominance
[~,TFw] = rmoutliers(w,'mean', 'ThresholdFactor',2); %gets rid of anything more than 2 SDs from mean 1/2 peak width
TF = TFp+TFw;
 if spotcheck == 1
    figure(1)
    hold on
    %plot(frames,data,'DisplayName',name);
    plot(frames,sdata,'Color',color,'DisplayName',[name,' smooth']);
    plot(locs(TF==0),pks(TF==0),'r*','DisplayName',[name,' peaks']);
    %refline(0, 0);
 end
meanprominance = mean(p(TF==0));stdprominance=std(p(TF==0));
meanwidth = mean(w(TF==0));stdwidth=std(w(TF==0));
end
%% fsto segmenting with x values

function [dutyfactor, stridelength, numcycles, strides, dutyfactorStd,stridelengthStd, byEye,fstos]=fsto(data,frames, directionfactor,spotcheck, graphfactor,color, name)
byEye="";
%% 
% find derivative and then smooth it.

dxdf = gradient(data(:));% ./ gradient(frames(:));
[dxdf] = rmoutliers([dxdf,frames],'movmedian',25);
rframes=dxdf(:,2);
dxdf=dxdf(:,1);
%% 
% find a good cutoff for the derivative
% 
% find peaks

[pks,~,~,p] = findpeaks(-directionfactor*(dxdf), 'MinPeakProminence',1,'MinPeakDistance',25,'Annotate','extents');
%% 
% The threshold for determining stance is determined by finding when the derivative 
% of the x coordinate (or -derivative depending on direction) drops below a certain 
% value. This thresholdis defined as the fraction (between 0 and 1) of the prominance 
% of the peak derivative in swing of that stride. A smaller value means a lower 
% threshold aka the foot has to be 'more still'. Default is 0.15. At 0, we tend 
% to get noise that makes the program think the foot is in swing when it's just 
% a little jitter.

threshs=pks-(1-0.15)*p; %tunable
threshi=NaN(length(threshs),1);
for i=1:length(threshs)%find the x values
    threshi(i,1)=find(-directionfactor*(dxdf)==pks(i));
end
if false%Check if your strides are very regular. I would guess this is typically checked in OA animals.
    P = polyfit(threshi,median(threshs)*ones(size(threshi)),1); %just fit to median - use this if you need to edit it because of outliers
    byEye = append(byEye,"stance/swing threshold is median");
else
    P = polyfit(threshi,threshs,1); %fit a line through the threshold values to make a little threshold line
end
threshline=((P(1)*(1:length(dxdf))+P(2)))';
%threshline=median(threshs)*ones(size(dxdf));
%% 
% define stride/stance based off this threshold to create a binary vector with 
% stance as 1 and swing as 0

stance1=NaN(size(dxdf));
stance1(-directionfactor*(dxdf)<=threshline)=1;
stance1(-directionfactor*(dxdf)>threshline)=0;%stance is 1
%% 
% We want to filter this binary vector to get rid of unintended stance or swing 
% points. The filter size should be around the size of stance, but if it is too 
% big it will make stance smaller, so make it the smallest reasonably expected 
% stance time (in indecies)

[stance1] = movmedian(stance1,60);%in frames - tunable
stance1(stance1==.5)=1; %the movmedian can result in .5 values when it goes from swing to stance or vice versa. I am assigning those to be stance because it kind of looks like they should be. Tunable, but only 1 frame
% find footstrike
% find points where the 'stance' binary vector changes (with certain limitations)

[pks,toestrike]= findpeaks([stance1;0],[rframes;rframes(end)+1],'MinPeakDistance',25,'MinPeakWidth',35);% frames in stance - tunable
%% 
% I added a point at 0 at the end because if I don't it doesn't see the last 
% stance time as a peak if it isn't followed by a trough.

%find toeoff
[~,toeoff]= findpeaks(-[stance1;0],[rframes;rframes(end)+1],'MinPeakDistance',25,'MinPeakWidth',35);% frames in swing - tunable
%% 
% %abort the trial if there are not footstrikes

if length(toestrike)<=1 
    dutyfactor=nan; stridelength=nan;strides=nan;
    numcycles=nan; 
    dutyfactorStd=nan;stridelengthStd=nan; byEye = 'toestrike<=1';
    return
end
% set up strides matrix

strides = [toestrike(1:end-1),toestrike(2:end)];
stridelengthF= strides(:,2)-strides(:,1);
%% 
% remove weird strides - determines which strides you will count in average. 
% decreases the stride count but can automatically get rid of strides with different 
% stridelengths relative to others in the trial. 

if false%Check here if you want to include all strides, even if the stridelength is an outlier (as judged by matlab's rmoutliers function)
else
    [stridelengthF,stridelengthTF] = rmoutliers(stridelengthF,'median'); %remove outliers in sliding window
    strides=strides(~stridelengthTF,:);
end

% find dutyfactor

fstos(:,1)=unique(strides);
dutyfactor=NaN(size(strides,1),1);
stridelength=NaN(size(strides,1),1);
ToeOff=NaN(size(strides,1),1);
for j=1:size(strides,1)
    stridestart=strides(j,1); %in frames
    strideend  =strides(j,2);
    strideToeOff=toeoff(toeoff>stridestart); %find toeoff after the toestrike
    ToeOff(j) = strideToeOff(1); %we just want the first one after toestrike
    fstos(j,2)=ToeOff(j);
    %temporal
    timestance =ToeOff(j)-stridestart; %in frames
    dutyfactor(j) = timestance/(strideend-stridestart); %duty factor
    %spatial
    stridelength(j) =abs(data(frames==stridestart)-data(frames==strideend)); 
end
% plotting

val = 0;%conditional value set up for editing (if you want to edit the FSTO values as you run it)
if spotcheck
    offsetforgraphing=rand(1);%I offset the stride/stance so that I can look at multiple stride/stance at once and they don't overlap
    while val ==0
        figure(1);hold on;
        %plot(frames,   data/10, 'Color',color,'DisplayName',['HT/10']);
        plot(rframes, dxdf*10, 'Color',color,'DisplayName',['|derivative*10|',name]);
        plot(rframes,threshline*10,'DisplayName','threshold line for stance/swing*10')
        plot(rframes, (stance1-.1)*100+offsetforgraphing,'.','Color',color*(2.5/3),'DisplayName',[name,' stance']);
        plot([toestrike,toestrike]',[(pks-1.1)*100-offsetforgraphing,(pks-.1)+100+offsetforgraphing]','Color',color*(2.2/3),'HandleVisibility','off'); 
            plot(0,0,'Color',color,'DisplayName',[name,' toestrike']); %just for the legend info
        plot(strides',[(ones(size(stridelengthF))-1.1)*100-offsetforgraphing*2,ones(size(stridelengthF))*100+offsetforgraphing*2]','Color',color*(2.2/3),'HandleVisibility','off');%toeoff 
        plot([ToeOff';ToeOff'],[(ones(length(ToeOff),1)-1.1)*100-offsetforgraphing,(ones(length(ToeOff),1)-.1)+50+offsetforgraphing]','Color',color*(2.2/3),'HandleVisibility','off');%toeoff 

%% 
% modify FSTOs for problem trials. Input will be requested in the command line.

        if false%check this box if you want to modify stride FSTOs for each trial as it runs.
            byEye = input('evaluation. enter e to edit.\n','s');
        end
        if byEye == 'e'
            disp(['strides =' [toestrike(1:end-1),toestrike(2:end)]])
            stridesToUse = input('choose stride indicies, use [] to contain them');
            strides = strides(stridesToUse, :);
            stridelengthF= strides(:,2)-strides(:,1);
            offsetforgraphing = offsetforgraphing+1;
        else 
            val =1;
        end
    end
end
%% 
% %define some things and average

numcycles = size(strides,1);
dutyfactorStd   = std(dutyfactor);    dutyfactor=mean(dutyfactor);    
stridelengthStd = std(stridelength);stridelength=mean(stridelength);
%% 
% %graph the average footstrike/toeoff, but we won't use the outputs

[~, ~, ~,~] = avgforstride(stance1*100,frames,strides, [], spotcheck,color*(1/3),[name,'Stance1 Swing0']);
end
%% Average strides together
% This function bins each of the strides into windows, then finds the stride 
% with the lowest number of points (minimum 1) within that window, and then decreses 
% the amount of points in every other stride to that value in that window. If 
% there is no data in that window, it tries to get 1 data point based on prior 
% and subsequent windows (sometimes happens because of likelyhood cutoffs).

function [avnwinstrides, avnwinfs, avwinstridesSTDEV,avwinfsSTDEV] = avgforstride(data,f2, strides, stridesXstart, spotcheck,color,name)
    winsize = 10;%size in percentage of stride. increase if you have faster strides.
% plot as percent of stride

    for j=1:size(strides,1)
        stridestartf=strides(j,1);
        strideendf  =strides(j,2);
        stridef = f2(f2>=stridestartf & f2<strideendf);
        if length(stridef)<1
            percentf=NaN; winfs{j}=nan;
            mystride=NaN; winstrides{j} =nan;
            break
        end        
        stridef = stridef - stridef(1);%need frames from begining of stride
        percentf = stridef/max(stridef)*100;%convert frames to % of all frames
        mystride = data(f2>=stridestartf& f2<strideendf);%data for stride j
        if length(stridesXstart)>1
        mystride = mystride-stridesXstart(j);%if this is x data, i want to normalize to the begining of the stride so that's the minimum
        end
        
        if spotcheck
            figure(2); hold on;
            %plot(percentf, mystride,':','Color',color/j,'DisplayName',[name,' stride ',num2str(j)])
            plot(percentf, mystride,':','Color',color,'HandleVisibility','off') %plot individual strides
        end
        winnum=1;%index for which window bin we are on
        for k = 0:winsize:100-winsize %window moving over data to get only minnumpoints for each winsize % of stride
            winstrides{j,winnum} = mystride(percentf>=k & percentf<k+winsize);%save strides
            winfs{j,winnum} = percentf(percentf>=k & percentf<k+winsize);%save f values for strides
            winnum=winnum+1;
        end
    end
    %j
%     if exist('winstrides', 'var')==0
%         pause
%     end
    % msg=['testing sizes line 420, should be true: ',num2str(size(winstrides,2)==100/winsize)] %troubleshooting
%% 
% %initialize some things

    nwinstrides=cell(size(winstrides,1),1);
    nwinfs=cell(size(winstrides,1),1);
    weights=cell(size(winstrides,1),1);
    
%% 
% %loop thru the number of bins that exist. should = 100/winsize

    for m=1:size(winstrides,2) %loop thru the number of bins that exist. should = 100/winsize
        wins=winstrides(:,m);%select only bin m for all strides
%% 
% %find the minimum number of points/stride in the window

        minps = min(cellfun('size',wins,1)); %find the minimum number of points/stride in the window
        if minps<=0
            disp(['need to increase winsize? winsize=',num2str(winsize)])
            %or maybe exclude the stride missing the data?
            %could decide based on stride length outliers?
            minps=1;
        end
        
%% 
% % loop thru the strides (n) for each bin (m)

        for n=1:size(winstrides,1) %loop thru the strides (n) for each bin (m).
            winnm=winstrides{n,m}; %select the bin in question
            winfnm=winfs{n,m}; %select frame % for the bin in question
            
%% 
% %if the stride has more points than the minimum, decrease to the minimum.

            if length(winnm)>minps
                newis = round(linspace(1,length(winnm),minps));%new indecies
                nwinnm = winnm (newis); %new window for this stride and bin
                nwinfnm = winfnm(newis); %new percentage window for this stride and bin
                nwinstrides{n} = [nwinstrides{n}, nwinnm'];%add this window to a large vector of windows
                nwinfs{n}      = [     nwinfs{n},nwinfnm'];%
                weights{n}     = [    weights{n},ones(size(nwinnm'))];%weight them all evenly                
%% 
% %if the stride has the minimum number of points, keep it the same.

            elseif length(winnm)==minps
                nwinstrides{n} = [nwinstrides{n}, winnm'];%
                nwinfs{n}      = [     nwinfs{n},winfnm'];%
                weights{n}     = [    weights{n},ones(size(winnm'))];%weight them all evenly  
%% 
% %if a stride has no data in a window, try to interpolate 

            else %if i don't have enough data i'm interpolating.
                disp(['no data for ',name,' Toeoff ',num2str(n),', percentf ',num2str(m)]);
                %initialize
                winmminus1 = nan; winmplus1  = nan; winfnminus1 = nan; winfnplus1  = nan;winnminus1=nan;winfallns=nan;
                %find data for previos and next bin or what we just
                %appended
                try    winmminus1 = winstrides{n,m-1};  end %same stride bin previous
                try    winmplus1  = winstrides{n,m+1};  end %same stride bin next
                try    winnminus1 = nwinstrides{n}(end);  end %the avg bin we just appended (this works, i just think it's better to nan these)
                for smalln=1:size(winfs,1) %all stride f for bin
                    winfallns = [winfallns;(winfs{smalln,m})];      
                end 
                %find the mean of our estimations, the interpolated value.
                nwinnm = nanmean([ winmminus1; winmplus1;winnminus1]);%what if mminus one and plus 1 are both empty??
                nwinfnm= nanmean(winfallns);
                
                %if isnan(mean([nwinnm,nwinfnm])) %catch error
                    % error('we do not got any data breh')
                %end
                
                nwinstrides{n} = [nwinstrides{n}, nwinnm];%
                nwinfs{n}      = [     nwinfs{n},nwinfnm];%
                weights{n}     = [    weights{n},0.5*ones(size(nwinnm))];%weigh these points less when averaging strides together
                %alternative - make winsize higher? or include 2 bins for this point?
            end
            if sum([size(nwinstrides{n})==size(nwinfs{n}),size(nwinstrides{n})==size(weights{n}),-4])
            error('yo our stuff is not lining up right')
            end
        end
    end
%% 
% %average strides together

    nwinstrides = cell2mat(nwinstrides); nwinfs = cell2mat(nwinfs); weights = cell2mat(weights); 
    % nanmean with weights for data
    avnwinstrides = nansum(nwinstrides.*weights,1)./nansum(weights,1);
    avwinstridesSTDEV = sqrt(nansum(weights.*((nwinstrides-avnwinstrides).^2))./(nansum(weights,1)-(nansum(weights.^2,1)/nansum(weights))));
    % nanmean with weights for percentf
    avnwinfs = nansum(nwinfs.*weights,1)./nansum(weights,1);
    avwinfsSTDEV =sqrt(nansum(weights.*((nwinfs-avnwinfs).^2))./(nansum(weights,1)-(nansum(weights.^2,1)/nansum(weights))));

    if spotcheck
        figure(2);hold on;
        plot(avnwinfs', avnwinstrides','Color',color,'LineWidth',2,'DisplayName',[name,' avg'])
    end

end
%% FSTO Y

function [dutyfactor, stridelength, numcycles, strides, dutyfactorStd,stridelengthStd, byEye,fstos]=fstoY(data,frames, directionfactor,spotcheck, graphfactor,color, name)
byEye="";
%% 
% find peaks - where foot is raised

[pks,locs,w,p] = findpeaks(data, 'MinPeakProminence',4,'MinPeakDistance',25,'Annotate','extents');
%% 
% define threshold below which foot is in stance. This thresholdis defined as 
% the fraction (between 0 and 1) of the prominance of the peak height of that 
% stride. A smaller value means a lower threshold aka the foot has to be lower. 
% Default is 0.25. At 0, we tend to get noise that makes the program think the 
% foot is in swing when it's just a little jitter.

threshs=pks-(1-0.25)*p; locs=[locs;length(data)]; 
toestrike=NaN(size(threshs));threshs=[threshs;data(end)];
stance=NaN(size(data));
for i =1:length(locs)-1
    window=data(locs(i):locs(i+1));%window from one peak in foot height to the next
    
    swindow = find(window<threshs(i));%strike window where I think the foot is below a certain height
    swindow = swindow(1)+locs(i);%add back the begining index so it's in the context of the full trial, not just window
    toestrike(i) = frames(swindow);%the first frame below my threshold value should be toestrike.
    
     owindow=find((window)<threshs(i+1));%segment where it is lower than the threshold for the next peak - used to judge toeoff
     if length(owindow)<1
         toestrike(i)=nan;
         continue
     end
     owindow = owindow(end)+locs(i);
     toeoff(i)=frames(owindow);
     stance(locs(i):locs(i+1))=0;%set to 0
     stance(swindow:owindow)=1;%set to 1 where I think it's in stance
end
toestrike = toestrike(~isnan(toestrike));
pks=ones(size(toestrike));
%% 
% %abort the trial if there are not footstrikes

if length(toestrike)<=1
    dutyfactor=nan; stridelength=nan;strides=nan;
    numcycles=nan; 
    dutyfactorStd=nan;stridelengthStd=nan;
    return
end
%% 
% %set up strides matrix

strides = [toestrike(1:end-1),toestrike(2:end)];
stridelengthF= strides(:,2)-strides(:,1);
%% 
% remove weird strides - determines which strides you will count in average. 
% decreases the stride count but can automatically get rid of strides with different 
% stridelengths relative to others in the trial. 

if false%Check here if you want to include all strides, even if the stridelength is an outlier (as judged by matlab's rmoutliers function)
else
    [stridelengthF,stridelengthTF] = rmoutliers(stridelengthF,'median'); %remove outliers in sliding window
    strides=strides(~stridelengthTF,:);
end
%% 
% %find dutyfactor

fstos(:,1)=unique(strides);
dutyfactor=NaN(size(strides,1));
stridelength=NaN(size(strides,1));
ToeOff=NaN(size(strides,1));
for j=1:size(strides,1)
    stridestart=strides(j,1);
    strideend  =strides(j,2);
    strideToeOff=toeoff(toeoff>stridestart); %find toeoff after the toestrike
    ToeOff(j) = strideToeOff(1); %we just want the first one after toestrike
    fstos(j,2)=ToeOff(j);
    %temporal
    timestance =ToeOff(j)-stridestart; %in frames
    dutyfactor(j) = timestance/(strideend-stridestart); %duty factor
    %spatial - can't find this because we have Y values, not X values
    %stridelength(j) =abs(data(frames==stridestart)-data(frames==strideend)); 
end
% plotting

if spotcheck
    val = 0;%conditional value set up for editing (if you want to edit the FSTO values as you run it)
    offsetforgraphing=rand(1);%I offset the stride/stance so that I can look at multiple stride/stance at once and they don't overlap
    while val ==0
        figure(1);hold on;
        plot(frames,   data/10, 'Color',color,'DisplayName',['HT/10']);
        plot(frames, (stance-.1)*10+offsetforgraphing,'.','Color',color*(2.5/3),'DisplayName',[name,' stance']);
        plot([toestrike,toestrike]',[(pks-1.1)*100-offsetforgraphing,(pks-.1)+100+offsetforgraphing]','Color',color*(2.2/3),'HandleVisibility','off'); 
            plot(0,0,'Color',color,'DisplayName',[name,' toestrike']); %just for the legend info
        plot(strides',[(ones(size(stridelengthF))-1.1)*100-offsetforgraphing*2,ones(size(stridelengthF))*100+offsetforgraphing*2]','Color',color*(2.2/3),'HandleVisibility','off');%toeoff 
        plot([ToeOff;ToeOff],[(ones(length(ToeOff),1)-1.1)*100-offsetforgraphing,(ones(length(ToeOff),1)-.1)+50+offsetforgraphing]','Color',color*(2.2/3),'HandleVisibility','off');%toeoff 

%% 
% modify FSTOs for problem trials. Input will be requested in the command line.

        if false%check this box if you want to modify stride FSTOs for each trial as it runs.
            byEye = input('evaluation. enter e to edit.\n','s');
        end
        if byEye == 'e'
            disp(['strides =' [toestrike(1:end-1),toestrike(2:end)]])
            stridesToUse = input('choose stride indicies, use [] to contain them');
            strides = strides(stridesToUse, :);
            stridelengthF= strides(:,2)-strides(:,1);
            offsetforgraphing = offsetforgraphing+1;
        else 
            val =1;
        end
    end
end
%% 
% %define some things and average

numcycles = size(strides,1);
dutyfactorStd   = std(dutyfactor);    dutyfactor=mean(dutyfactor);    
stridelengthStd = std(stridelength);stridelength=mean(stridelength);
%% 
% %graph the average footstrike/toeoff, but we won't use the outputs:

[~, ~, ~,~] = avgforstride(stance*100,frames,strides, 'y', spotcheck,color*(1/3),[name,'Stance1 Swing0']);
end
%% a basic smoothing function
% filters the data so we have smoother lines, but also removes frame values 
% for the x/y data points that were removed.

function [sdatax,sdatay,sdataf] = basiccmooth(datax,datay,dataf, d1)
    if exist('d1','var')
        %do nothing
    else %if i do not pass it a filter (d1), I make one here
        d1 = designfilt('lowpassiir','FilterOrder',3, ...
        'HalfPowerFrequency',0.1,'DesignMethod','butter');
    end
    [sdata,TFdata] = rmoutliers([datax,datay],'movmedian',20);
    sdatax=sdata(:,1);
    sdatay=sdata(:,2);
    sdataf=  dataf(sum(TFdata,2)==0);
    try
    sdatax = filtfilt(d1, sdatax); %rmoutliers isn't a 0 shift filter
    sdatay = filtfilt(d1, sdatay); %rmoutliers isn't a 0 shift filter
    catch
        disp('issue filtering x or y data of a point')
    end
end
%% import function
% assumes the data is in the format output by DLC trained with these specific 
% points in this order!

function [frame, nosex,nosey,nosel,RightForex,RightForey,RightForel,LeftForex,LeftForey,LeftForel,RightHindToex,RightHindToey,RightHindToel,RightHindMidx, RightHindMidy, RightHindMidl, RightHindHeelx,RightHindHeely,RightHindHeell,LeftHindToex,LeftHindToey,LeftHindToel,LeftHindMidx, LeftHindMidy, LeftHindMidl, LeftHindHeelx,LeftHindHeely,LeftHindHeell,MirrorRightHindx,MirrorRightHindy,MirrorRightHindl,MirrorLeftHindx,MirrorLeftHindy,MirrorLeftHindl,Backx,Backy,Backl,TailBasex,TailBasey,TailBasel,ABforcePlateCenterx,ABforcePlateCentery,ABforcePlateCenterl,CDforcePlateCenterx,CDforcePlateCentery,CDforcePlateCenterl] = sciimport(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as column vectors.
%   [SCORER,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_180000,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_1,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_2,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_3,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_4,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_5,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_6,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_7,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_8,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_9,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_10,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_11,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_12,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_13,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_14,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_15,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_16,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_17,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_18,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_19,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_20,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_21,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_22,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_23,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_24,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_25,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_26]
%   = IMPORTFILE(FILENAME) Reads data from text file FILENAME for the
%   default selection.
%
%   [SCORER,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_180000,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_1,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_2,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_3,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_4,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_5,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_6,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_7,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_8,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_9,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_10,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_11,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_12,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_13,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_14,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_15,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_16,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_17,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_18,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_19,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_20,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_21,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_22,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_23,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_24,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_25,DLC_RESNET50_SCIATICAMAR25SHUFFLE1_26]
%   = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from rows STARTROW
%   through ENDROW of text file FILENAME.
%
% Example:
%   [scorer,DLC_resnet50_SciaticaMar25shuffle1_180000,DLC_resnet50_SciaticaMar25shuffle1_1,DLC_resnet50_SciaticaMar25shuffle1_2,DLC_resnet50_SciaticaMar25shuffle1_3,DLC_resnet50_SciaticaMar25shuffle1_4,DLC_resnet50_SciaticaMar25shuffle1_5,DLC_resnet50_SciaticaMar25shuffle1_6,DLC_resnet50_SciaticaMar25shuffle1_7,DLC_resnet50_SciaticaMar25shuffle1_8,DLC_resnet50_SciaticaMar25shuffle1_9,DLC_resnet50_SciaticaMar25shuffle1_10,DLC_resnet50_SciaticaMar25shuffle1_11,DLC_resnet50_SciaticaMar25shuffle1_12,DLC_resnet50_SciaticaMar25shuffle1_13,DLC_resnet50_SciaticaMar25shuffle1_14,DLC_resnet50_SciaticaMar25shuffle1_15,DLC_resnet50_SciaticaMar25shuffle1_16,DLC_resnet50_SciaticaMar25shuffle1_17,DLC_resnet50_SciaticaMar25shuffle1_18,DLC_resnet50_SciaticaMar25shuffle1_19,DLC_resnet50_SciaticaMar25shuffle1_20,DLC_resnet50_SciaticaMar25shuffle1_21,DLC_resnet50_SciaticaMar25shuffle1_22,DLC_resnet50_SciaticaMar25shuffle1_23,DLC_resnet50_SciaticaMar25shuffle1_24,DLC_resnet50_SciaticaMar25shuffle1_25,DLC_resnet50_SciaticaMar25shuffle1_26] = importfile('Sciatica_01152019_R05_on_T09DLC_resnet50_SciaticaMar25shuffle1_180000filtered.csv',4, 1205);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2020/04/06 16:04:57

%% Initialize variables.
delimiter = ',';
if nargin<=2
    startRow = 4;
    endRow = inf;
end

%% Format for each line of text:
%   column1: double (%f)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
%	column6: double (%f)
%   column7: double (%f)
%	column8: double (%f)
%   column9: double (%f)
%	column10: double (%f)
%   column11: double (%f)
%	column12: double (%f)
%   column13: double (%f)
%	column14: double (%f)
%   column15: double (%f)
%	column16: double (%f)
%   column17: double (%f)
%	column18: double (%f)
%   column19: double (%f)
%	column20: double (%f)
%   column21: double (%f)
%	column22: double (%f)
%   column23: double (%f)
%	column24: double (%f)
%   column25: double (%f)
%	column26: double (%f)
%   column27: double (%f)
%	column28: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
frame = dataArray{:, 1};

nosex = dataArray{:, 2};
nosey = dataArray{:, 3};
nosel = dataArray{:, 4};

RightForex = dataArray{:, 5};
RightForey = dataArray{:, 6};
RightForel = dataArray{:, 7};
LeftForex = dataArray{:, 8};
LeftForey = dataArray{:, 9};
LeftForel = dataArray{:, 10};

RightHindToex = dataArray{:, 11};
RightHindToey = dataArray{:, 12};
RightHindToel = dataArray{:, 13};
    RightHindMidx = dataArray{:, 14};
    RightHindMidy = dataArray{:, 15};
    RightHindMidl = dataArray{:, 16};
        RightHindHeelx = dataArray{:, 17};
        RightHindHeely = dataArray{:, 18};
        RightHindHeell = dataArray{:, 19};
LeftHindToex = dataArray{:, 20};
LeftHindToey = dataArray{:, 21};
LeftHindToel = dataArray{:, 22};
    LeftHindMidx = dataArray{:, 23};
    LeftHindMidy = dataArray{:, 24};
    LeftHindMidl = dataArray{:, 25};
        LeftHindHeelx = dataArray{:, 26};
        LeftHindHeely = dataArray{:, 27};
        LeftHindHeell = dataArray{:, 28};

MirrorRightHindx =	dataArray{:, 29};
MirrorRightHindy =	dataArray{:, 30};
MirrorRightHindl =  dataArray{:, 31};
MirrorLeftHindx	= dataArray{:, 32};
MirrorLeftHindy	= dataArray{:, 33};
MirrorLeftHindl = dataArray{:, 34};

Backx = dataArray{:, 35};
Backy = dataArray{:, 36};
Backl = dataArray{:, 37};

TailBasex = dataArray{:, 38};
TailBasey = dataArray{:, 39};
TailBasel = dataArray{:, 40};

ABforcePlateCenterx = dataArray{:, 41};
ABforcePlateCentery = dataArray{:, 42};
ABforcePlateCenterl = dataArray{:, 43};
    CDforcePlateCenterx = dataArray{:, 44};
    CDforcePlateCentery = dataArray{:, 45};
    CDforcePlateCenterl = dataArray{:, 46};


end