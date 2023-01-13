
for i=1:length(output)
    try
    table(i).fileName = output(i).filname;
    table(i).note = output(i).byEye;
    table(i).numberOfCycles = output(i).numcycles;
    table(i).direction = output(i).direction;
     output(i).velocityNose =  ( abs(gradient(output(i).N.x.avg)) ) ./ (gradient(output(i).N.f.avg));%pixels per frame
    table(i).avgVnose = mean(output(i).velocityNose);
     output(i).velocityBack =  ( abs(gradient(output(i).B.x.avg)) ) ./ (gradient(output(i).B.f.avg));%pixels per frame
    table(i).avgVback = mean(output(i).velocityBack);
    table(i).AverageVelocityN = output(i).AverageVelocity;
    table(i).maxHTy = max(output(i).HT.y.avg);
    table(i).df = output(i).Hdutyfactor.avg;
    table(i).stridelength = output(i).Hstridelength.avg;
    table(i).toeProm = output(i).ToeProminance.avg;
    table(i).heelProm = output(i).HeelProminance.avg;
    table(i).FootAngleToeOff = output(i).FootAngleToeOff.avg; 
    table(i).stanceDuration = output(i).stanceDuration.avg; 
    if output(i).stepwidth.variability<.5*output(i).stepwidth.avg
        table(i).stepwidth = output(i).stepwidth.avg;
    else
        table(i).stepwidth = NaN;
    end
    end
end
%% velocity
%VARx table
for i=1:length(output)
    try
    mylength(i)=length(output(i).velocity);
    end
end
mylength=max(mylength);
varx=nan(length(output),mylength);
for i=1:length(output)
    try
    myvel=output(i).velocity;
    varx(i,1:length(myvel))=myvel;
    end
end
%filename = 'gait3_1.xlsx';
%writematrix(varx,filename,'Sheet',2)    
%VARf table
for i=1:length(output)
    try
    mylength(i)=length(output(i).N.f.avg);
    end
end
mylength=max(mylength);
varf=nan(length(output),mylength);
for i=1:length(output)
    try
    myvel=output(i).N.f.avg;
    varf(i,1:length(myvel))=myvel;
    end
end
%writematrix(varf,filename,'Sheet',3)  
%% var of interest

%VARx table
for i=1:length(output)
    try
    mylength(i)=length(output(i).B.x.avg);
    end
end
mylength=max(mylength);
varx=nan(length(output),mylength);
for i=1:length(output)
    try
    myvel=output(i).B.x.avg;
    varx(i,1:length(myvel))=myvel;%/output(i).pixpercm;
    end
end
%writematrix(varx,filename,'Sheet',4)  

%VARy table
for i=1:length(output)
    try
    mylength(i)=length(output(i).HT_TB.y.avg);
    end
end
mylength=max(mylength);
vary=nan(length(output),mylength);
for i=1:length(output)
    try
    myvel=output(i).HT_TB.y.avg;
    vary(i,1:length(myvel))=myvel;%/output(i).pixpercm;
    end
end
%writematrix(vary,filename,'Sheet',5)

%VARf table
for i=1:length(output)
    try
    mylength(i)=length(output(i).HT_TB.f.avg);
    end
end
mylength=max(mylength);
varf=nan(length(output),mylength);
for i=1:length(output)
    try
    myvel=output(i).HT_TB.f.avg;
    varf(i,1:length(myvel))=myvel;
    end
end
%writematrix(varf,filename,'Sheet',6)  