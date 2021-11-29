%% now let's normalize - have to add pixpercm 1st tho
for i=1:length(output)
    try
    output(i).velocity =  ( abs(gradient(output(i).N.x.avg)) ./ output(i).pixpercm ) ./ (gradient(output(i).N.f.avg) ./500);
    output(i).velocityB =  ( abs(gradient(output(i).B.x.avg)) ./ output(i).pixpercm ) ./ (gradient(output(i).B.f.avg) ./500);
    output(i).velocitystd = std(output(i).velocity);
    output(i).avgV = mean(output(i).velocity);
    
    output(i).Hstridelength=output(i).Hstridelength.avg/output(i).pixpercm;
    end
end

%VARx table
for i=1:length(output)
    try
    mylength(i)=length(output(i).HT.x.avg);
    end
end
mylength=max(mylength);
varx=nan(length(output),mylength);
for i=1:length(output)
    try
    myvel=output(i).HT.x.avg;
    varx(i,1:length(myvel))=myvel/output(i).pixpercm;
    end
end
    
%VARy table
for i=1:length(output)
    try
    mylength(i)=length(output(i).HT.y.avg);
    end
end
mylength=max(mylength);
vary=nan(length(output),mylength);
for i=1:length(output)
    try
    myvel=output(i).HT.y.avg;
    vary(i,1:length(myvel))=myvel/output(i).pixpercm;
    end
end

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