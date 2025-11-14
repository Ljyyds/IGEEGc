function [ FTIG,bps ] = merging( data,bps,parameter)

n = length(bps);

minNum = parameter.ftigs;
while (n-1)>minNum
    tempLI = zeros(1,n-2);
    for i = 1:n-2
        [k,~] = lineRe(data,bps(i),bps(i+2));
        sigma = std(data(bps(i):bps(i+2)));
        len  = bps(i+2)-bps(i);
        tempLI(i) = sigma*sqrt(len/(1+k.^2));
    end
    [~,j] = min(tempLI);
    bps(j+1) = [];
    n = n-1;
end

i = 2;
while i <= n-1
    if bps(i)-bps(i-1)<=12
        bps(i) = [];
        n = n-1;
        i = i-1;
    end
    i = i+1;
    if i ==n 
        break;
    end
end


if bps(n)-bps(n-1) <=12
    bps(n-1) = [];
    n = n-1;
end

FTIG = struct('aopt',{},'bopt',{},'ki',{},'bi',{},'W',{});
U = cell(1,n-1);

Fs = parameter.fs;

for i=1:n-1
    [FTIG(i).ki,FTIG(i).bi] = lineRe(data,bps(i),bps(i+1)); 
    Leni = bps(i+1)-bps(i);
    U{i} = zeros(1,Leni);
    for t = 1:Leni
        U{i}(t) =data(bps(i)+t-1)- t*FTIG(i).ki-FTIG(i).bi;
    end
    [aopt,bopt] = au_bu(U,i);
    FTIG(i).aopt = aopt+FTIG(i).bi;
    FTIG(i).bopt = bopt+FTIG(i).bi;
    
    eeg_signal = data(bps(i):bps(i+1));
    if Fs == 250 || Fs==400 || Fs==256
        upw = 100;
    elseif Fs == 100 || Fs == 173.61
        upw = 47.8;
    end

    FTIG(i).W = bandpower(eeg_signal, Fs, [0.5,upw]);
end
end

