function [ FTGs ] = granulations( data0,parameter )

    FTGs = cell(1,parameter.n); 
    parfor i =1:parameter.n
        data = data0(i,:);
        data = data.';   
        lambda_max = l1tf_lambdamax(data);  

        [z1,~] = l1tf(data, parameter.y*lambda_max);
        breakpoints = findbreakpoints(z1);
        %     maxnum = length(breakpoints)-1;
        FTG = struct();
        FTG.FTIGs = [];
        FTG.bps = breakpoints;
        [FTG.FTIGs,FTG.bps] = merging(data.',FTG.bps,parameter);
        FTGs{i}(end+1) = FTG;
    end
    disp('granualar finish');
end

