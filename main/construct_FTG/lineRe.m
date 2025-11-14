function [ k,c ] = lineRe( data,left,right )
time_points = 1:numel(data(left:right));

mean_time_series = mean(data(left:right));
mean_time_points = (right - left)/2;

k = sum((time_points - mean_time_points) .* (data(left:right) - mean_time_series)) / sum((time_points - mean_time_points) .^ 2);
c = mean_time_series - k * mean_time_points;
end

