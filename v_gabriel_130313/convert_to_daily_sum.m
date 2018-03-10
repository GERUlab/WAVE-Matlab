function [timeday daily_var] = convert_to_daily_sum(var,time)

%Convert data to daily data by summation of all data of the same day
%Missing values are not noted as NaN, but are skipped in this conversion,
%for the reason that NaN can not be used as a boundary condition is wavemat
daynumber = floor (time);
minday = min(daynumber);
maxday = max(daynumber);
i=minday;
linenumber =1;
while i <= maxday
index = find(daynumber==i);
if isempty(index)==0
    if size(index,1)>=23;
    s = cumsum(var(index));
    daily_var(linenumber) = s(end);
    timeday(linenumber)=i;
    linenumber = linenumber +1;
    end
end
i=i+1;
end

