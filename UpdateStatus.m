function UpdateStatus(src,event)
% Update Temperature Measurement and Update Status

global temp;

if isempty(temp)
    temp = mean(event.Data);
else
    temp = (temp + mean(event.Data)) / 2;
end

end