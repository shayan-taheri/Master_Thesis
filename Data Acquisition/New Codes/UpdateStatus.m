function UpdateStatus(src,event)
% Update Temperature Measurement and Update Status

global temp;
global cond;
global s;

temp(1,:) = (temp(1,:) + mean(event.Data)) / 2;

if (cond ~= 0)
    stop(s);
end

end