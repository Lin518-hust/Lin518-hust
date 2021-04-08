function [data_x] = derWENOr2_minus(data, dx)



ep = 1.d-6;


% output variable
data_x = zeros(length(data)-4,1);


% absorbing BC
data(2)=0; data(1)=0;
data(end-1)=0; data(end)=0;


D1 = zeros(size(data)); D2=zeros(size(data));
D1(1:end-1) = (data(2:end)-data(1:end-1));
D2(2:end-1) = data(3:end)-2*data(2:end-1)+data(1:end-2);

for j=1:(length(data)-4)

    i = j+2;
    rm = (ep + D2(i-1)^2)/(ep + D2(i)^2);
    wm = 1/(1+2*rm^2);
    data_x(j) = ((D1(i-1)+D1(i)) - wm*(D1(i-2)-2*D1(i-1)+D1(i)))/2;

end

