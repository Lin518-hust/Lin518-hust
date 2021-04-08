function [data_x] = derWENOr2_plus(data, dx)


% constant parameters
ep = 1.d-6;


% output variable
data_x = zeros(length(data)-4,1);


% absorbing BC: p(x)=0, for |x|>=1;
data(2)=0; data(1)=0;
data(end-1)=0; data(end)=0;

%Generate the divided difference tables, D1(j)=Delta+f[j],
%            D2(j) = Delta-(Delta+ f[j]))
D1 = zeros(size(data)); D2=zeros(size(data));
D1(1:end-1) = (data(2:end)-data(1:end-1));
D2(2:end-1) = data(3:end)-2*data(2:end-1)+data(1:end-2);

for j=1:(length(data)-4)

    i = j+2;
    rp = (ep + D2(i+1)^2)/(ep + D2(i)^2);
    wp = 1/(1+2*rp^2);
    data_x(j) = ((D1(i-1)+D1(i)) - wp*(D1(i+1)-2*D1(i)+D1(i-1)))/2;

end

