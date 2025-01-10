function obj = ddr_obj(x, datas, ybmls, N, lenreading)
sz = N*lenreading;

xx = reshape(x(N+1:N+sz),N,lenreading);

dev_data = datas - xx;
dev_model =  x(1:N).\ybmls - xx;

dev_data2 = dev_data.*dev_data;
dev_model2 = dev_model.*dev_model;

clear dev_data dev_model

sse_data = sum(dev_data2, 2);
sse_model = sum(dev_model2, 2);

obj = sum(sse_model) + sum(sse_data);

return

