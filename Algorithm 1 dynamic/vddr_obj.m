function obj = vddr_obj(x, ybml, N, w)

dev_model =  w.\ybml - x(1:N);

dev_model2 = dev_model.*dev_model;

sse_model = sum(dev_model2, 2);

obj = sum(sse_model) ;

return

