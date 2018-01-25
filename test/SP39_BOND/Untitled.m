Data = importdata('lable.txt');
index = [];
for i=1:1:size(Data,1)
    name = sprintf('%s',Data{i});
    filename = importdata(name);
    for j=2:1:size(filename,1)
        index = [index;i filename(j)];
    end
end


% D2 = importdata('XEN1_XENLA');
