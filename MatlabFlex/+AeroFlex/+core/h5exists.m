function flag = h5exists(h5file, ds)
info = h5info(h5file);
allp = recList(info,'');
flag= ismember(ds, allp);
end
function lst = recList(info, prefix)
lst= {};
for i=1:length(info.Datasets)
   nm= [prefix '/' info.Datasets(i).Name];
   lst{end+1}= nm; 
end
for ig=1:length(info.Groups)
   sp= [prefix '/' info.Groups(ig).Name];
   subn= recList(info.Groups(ig), sp);
   lst= [lst, subn]; 
end
end