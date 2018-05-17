
filename ='origin.xlsx';
filename1 = 'data.xlsx';
load('dict.mat','dictionary');
dictionary = cellstr(dictionary);
[num,txt,raw] = xlsread(filename);
[sig,tstr,raw1] = xlsread(filename,1,'I2:I32727');
[num1,txt1,raw2] = xlsread(filename1);
name = txt(:,3);
value = num(:,13:14);
c= 1;

name(1,:) = [];
tdmatrix = sparse(length(dictionary),32726);

for i = 1:32726
    nam = txt1(i,:);
    x = find(~cellfun(@isempty,nam));
    for j = 1:size(x)
        
         na = nam(x(j));
         naa =na{:};
         for k = 1: size(dictionary)
             dic = dictionary(k);
             di = dic{:};
             if(size(di,2)==size(naa,2))
                if(di==naa)
                    tdmatrix(k,i) = 1;
                end
             end
         end
        disp(i);
    end
end


for i = 1:32726
    tr= strsplit(tstr{i},' ');

    time(i) = tr(1,3);
    tr1= strsplit(tr{3},':');
    time1(i,1) = tr1(1,1);
    time1(i,2) = tr1(1,2);

%tim(i)= datenum(time(i),'HH:MM:SS');
end
timeh = str2num(cell2mat(time1(:,1)));
timem = str2num(cell2mat(time1(:,2)));
timem = timem + timeh*60;


%timem = gpuArray(timem);
tdm_tfidf = tmg_tdm(tdmatrix);
%tdm_tfidf = gpuArray(tdm_tfidf);
% setup the 2-d space

xmin = min(value(:,1));
xmax = max(value(:,1));
ymin = min(value(:,2));
ymax = max(value(:,2));
tmin = min(timem);
tmax = max(timem);


N = size(tdm_tfidf,2);


corr_nan = cell(1,11);
wav_init_res_time = 60;
wav_init_res_space_x = 100;
wav_init_res_space_y = 0.5;
%%%%%%%%%%%% Wavelet-based clustering %%%%%%%%%%%%
doc_per_word = sum(tdmatrix>0,2);
frequentwordidx = find(doc_per_word>=10); % idx of frequent words in dictionary
frequentdocidx = find(sum(tdmatrix(frequentwordidx,:),1))'; % idx of docs containing frequent words

fdict = dictionary(frequentwordidx);
ftdm = tdmatrix(frequentwordidx,frequentdocidx)>0;
df = sum(ftdm,2);
ftweetgeo = value(frequentdocidx,:);
ftweettime = timem(frequentdocidx,:);

% set initial time interval
ftweettimeslot = ceil(ftweettime/wav_init_res_time);

% set initial spatial interval
lat = ymax:-wav_init_res_space_y:ymin;
if lat(end)~=ymin
    lat = [lat ymin];
end
lon = xmin:wav_init_res_space_x:xmax;
if lon(end)~=xmax
    lon = [lon xmax];
end
latgrid = zeros(length(lat)-1,1);
longrid = zeros(length(lon)-1,1);

% associate each grid a gps coordinates which is the center of the grid
for i = 1:length(latgrid)
    latgrid(i) = (lat(i)+lat(i+1))/2;
end
for i = 1:length(longrid)
    longrid(i) = (lon(i)+lon(i+1))/2;
end
tweetindex1 = cell(length(latgrid),length(longrid),max(ftweettimeslot));
tweetindex2 = zeros(length(frequentdocidx),3);

for i = 1:length(latgrid)
    for j = 1:length(longrid)
        idx = find( ftweetgeo(:,2)<=lat(i) & ftweetgeo(:,2)>=lat(i+1) & ftweetgeo(:,1)>=lon(j) & ftweetgeo(:,1)<=lon(j+1) );
        for k = 1:max(ftweettimeslot)
            tmpidx = idx(ftweettimeslot(idx)==k);
            tweetindex1{i,j,k} = tmpidx;
            if ~isempty(tmpidx)
                tweetindex2(tmpidx,:) = repmat([i j k],length(tmpidx),1);
            end
            
        end
    end
    disp(i);
end

% compute spatial scales
maxdist = norm([latgrid(1) longrid(1)]-[latgrid(end) longrid(end)],2)+0.01;
mindist = norm([latgrid(1) longrid(1)]-[latgrid(2) longrid(1)],2)-0.01;
nscales = 4;
distscales = exp(linspace(log(mindist),log(maxdist),nscales+1));
sim = sparse(length(frequentdocidx),length(frequentdocidx));
[eidx1,eidx2] = find(triu(double(ftdm)'*double(ftdm),1)); % get tweet pairs that share common keywords

% compute time series for each keyword
for kkk = 1:size(ftdm,1)
    ts{kkk} = accumarray(gather(tweetindex2),full(ftdm(kkk,:)'),[],@(x)sum(x));
    ts{kkk} = cat(3,ts{kkk},zeros(size(ts{kkk},1),size(ts{kkk},2),ceil(tmax/wav_init_res_time))); % zero padding to ensure each time series has length 128
    cc{kkk} = ts{kkk};
    disp(kkk);
end

g=0;
for m = 1:length(eidx1)
    
    tmppos1 = tweetindex2(eidx1(m),1:2);
    tmppos2 = tweetindex2(eidx2(m),1:2);
    
    if sum(tmppos1==tmppos2) == 2
        tmpsim = 1;
        % two tweets come from the same cell, two time series are the same
        % spatiotemporal similarity of any common word is 1
        % equivalent to only keeping the spatial constraint but relaxing
        % the temporal constraint in LED
    else
        tmpdist = norm([latgrid(tmppos1(1)) longrid(tmppos1(2))] - [latgrid(tmppos2(1)) longrid(tmppos2(2))],2);
        % set the temporal scale according to the spatial scale
        tmpscale = find((tmpdist-distscales)<=0,1)-1; % 1,2,3,4
        tmpscale = nscales+1-tmpscale; % 4,3,2,1
        %tmpscale = min(tmpscale,maxwaveletscale);
        tmpscale = gather(tmpscale);
        % for each of the common words, compute correlations of wavelet packet
        % coefficients at the determined temporal scale
        commonwords = intersect(find(ftdm(:,eidx1(m))),find(ftdm(:,eidx2(m))));
        for n = 1:length(commonwords)

            tmpts1 = ts{commonwords(n)}(tmppos1(1),tmppos1(2),:);
            tmpts1 = reshape(tmpts1,[],size(tmpts1,2),1);
            tmpts2 = ts{commonwords(n)}(tmppos2(1),tmppos2(2),:);
            tmpts2 = reshape(tmpts2,[],size(tmpts2,2),1);
            Wcoeff1 = wavedec(tmpts1,tmpscale,'haar'); Wcoeff1 = Wcoeff1(1:num_dwt_coeff(length(tmpts1),tmpscale));
            Wcoeff2 = wavedec(tmpts2,tmpscale,'haar'); Wcoeff2 = Wcoeff2(1:num_dwt_coeff(length(tmpts2),tmpscale));
            
            if isnan(corr(Wcoeff1,Wcoeff2)) == 1
                corr_nan{1} = [corr_nan{1};m];
            end
            tmpsim(n) = max(corr(Wcoeff1,Wcoeff2),0);
        end
    end
    % the overall similarity will be maximum correlation times textual similarity
    sim(eidx1(m),eidx2(m)) = max(tmpsim)*real(dot(tdm_tfidf(:,frequentdocidx(eidx1(m))),tdm_tfidf(:,frequentdocidx(eidx2(m)))));
    sim(eidx2(m),eidx1(m)) = sim(eidx1(m),eidx2(m));
    disp(g);
    g = g+1;
end


% map the similarity matrix of frequent docs to similarity matrix of original docs
multiscalegraph = sparse(N,N);
[a,b] = find(sim);
for i = 1:length(a)
    multiscalegraph(frequentdocidx(a(i)),frequentdocidx(b(i))) = sim(a(i),b(i));
    disp(i);
end

[COMTY2, ending2] = cluster_jl(multiscalegraph,1,1,1,1);
cluster_idx2 = COMTY2.COM{1,1}';
clustercount2 = accumarray(cluster_idx2,1);
