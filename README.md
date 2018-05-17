# multi-scale-event-detection
Dataset: 
The data size is 32736 and there are 15 features for each data. We are going to take four features, “name”, “number”, “rating”, “add_time” to implement multi-scale event detection. 
<h2>How to use:</h2> 1. Run data_seg.py <br>2. Run run.mat First, We run the data_seg.py to separate the Chinese words with Jieba, so we generate the segmentation words dictionary (“data.xlsx”) . Later, we just need to run the “run.mat” whenever the input data file is named as “origin.xlsx”. The result is the clusters that represent different events. 
 
