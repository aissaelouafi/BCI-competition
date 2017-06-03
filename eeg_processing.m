% LOAD EEG signals data 
% fs = 100 Hz

mat = dir('./data/*.mat');
for K = 1 : length(mat)
  thisfile = mat(K).name;
  %disp(strcat('./data/',thisfile));
  thisdata = load(strcat('./data/',thisfile));
end
