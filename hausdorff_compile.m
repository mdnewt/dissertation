extra  = 'Excel Output\Extras\Extras Compiled.xlsx';
HD = struct;

HD.g1 = cell(49*10,6);
for ii = 1:10
   [hd,names] = xlsread(extra,num2str(ii))         ;
    hd        = hd(:,[3,5,7])                      ;
    HD.g1(49*(ii-1)+1:49*(ii),1  ) = {1}           ;
    HD.g1(49*(ii-1)+1:49*(ii),2  ) = {ii}          ;
    HD.g1(49*(ii-1)+1:49*(ii),3  ) = names(2:end,1);
    HD.g1(49*(ii-1)+1:49*(ii),4:6) = num2cell(hd)  ;
end

HD.g3 = cell(47*10,6);
for ii = 1:10
   [hd,names] = xlsread(extra,num2str(ii+10))      ;
    hd        = hd(:,[3,5,7])                      ;
    HD.g3(47*(ii-1)+1:47*(ii),1  ) = {3}           ;
    HD.g3(47*(ii-1)+1:47*(ii),2  ) = {ii+10}       ;
    HD.g3(47*(ii-1)+1:47*(ii),3  ) = names(2:end,1);
    HD.g3(47*(ii-1)+1:47*(ii),4:6) = num2cell(hd)  ;
end

HD.g5 = cell(45*10,6);
for ii = 1:10
   [hd,names] = xlsread(extra,num2str(ii+20))      ;
    hd        = hd(:,[3,5,7])                      ;
    HD.g5(45*(ii-1)+1:45*(ii),1  ) = {5}           ;
    HD.g5(45*(ii-1)+1:45*(ii),2  ) = {ii+20}       ;
    HD.g5(45*(ii-1)+1:45*(ii),3  ) = names(2:end,1);
    HD.g5(45*(ii-1)+1:45*(ii),4:6) = num2cell(hd)  ;
end

HD.g10 = cell(40*10,6);
for ii = 1:10
   [hd,names] = xlsread(extra,num2str(ii+30))       ;
    hd        = hd(:,[3,5,7])                       ;
    HD.g10(40*(ii-1)+1:40*(ii),1  ) = {10}          ;
    HD.g10(40*(ii-1)+1:40*(ii),2  ) = {ii+30}       ;
    HD.g10(40*(ii-1)+1:40*(ii),3  ) = names(2:end,1);
    HD.g10(40*(ii-1)+1:40*(ii),4:6) = num2cell(hd)  ;
end

HD.g20 = cell(30*10,6);
for ii = 1:10
   [hd,names] = xlsread(extra,num2str(ii+40))       ;
    hd        = hd(:,[3,5,7])                       ;
    HD.g20(30*(ii-1)+1:30*(ii),1  ) = {20}          ;
    HD.g20(30*(ii-1)+1:30*(ii),2  ) = {ii+40}       ;
    HD.g20(30*(ii-1)+1:30*(ii),3  ) = names(2:end,1);
    HD.g20(30*(ii-1)+1:30*(ii),4:6) = num2cell(hd)  ;
end

% Compile and save as Excel data
HD_all = [{'G','It','Name','HD_AC','HD_Bone','HD_Fill'};...
    HD.g1;HD.g3;HD.g5;HD.g10;HD.g20];
xlswrite(extra,HD_all,'HD List');