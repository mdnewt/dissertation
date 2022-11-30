 % Perform initial setup
 xlsx      =     [pwd,'..\..\..\Scan Log.xlsx']            ;
 p         =     '1 - Images\'                             ;
 fCORCHK   = @(x) isempty(whos('-file',x,'correction_log'));
[~,cor]    =      xlsread(xlsx,'Correction Log','A:C')     ;
 nirfields =    {'x700','x800','W'}                        ;
 %%
 for ii = 3:size(cor,1)
     %%
     name = cell(1,2); fname = name; file = name;
     
     % Which two samples are flipped?
     name{1} = cor{ii,1}; fname{1} = [p,name{1},'.mat'];
     name{2} = cor{ii,2}; fname{2} = [p,name{2},'.mat'];
     
     % Define correction in log; ensure it has not already been applied
     correction_log = {strjoin(cor(ii,:),'_')}; JJ = [];
     for jj = 1:2
         if ~fCORCHK(fname{jj})
              c  =  load(fname{jj},'correction_log');
              if ~any(ismember(c.correction_log,correction_log))
                  JJ = [JJ,jj];
              end
         else JJ = [JJ,jj];
         end
     end
     JJ = [1,2]; % Overrule

     % Do both samples have image files? Does either still need to be completed?
     if exist(fname{1},'file') && exist(fname{2},'file') && ~isempty(JJ)
         %%
         
         % Load files
         tp = lower(cor{ii,3});
         switch tp
             case 'pre'
                 file{1} = load(fname{1},'nir'      );
                 file{2} = load(fname{2},'nir'      );
             case 'post'
                 file{1} = load(fname{1},'nir','vct');
                 file{2} = load(fname{2},'nir','vct');
             case 'ceuct'
                 file{1} = load(fname{1},'ceuct'    );
                 file{2} = load(fname{2},'ceuct'    );
         end         

         %% Overall loop for applying corrections to these two files
         for jj = JJ
             %%
             % Grab set of images to be corrected
             switch tp
                 case 'pre'
                     nir   =   file{jj}.nir  ; fields =  fieldnames(nir);
                 case 'post'
                     nir   =   file{jj}.nir  ; fields =  fieldnames(nir);
                     vct   =   file{jj}.vct  ;
                 case 'ceuct'
                     ceuct =   file{jj}.ceuct; fields = {'stack'}       ;
             end
             
             % Loop for each scan orientation (condyles v. trochlea)
             for kk = 1:length(fields)
                 %%
                 % Apply corrections based on correction type
                 switch tp
                     case 'pre'
                         % Correct NIR
                         nirfields = fieldnames(nir.(fields{kk}).(tp));
                         for mm = 1:length(nirfields)
                             nir.(fields{kk}).(tp).(nirfields{mm}) = ...
                                 file{rem(jj,2)+1}.nir.(fields{kk}).(tp).(nirfields{mm});
                         end
                     case 'post'
                         % Correct NIR
                         nirfields = fieldnames(nir.(fields{kk}).(tp));
                         for mm = 1:length(nirfields)
                             nir.(fields{kk}).(tp).(nirfields{mm}) = ...
                                 file{rem(jj,2)+1}.nir.(fields{kk}).(tp).(nirfields{mm});
                         end
                         vctfields = fieldnames(vct.(fields{kk})     );
                         for mm = 1:length(vctfields)
                             vct.(fields{kk}).(vctfields{mm}) = ...
                                 file{rem(jj,2)+1}.vct.(fields{kk}).(vctfields{mm});
                         end
                     case 'ceuct'
                         ceuct.(fields{kk}) = file{rem(jj,2)+1}.ceuct.(fields{kk});
                 end
             end
             
             % Save
             switch tp
                 case 'pre'
                     if ~fCORCHK(fname{jj})
                          c  = load(fname{jj},'correction_log');
                          c2 = correction_log;
                          correction_log = [c.correction_log;correction_log];
                          save(fname{jj},'nir'      ,'correction_log','-append');
                          correction_log = c2; clear c c2
                     else save(fname{jj},'nir'      ,'correction_log','-append');
                     end
                 case 'post'
                     if ~fCORCHK(fname{jj})
                          c  = load(fname{jj},'correction_log');
                          c2 = correction_log;
                          correction_log = [c.correction_log;correction_log];
                          save(fname{jj},'nir','vct','correction_log','-append');
                          correction_log = c2; clear c c2
                     else save(fname{jj},'nir','vct','correction_log','-append');
                     end
                 case 'ceuct'
                     if ~fCORCHK(fname{jj})
                          c  = load(fname{jj},'correction_log');
                          c2 = correction_log;
                          correction_log = [c.correction_log;correction_log];
                          save(fname{jj},'ceuct','correction_log','-append');
                          correction_log = c2; clear c c2
                     else save(fname{jj},'ceuct','correction_log','-append');
                     end
             end
         end
         
     fprintf(['Corrections applied: ',correction_log{end},'\n']);
     
     end
 end