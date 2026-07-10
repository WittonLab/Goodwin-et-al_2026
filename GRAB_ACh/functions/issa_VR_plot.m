function [vars, v, y, t, licks, reward] = issa_VR_plot(exp_name,exp_date,update_flag,ni,get_daq,plt);
% function [vars,data,data1,vs,dd] = issa_VR_plot(exp_name,exp_date,update_flag,ni,get_daq,plt)

% ISSA_VR_PLOT: It's a VR plot.
%
% [vars,data,data1] = ISSA_VR_PLOT(exp_name,exp_date,update_flag,ni)
%
% exp_name:     string that is the name of that experiment
% exp_date:     date of experiment as a string, in YYYYMMDD format
% update_flag:  0 (regular mode) or 1 (updates plot as data is acquired)
% ni:           number of experiment, if multiple runs for that same day
%
% vars:         structure with all the experiment data from Virmen
% data:         frame-by-frame data (see vars.si for column labels)
% data1:        DAQ data, typically at 1 kHz (lick and frame times)
%
% example usage: [vars] = ISSA_VR_PLOT('Mouse23','20190401',0);
%
% John Issa, Northwestern University, 2019

% 10/3/19 - should try to incorporate plot_traces.m as the plotting
% function

%% plot data
% possible file locations[
path_names = {[userpath filesep] ,'C:\Users\dgood\OneDrive - University of Exeter\Desktop\VR Data\Cohort 6\Behaviour'};
% path_names = {[userpath filesep] ,'C:\Users\dg453\OneDrive - University of Exeter\Desktop\Laura_behaviour\Behaviour'};
% path_names = {[userpath filesep] ,'C:\Users\dg453\OneDrive - University of Bristol\Desktop\230927\DG9\Behaviour'};

if ~exist('plt','var'), plt = true; end

%% get full path names
if ~exist('exp_date','var') || isempty(exp_date)
    exp_date = datestr(now,'yyyymmdd');
end
if ~exist('update_flag','var')
    update_flag = false;
end
if ~exist('ni','var') || isempty(ni)
    nis = '';
else nis = ['_' num2str(ni)]; end
if ~exist('get_daq','var') || isempty(get_daq)
    get_daq = 1;
end

file_name = [exp_name '_' exp_date nis '_virmenLog'];

n_path_names = length(path_names);
for i = 1:n_path_names
    dir_root_i = [path_names{i} '\virmenLogs' filesep];
    dir_info = dir([dir_root_i file_name '*.dat']);
    if ~isempty(dir_info)
        dnames = dir_info;
        dir_root = dir_root_i;
    end
    
    % also look in subfolder that uses experiment name
    dir_root_i = [path_names{i} '\virmenLogs' filesep exp_name filesep];
    dir_info = dir([dir_root_i file_name '*.dat']);
    if ~isempty(dir_info)
        dnames = dir_info;
        dir_root = dir_root_i;
    end
end


%% load files
% for now just use first file
if ~exist('dnames','var'), [vars,data,data1,vs,dd] = deal([]); return; end
name_root = dnames(1).name;
name_root(end-3:end) = [];
name_root = [dir_root name_root];

% load 'vars' (mat file)
if exist([name_root '3.mat'],'file')
    load([name_root '3.mat']); % load vars
    if iscell(vars)
        vars_temp = vars; vars = []; vars.worlds = vars_temp;
        vars.save_size = 5;
    end
    save_size = vars.save_size;
    
    if isfield(vars,'t') || (isfield(vars,'exper') & isfield(vars.exper,'t0'))
        tflag = 0; % all times in seconds (newer approach)
    else tflag = 1; % all times as date number format
    end
    if ~isfield(vars,'si')
        vars.si = [];
        vars.si.t = 1; % timestamp
        vars.si.p = 2:3; % x,y position
        vars.si.v = 4:5; % x,y velocity
        if save_size > 5
            vars.si.wi = 6; % current world index'
        end
        if save_size > 6
            vars.si.licks = 7; % number of licks
        end
    end
else
    save_size = 7; tflag = 0;
end

% load data (dat file)
fid = fopen([name_root '.dat'],'r');
data = fread(fid,[save_size inf],'double');
if ~update_flag, fclose(fid); % leave open for updating
end

%% Name variables
t = data(vars.si.t,:); t = t';
v = data(vars.si.v(2),:).*data(vars.si.vs,:)*100; v = v'; % velocity in cm/s
y = data(vars.si.p(2),:).*data(vars.si.vs,:)*100; y = y'; % y-position in cm
licks = data(vars.si.licks,:)>.5; licks = licks';

reward = zeros(length(t),1);
%using events1.reward_times to approximate reward timepoint (min distance)
for i = 1:length(vars.times.reward_times(:,2))
    [mi,id] = min(abs(t - vars.times.reward_times(i,2)),[],1,'linear');
    reward(id,1) = 1;
end
clearvars mi id

plot(licks)



% if exist, load daq digital data
if isfield(vars,'daq') && vars.daq.lick_times_flag && get_daq
    NCh = vars.daq.lick_times_flag;
    
    fid1 = fopen([name_root '1.dat'],'r');
    if fid1 ~= -1
        data1 = fread(fid1,inf,'ubit1');
        data1 = int8(data1);
        L = size(data1,1);
        data1 = data1(1:L-mod(L,NCh)); L = size(data1,1);
        data1 = reshape(data1,NCh,L/NCh)';
        if ~update_flag, fclose(fid1); % leave open for updating
        end
    else data1 = [];
    end
else fid1 = -1; data1 = [];
end




% % get imaging frame times
% i = 3;
% 
% ii = find(t2>vars.times.start_times(i,2) & t2<vars.times.end_times(i,2));
% tic
% [pks,locs] = findpeaks(data1(ii,2));
% toc
% 
% disp(size(locs))
% 
% plot(t2,data1(:,2),'k')
% hold on
% plot(t2(min(ii)+(locs-1)),data1(min(ii)+(locs-1),2),'g.')
% hold off


% load lap performance (txt file)
try
    dd=dlmread([name_root '2.txt']);
    dd_ind = 6:3:size(dd,2); % find all timestamps
    if tflag
        dd(:,3) = datevec(dd(:,3)-data(1,1))*[0,0,24*60*60,60*60,60,1]';
        dd(:,dd_ind) = datevec(dd(:,dd_ind)-data(1,1))*[0,0,24*60*60,60*60,60,1]';
    end
catch ME
    dd = [];
    dd_ind = [];
end




%% plot data
if tflag
t = datevec((data(vars.si.t,:)-data(vars.si.t,1)))*[0,0,24*60*60,60*60,60,1]';
else t = data(vars.si.t,:);
end
Nw = length(vars.worlds);
if isfield(vars.si,'wi')
    vs_w = NaN(Nw,1);
    for wi = 1:Nw
        vs_w(wi) = vars.worlds{wi}.virmenScale;
    end
%     vs = vs_w(data(vars.si.wi,:)); vs = vs(:)';
% else
    vs = ones(size(t));
end

% return here if not plotting
if ~plt
    return
end


if isfield(vars.si,'vs')
%                 v = data(vars.si.v(2),:)./data(vars.si.vs,:); % scale velocity vector
    v = data(vars.si.v(2),:).*data(vars.si.vs,:)*100;% ./data(vars.si.dt,:); % scale velocity vector to cm/s
    p = data(vars.si.p(2),:).*data(vars.si.vs,:)*100;
else
    if isfield(vars.worlds{1},'virmenScale')
        v = data(vars.si.v(2),:).*vs*100;
        p = data(vars.si.p(2),:).*vs*100;
    else
        v = data(vars.si.v(2),:);
        p = data(vars.si.p(2),:);
    end
end
% v = smooth(v,3);
ax = [];
lh = [];

% ax(1)=subplot(2,1,1);
% % plot(repmat(15:15:max(30,max(t)),2,1),[0 ceil(max(data(vars.si.v(2),:)))],'color',[.7 .7 .7])
% % plot(repmat(15:15:max(30,max(t)),2,1),[0 ceil(max(v))],'color',[.7 .7 .7])
% plot(repmat(15:15:max(30,max(t)),2,1),[0 50],'color',[.7 .7 .7])
% hold on
% plot(t([1 end]),[0 0],'color',[.5 .5 .5])
% % lh(1) = plot(t,data(vars.si.v(2),:),'k');
% lh(1) = plot(t,v,'k');
% 
% if isfield(vars,'times')
%     if isfield(vars.times,'reward_times')
%         if ~isempty(vars.times.reward_times)
%             plot(repmat(vars.times.reward_times(:,2),1,2),[0 50],'-','linewidth',.75,'color',[0 .35 .8])
%         end
%     end
% end
% hold off
% axis off
% % ylim([-1 5])
% ylim([-10 90])
% 
% ax(2)=subplot(2,1,2);
% if isfield(vars.si,'p')
% % plot(repmat(15:15:max(30,max(t)),2,1),[0 ceil(max(data(vars.si.p(2),:)))],'color',[.8 .8 .8])
% plot(repmat(15:15:max(30,max(t)),2,1),[0 ceil(max(p))],'color',[.8 .8 .8])
% hold on
% % lh(2) = plot(t,data(vars.si.p(2),:),'k');
% lh(2) = plot(t,p,'k');
% end
% 
% % if ~isempty(dd), lh(4) = plot(dd(:,3),dd(:,4)*max(data(vars.si.p(2),:)),'go'); end
% % if ~isempty(dd), lh(4) = plot(dd(:,3),dd(:,4)*max(p),'go'); end
% if isfield(vars,'times')
%     if isfield(vars.times,'reward_times')
%         if ~isempty(vars.times.reward_times)
%             %         plot(vars.times.reward_times(:,2),10,'go')
%             plot(repmat(vars.times.reward_times(:,2),1,2),[15 50],'-','linewidth',.75,'color',[0 .35 .8])
%         end
%     end
% end
% % for i_d = 1:length(dd_ind), plot(dd(:,dd_ind(i_d)),dd(:,dd_ind(i_d)-1),'rs'); end
% for i_d = 1:length(dd_ind), plot(dd(:,dd_ind(i_d))*[1 1],[15 50],'r-'); end
% % if ~isempty(vars.worlds{1}.loom_times)
% %     plot(datevec((vars.worlds{1}.loom_times-data(1,1)))*[0,0,24*60*60,60*60,60,1]',0,'rx') % loom times
% % end
% if isfield(vars.si,'licks')
% %     lh(3) = plot(t,data(vars.si.licks,:),'color',[.4 .65 1]);
%     lh(3) = plot(t,(data(vars.si.licks,:)>2)*10,'color',[.2 .5 1]);
% end
% hold off
% axis off
% ylim([0 420])
% set(gcf,'color','w')
% linkaxes(ax,'x')
% % xlim(t(end)+[-300 10])
% fprintf('%d minutes, %d laps\n',[round(t(end)/60) size(dd,1)])

%% update
if update_flag
    cc = ''; % current character
    xrange = 200; % seconds to look back
    set(gcf,'CurrentCharacter','a')
    while ~strcmp(cc,'q')
        data = [data fread(fid,[save_size inf],'double')];
        
        if fid1 ~= -1
            data1new = fread(fid1,inf,'ubit1');
            L = size(data1new,1);
            data1new = data1new(1:L-mod(L,NCh)); L = size(data1new,1); % this may not work properly
            data1new = reshape(data1new,NCh,L/NCh)';
            data1 = [data1; data1new];
        end

        if tflag
            t = datevec((data(vars.si.t,:)-data(vars.si.t,1)))*[0,0,24*60*60,60*60,60,1]';
        else t = data(vars.si.t,:);
        end
        if isfield(vars.si,'vs')
%             v = data(vars.si.v(2),:)./data(vars.si.vs,:); % scale velocity vector
    v = data(vars.si.v(2),:).*data(vars.si.vs,:)*100;% ./data(vars.si.dt,:); % scale velocity vector to cm/s
    p = data(vars.si.p(2),:).*data(vars.si.vs,:)*100;
        else 
            v = data(vars.si.v(2),:);
            p = data(vars.si.p(2),:);
        end
        
        try
            dd2=dlmread([name_root '2.txt'],'',size(dd,1)-1,0); % read just new data
            dd_ind = 6:3:size(dd2,2); % find all timestamps
            if tflag
                dd2(:,3) = datevec(dd2(:,3)-data(1,1))*[0,0,24*60*60,60*60,60,1]';
                dd2(:,dd_ind) = datevec(dd2(:,dd_ind)-data(vars.si.t,1))*[0,0,24*60*60,60*60,60,1]';
            end
            
            %             dd2(:,3) = datevec(dd2(:,3)-data(1,1))*[0,0,24*60*60,60*60,60,1]';
            %             dd_ind = 6:3:size(dd2,2); % find all timestamps
            %             dd2(:,dd_ind) = datevec(dd2(:,dd_ind)-data(1,1))*[0,0,24*60*60,60*60,60,1]';
        catch ME
            dd2 = [];
            dd_ind2 = [];
        end
        dd = [dd; dd2];
        dd_ind = 6:3:size(dd,2); % find all timestamps
        if ~isempty(dd2), fprintf('%d minutes, %d laps\n',[round(t(end)/60) size(dd,1)]); end
        
%         set(lh(1),'xdata',t,'ydata',data(vars.si.v(2),:))
        set(lh(1),'xdata',t,'ydata',v)
%         if isfield(vars.si,'p'), set(lh(2),'xdata',t,'ydata',data(vars.si.p(2),:)); end
        if isfield(vars.si,'p'), set(lh(2),'xdata',t,'ydata',p); end
%         if isfield(vars.si,'licks'), set(lh(3),'xdata',t,'ydata',data(vars.si.licks,:)); end
        if isfield(vars.si,'licks'), set(lh(3),'xdata',t,'ydata',(data(vars.si.licks,:)>2)*10); end
        %         if size(data,1)>=7, set(lh(3),'xdata',t,'ydata',data(7,:)); end
%         if ~isempty(dd), set(lh(4),'xdata',dd(:,3),'ydata',dd(:,4)*max(data(vars.si.p(2),:))); end
%         xlim(t(end)+[-200 10])
        xlim(t(end)+[-1 .05]*xrange)
        pause(.15)
        cc = get(gcf,'CurrentCharacter');
        if strcmp(cc,'+'), xrange = xrange/2; set(gcf,'CurrentCharacter','a'); end
        if strcmp(cc,'-'), xrange = xrange*2; set(gcf,'CurrentCharacter','a'); end
    end
end
%%
if update_flag, fclose(fid); if fid1 ~= -1, fclose(fid1); end; end

%%
% v=data(vars.si.v(2),1:end-1)./diff(t);
% t2 = t(2:end);
% % t2(isnan(v))=[];
% v(isnan(v))=0;
% plot(t2,smooth(v,3))
% % xlim([3000 3500])
%% look at licks triggered on rewards
% rt = vars.worlds{1}.reward_times; % reward times
% L = length(rt); % number of rewards
% b = -100:800; % bin indices
% B = length(b); % number of bins
% h = zeros(L,B); % histogram data
% for i = 1:L
% [~,ti] = min(abs(data(vars.si.t,:)-rt(i))); % find index of reward time
% h(i,:) = data(vars.si.licks,ti+b); % get licks before and after
% end
% imagesc(h)
% % plot(b*vars.dt,mean(h(30:end,:)))

%     %% put everything into output structure
%     analysis.v = v;
%     analysis.licks = licks;
%     analysis.p = p;
%     analysis.y = y;
%     analysis.t = t;
%     analysis.reward = reward;
%     analysis.vs = vs;
    
    
%     Cell_no = genvarname(['cell_' cell_num]);
%     
%     %%% metainformation
%     cells(i).filename = file_name;
%     cells(i).cell_num = Cell_no;
%     
%     %%% raw traces
%     cells(i).data           = data;
%     cells(i).t              = t;
%     cells(i).Fs             = Fs;
%     
%     %%% firing properties
%     cells(i).ISI                    = ISI;
%     cells(i).AI                     = AI;
%     cells(i).burst_idx              = burst_idx;
%     cells(i).spike_times            = spike_times;
%     cells(i).latency_to_first_AP    = latency;
%     cells(i).number_spikes          = number_spikes;
%     cells(i).meanspikefreq          = meanspikefreq;
%     
%     %%% tag comments
%     cells(i).tag_sweeps     = tagSw;
%     cells(i).tag_comments   = tagStr;