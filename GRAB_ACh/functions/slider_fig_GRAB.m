function [threshold,regions] = slider_fig_GRAB(img)

H.img = img;                                % input image
prct = 20;                                  % upper N percent of signal used for detection region (defauly 20%) 
H.thresh = min(H.img(:)) + (max(H.img(:))-min(H.img(:)))*(prct/100); % initial threshold 

H.fig = figure;                             % figure
H.topLabel = annotation('textbox',[0.18,0.87,0.8,0.1],...
    'string','Set the threshold for GRAB-ACh3.0 detection','fontsize',12,'fontweight','bold','LineStyle','none');

% --- display input image --- %
H.ax1 = axes('Parent',H.fig,'Units','Normalized','Position',[0.05 0.35 0.5 0.5]); % axis 1
H.plt1 = imagesc(H.img);                    % plot 1
colormap(hot)
colorbar
axis image
axis off
title('Mean projection')
hold on

% --- display thresholded area on input image --- % 
[H.rowPix,H.colPix] = size(H.img);          % number of row and column pixels
H.B = bwboundaries(H.img>H.thresh);         % blob detection (thresholded areas
H.blobAx = false(H.rowPix,H.colPix);        % pre-allocated blob axis
for k = 1:numel(H.B)                        % cycle through all thresholded areas
    H.blobInd = sub2ind([H.rowPix H.colPix],H.B{k}(:,1),H.B{k}(:,2)); % integer index of kth blob perimeter 
    H.blobAx(H.blobInd) = true;             % set kth blob perimeter values as true on blob axis
end
H.overlay = repmat(ones(H.rowPix,H.colPix),[1 1 3]); % overlay: white image
H.lines = imshow(H.overlay);                % display overlay
set(H.lines,'alphadata',H.blobAx)           % set overlay as blob perimeters

% --- display thresholded area as binary image --- %   
H.ax2 = axes('Parent',H.fig,'Units','Normalized','Position',[0.5 0.35 0.5 0.5]);
H.plt2 = imagesc(H.img>H.thresh);           % binarised image
axis image
axis off
title('Thresholded')

% --- make thresholding slider --- %   
H.loThresh = min(H.img(:));                 % slider minimum value 
H.hiThresh = max(H.img(:));                 % slider maximum value
H.slider1 = uicontrol('Parent',H.fig,...    
                      'Units','Normalized',...
                      'Position',[0.2 0.18 0.6 0.05],...
                      'Style','Slider',...
                      'BackgroundColor',[1 1 1],...
                      'Min',H.loThresh,'Max',H.hiThresh,...
                      'Value',H.thresh,...
                      'Callback',@slider1Callback);
sliderLabel = 'Move slider to adjust image threshold';                                                           % text label for slider
H.mainLabel = annotation('textbox',[0.28,0.20,0.6,0.1],'string',sliderLabel,'LineStyle','none');                 % display main slider label
H.loLabel = annotation('textbox',[0.17,0.08,0.2,0.1],'string',num2str(round(H.loThresh,2)),'LineStyle','none');  % display slider low threshold label    
H.hiLabel = annotation('textbox',[0.75,0.08,0.2,0.1],'string',num2str(round(H.hiThresh,2)),'LineStyle','none');  % display slider high threshold label
H.midThresh = H.loThresh+((H.hiThresh-H.loThresh)/2);
H.midLabel = annotation('textbox',[0.445,0.08,0.2,0.1],'string',num2str(round(H.midThresh,2)),'LineStyle','none');

% --- make OK button --- % 
H.okButton = uicontrol('Parent',H.fig,...    
                      'Units','Normalized',...
                      'Position',[0.87 0.03 0.07 0.07],...
                      'Style', 'Pushbutton',...
                      'Backgroundcolor', [0 1 0],... 
                      'Userdata', 0,... 
                      'String', 'OK',... 
                      'Fontsize', 12,...
                      'Fontweight','bold',...
                      'Callback', @okButtonCallback);  % OK button 

% --- wait for user to interact with GUI --- %                
guidata(H.fig,H)                            % set GUI data as H object   
waitfor(H.okButton,'UserData',1)            % wait until OK button is pressed to continue 

% --- generate outputs --- % 
guidata(H.fig,H)                            % latest GUI data from figure
threshold = H.slider1.Value;                % get slider value
regions = H.B;                              % get final blob values

close(H.fig)                                % close figure    

end

% ---- callback function for slider ---- %    
function slider1Callback(hObject,eventdata)
    H = guidata(hObject);
    
    % --- Binarised image --- %
    frame = H.plt1.CData;                   % extract colour data from displayed input image
    H.plt2.CData = frame>H.slider1.Value;   % re-threshold plot 2 (binarised image) based on slider position  

    % --- Input image overlay lines --- %
    H.B = bwboundaries(H.img>H.slider1.Value); % redo blob analysis based on slider position
    H.blobAx = false(H.rowPix,H.colPix);    % redefine blob plot
    for j = 1:numel(H.B)                    % cycle through new blobs
        H.blobInd = sub2ind([H.rowPix H.colPix],H.B{j}(:,1),H.B{j}(:,2)); % integer indices of jth blob perimeter
        H.blobAx(H.blobInd) = true;         % add jth blob perimeter to blob plot
    end
    set(H.lines,'alphadata',H.blobAx)       % reset overlay as updated blob plot
end

% ---- callback function for OK button ---- %  
function okButtonCallback(source, eventdata)
    % Callback function for 'OK' button
    set(source, 'userdata', 1);             % when press button set 
    uiresume
end

