function [roiMask, roiX, roiY] = drawROI_GRAB(img)  

%% Make gui figure
scrSz = get(0, 'screensize');                   % find screen dimensions in pixels
fHan = figure('units', 'pixels', 'position', [scrSz(3)*0.2 scrSz(4)*0.2 scrSz(3)*0.6 scrSz(4)*0.6],...
    'visible', 'off');                          % make gui figure
movegui(fHan, 'center')                         % move gui to centre of screen
fSz = get(fHan, 'position');                    % get figure dimensions in pixels

%% Make pushbuttons
okButton = uicontrol('style', 'pushbutton', 'position', [fSz(3)*0.75 fSz(4)*0.5 fSz(3)*0.1 fSz(4)*0.05],...
    'backgroundcolor', [0 1 0], 'userdata', 0, 'string', 'OK', 'fontname', 'arial', 'fontsize', 10, 'fontweight', 'bold',...
    'callback', @okButtonCallback);             % OK button 

resetButton = uicontrol('style', 'pushbutton', 'position', [fSz(3)*0.75 fSz(4)*0.4 fSz(3)*0.1 fSz(4)*0.05],...
    'backgroundcolor', [1 0.5 0], 'userdata', 0, 'string', 'Reset', 'fontname', 'arial', 'fontsize', 10, 'fontweight', 'bold',...
    'callback', @resetButtonCallback);          % reset button 

cancelButton = uicontrol('style', 'pushbutton', 'position', [fSz(3)*0.75 fSz(4)*0.3 fSz(3)*0.1 fSz(4)*0.05],...
    'backgroundcolor', [1 0 0], 'userdata', 0, 'string', 'Cancel', 'fontname', 'arial', 'fontsize', 10, 'fontweight', 'bold',...
    'callback', @cancelButtonCallback);         % cancel button 

%% Write  text instructions
annotation('textbox',[0.12,0.85,0.8,0.1],...
    'string','Define ROI for background subtraction','fontsize',12,'fontweight','bold','LineStyle','none'); % instruction label

txtInfo = uicontrol('style', 'text', 'position', [fSz(3)*0.65 fSz(4)*0.55 fSz(3)*0.3 fSz(4)*0.2],...
    'fontname', 'arial', 'fontsize', 11, 'horizontalAlignment', 'left',...
    'string', ['(1) Manually draw ROI on figure' newline '(2) Double-click to finish']);

%% Plot image
axHan = axes('units', 'pixels', 'position', [fSz(3)*0.1 fSz(4)*0.05 fSz(4)*0.9 fSz(4)*0.9]);  % create axes in gui
imagesc(img)                                    % plot figutre
colormap(hot)                                   % alter colormap 
colorbar
axis image
set(axHan, 'xtick', [], 'ytick', [])            % remove ticks from axes

%% Make gui visible
set(fHan, 'visible', 'on')

%% Use recursive function to repeatedly define the ROI until OKed by user input
while true    
    [roiMask, roiX, roiY] = roipoly;                    % define the roi
    
    maskOverlay = alphamask(roiMask, [1 1 1], 0.3); axis(gca, 'square', 'tight')  % overlay the roi mask on the image
    lineOverlay = line(roiX, roiY, 'color', [1 1 1]);  % show the boundary of the roi
    
    set(txtInfo, 'string', ['ROI done!' newline '(1) Press OK to continue',...
        newline '(2) Press Reset to redo' newline '(3) Press Cancel to to cancel'], 'foregroundcolor', 'red') % update the text instructions      
        
    uiwait(fHan)                                        % pause function execution until user presses one of gui push buttons
    
    if get(okButton, 'userdata') == 1                   % if press okButton
        break                                           % terminate while loop  
    elseif get(resetButton, 'userdata') == 1            % if press resetButton
        clear('roiMask', 'roiX', 'roiY')                % clear roipoly outputs  
        delete(maskOverlay)                             % delete mask overlay 
        delete(lineOverlay)                             % delete line overlays 
        set(txtInfo, 'string', ['(1) Manually draw ROI on figure' newline '(2) Double-click to finish'],...
            'foregroundcolor', 'black');                % reset text instructions     
        set(resetButton, 'userdata', 0);                % reset resetButton userdata back to 0
    elseif get(cancelButton, 'userdata') == 1           % if press cancelButton
        roiMask = []; roiX = []; roiY = [];             % output roipoly variables as empty matrices
        break                                           % terminate while loop
    end
end

close(fHan)                                             % close gui

end

%% Callback functions
function okButtonCallback(source, eventdata)
% Callback function for 'OK' button
set(source, 'userdata', 1);                             % when press button set 
uiresume
end

function resetButtonCallback(source, eventdata)
% Callback function for 'Reset' button
set(source, 'userdata', 1);                             % when press button set 
uiresume
end

function cancelButtonCallback(source, eventdata)
% Callback function for 'cancel' button
set(source, 'userdata', 1);                             % when press button set 
uiresume
end

