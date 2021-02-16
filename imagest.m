function fig = imagst(propertyArray, dataOptions, plotOptions)
% IMAGEST  create an image stack figure out of a 3D array. Slices up the
%          3D array and lays out each slice in an optimized grid. Add
%          masks and contours to make your images tell a story.
%
%    fig = IMAGEST(propertyArray) creates an image stack out of 
%            propertyArray
%    fig = IMAGEST(propertyArray, dataOptions, plotOptions) creates 
%            an image stack using extra dataOptions and plotOptions
%   
%    dataOptions and plotOptions are (name,value) pairs. The function
%    returns a FIGURE object whose 1st child, fig.Children(1), is the
%    COLORBAR, and the remaining children are the AXES. One can fine tune
%    the behavior of the plot by directly accessing the properties within
%    the children. All titles support LaTeX in math mode ex. $\phi$ 
%
%    DATAOPTIONS
%      contourmask      3d logical array that defines where contours are
%                       drawn
%      numSlices        integer that determines how many data slices are
%                       used. Note that if numSlices does not
%                       divide propertyArray's 3rd dimension, then linear
%                       interpolation is used (default: 12)
%      mask             3d logical array used to mask unwanted regions 
%    PLOTOPTIONS
%      background       RGB triplet defining the background color
%      clip             logical value that clips the border based on the 
%                       size of a mask (default: true) 
%      colorbarFontSize positive integer changes the font size in the
%                       colorbar (default : 18)
%      colorbarTitle    string title for the colorbar (default: none)
%      colorbarTitleOffset
%                       offsets the title if the rendering is weird. Sorry
%                       I couldn't predict your every need. 
%      colorMap         colormap used for the figure (default : parula) 
%                       Go Bears!
%      contourColor     RGB triplet for the countour color (default : red)
%      contourLineWidth make yo' contours thicc
%      dim              (r,c) makes the motage stack have r rows and c
%                       columns. If r*c is larger than the number of slices
%                       in the data then they will get optimized. (default
%                       : (3,4))
%      monitor          integer that says the monitor to plot on
%      dataRange        1. assumes data is Gaussian, so will display 
%                        data out to a certain number of standard devs.
%                        choosing a smaller number will give more contrast
%                        at the risk of clipping the data. (default : 3)
%                       2. if you input [a,b] with 0<=a,b<=1, are 
%                        percentages then the color range will be from 
%                        100*a% of the min propValue to 100*b% of the max
%                        propValue.
%      title            string title for the entire montage
%      titleFontSize    change the font size of the title
%
%    Written by yours truly, Damian R. Sowinski, in the year of Corona 2020
%    Last updated - December 20, 2020
%    Questions or comments: drsowinski@gmail.com

  arguments
    propertyArray (:,:,:) double
    dataOptions.mask (:,:,:) {mustBeNumericOrLogical,mustBeInteger,...
      mustBeSameSize(dataOptions.mask,propertyArray)} ...
      = ones(size(propertyArray))
    dataOptions.contourmask (:,:,:) {mustBeNumericOrLogical, ...
      mustBeInteger,mustBeSameSize(dataOptions.contourmask,...
      propertyArray)} = zeros(size(propertyArray))
    dataOptions.numSlices {mustBeInteger,mustBePositive} ...
      = size(propertyArray,3)
    plotOptions.dim (1,2) double {mustBeInteger,mustBePositive} = [3 4]
    plotOptions.colorMap (:,3) double {mustBeNonempty,mustBeNonnegative,...
      mustBeLessThanOrEqual(plotOptions.colorMap,1)} = parula
    plotOptions.contourColor (1,3) double {mustBeNonempty,mustBeNonnegative,...
      mustBeLessThanOrEqual(plotOptions.contourColor,1)} = [1 0 0]
    plotOptions.contourLineWidth double {mustBePositive} = 1
    plotOptions.title char = 'none'
    plotOptions.titleFontSize double ...
      {mustBePositive,mustBeInteger} = 18
    plotOptions.background (1,3) double {mustBeNonempty, ...
      mustBeNonnegative,mustBeLessThanOrEqual(plotOptions.background,1)}...
      =[1 1 1]
    plotOptions.dataRange double {mustBeNonnegative} = 3
    plotOptions.clip logical = true
    plotOptions.zclip logical = false
    plotOptions.colorbarFontSize double ...
      {mustBeInteger,mustBePositive} = 14
    plotOptions.colorbarTitle char = []
    plotOptions.colorbarTitleOffset double = 1
    plotOptions.monitor double = 1
  end
  
  %turns on the contours if detects a non-default contour map
  if sum(dataOptions.contourmask(:))~=0
    contourTag = true;
  else
    contourTag = false;
  end
  
  %gets rid of slices full of 0s or nans
  
  %clip the 0s off of the mask and array based on the mask
  if plotOptions.clip
    dataOptions.mask(isnan(propertyArray)) = 0;
    dataOptions.mask(propertyArray==0) = 0;
    
    xL = find(diff(sum(dataOptions.mask,[2,3]))~=0,1);
    if isempty(xL)
      xL = 1;
    else
      xL = max(1,xL-1);
    end
    xR = find(diff(sum(dataOptions.mask,[2,3]))~=0,1,'last');
    if isempty(xR)
      xR = size(propertyArray,1);
    else
      xR = min(xR+1,size(propertyArray,1));
    end
    
    yL = find(diff(sum(dataOptions.mask,[1,3]))~=0,1);
    if isempty(yL)
      yL = 1;
    else
      yL = max(1,yL-1);
    end
    yR = find(diff(sum(dataOptions.mask,[1,3]))~=0,1,'last');
    if isempty(yR)
      yR = size(propertyArray,2);
    else
      yR = min(yR+1,size(propertyArray,2));
    end
    
    propertyArray = propertyArray(xL:xR,yL:yR,:);
    dataOptions.mask = dataOptions.mask(xL:xR,yL:yR,:);
    dataOptions.contourmask = dataOptions.contourmask(xL:xR,yL:yR,:);
  end
  
  if plotOptions.zclip
    zL = find(diff(reshape(sum(dataOptions.mask,[1,2]),...
      [1,dataOptions.numSlices]))~=0,1);
    if isempty(zL)
      zL = 1;
    else
      zL = max(1,zL-1);
    end
    zR = find(diff(reshape(sum(dataOptions.mask,[1,2]),...
      [1,dataOptions.numSlices]))~=0,1,'last');
    if isempty(zR)
      zR = size(propertyArray,3);
    else
      zR = min(zR+1,size(propertyArray,3));
    end
    
    propertyArray = propertyArray(:,:,zL:zR);
    dataOptions.mask = dataOptions.mask(:,:,zL:zR);
    dataOptions.contourmask = dataOptions.contourmask(:,:,zL:zR);
    
    dataOptions.numSlices = size(propertyArray,3);
  end
  
  % sets up the number of plots in the montage
  nRow = plotOptions.dim(1);
  nCol = plotOptions.dim(2);
  [imHeight,imWidth,imSlices] = size(propertyArray);
  
  % check to make sure user is not asking to see more slices than exist
  % in the property array.
  if dataOptions.numSlices > imSlices
    dataOptions.numSlices = imSlices;
    %fprintf(['WARNING Data has %d slices. Consider changing DIM or '...
    %  ,'NUMSLICES options.\n'],imSlices)
  end
  deltaSlice = imSlices/(dataOptions.numSlices);
  if dataOptions.numSlices > nRow*nCol
    %fprintf(['WARNING %d by %d plot cannot fit %d slices. More ' ...
    %  ,'subplots added.\n\tConsider changing DIM or NUMSLICES ' ...
    %  ,'options.\n'], nRow,nCol,dataOptions.numSlices)
    while dataOptions.numSlices > (nRow*nCol)
      nCol = nCol + ceil(imHeight/imWidth);
      if dataOptions.numSlices > (nRow*nCol)
        nRow = nRow + ceil(imWidth/imHeight);
      end
    end
  end
  while dataOptions.numSlices <= nRow*(nCol -1)
    nCol = nCol - 1;
  end
  
  %check to make sure the plot size doesn't exceed the screen size,
  % and shrink it if it does.
  MP = get(0, 'MonitorPositions');
  plotWidth  = 200;
  plotHeight = ceil(plotWidth*imHeight/imWidth);
  border     = ceil(max(plotWidth,plotHeight)/4);
  while ((nRow*plotHeight+5*border)>...
      MP(plotOptions.monitor,4))||...
      ((nCol*plotWidth+5.5*border)>...
      MP(plotOptions.monitor,3))
    plotWidth  = plotWidth - 1;
    plotHeight = ceil(plotWidth*imHeight/imWidth);
    border     = ceil(max(plotWidth,plotHeight)/4);
  end
  figHeight = 2*border + nRow*plotHeight + nRow;
  figWidth  = 4*border + nCol*plotWidth + nCol;
  fig       = figure(...
    'position',...
    [MP(plotOptions.monitor,1) MP(plotOptions.monitor,2) ...
      figWidth figHeight],...
    'color',plotOptions.background,...
    'resize','on');
  plotOptions.colorMap(1,:) = plotOptions.background;
  
  %assume data is Gaussian, and get the data limits for plotting. The fewer
  %the number of standard deviations used the more contrast in the color 
  %images, but then you risk clipping the data.
  if length(plotOptions.dataRange) == 1
    propMean = mean(propertyArray(dataOptions.mask==1));
    propSTD  = std(propertyArray(dataOptions.mask==1));
    propMax  = propMean + plotOptions.dataRange*propSTD;
    propMin  = propMean - plotOptions.dataRange*propSTD;
  elseif length(plotOptions.dataRange) == 2
    propMax  = max(propertyArray(dataOptions.mask==1));
    propMin  = min(propertyArray(dataOptions.mask==1));
    propRange = propMax - propMin;
    propMax  = propMin + plotOptions.dataRange(2)*propRange;
    propMin  = propMin + plotOptions.dataRange(1)*propRange;
  end
  
  %generates the axes objects for each of the subplots
  ax = axes(fig,...
    'Position',[0 0 1 1-border/figHeight],...
    'xtick',[],...
    'ytick',[]...
    );
  for r = 1:nRow
    for c = 1:nCol
      regNum = c + nCol*(r-1);
      if regNum <= dataOptions.numSlices
        slicePlot(regNum) = axes(fig,...
          'Position',...
          [ (border + (c-1)*plotWidth+c)/figWidth,...
            1 - (border + r*plotHeight+r)/figHeight,...
            plotWidth/figWidth,...
            plotHeight/figHeight ],...
          'colormap',plotOptions.colorMap,...
          'climmode','manual',...
          'clim',[propMin,propMax],...
          'visible','off',...
          'YDir','reverse');
        hold on
        axis tight
      end
    end
  end
  
  %the actual plots get made here
  for regNum = 1:dataOptions.numSlices
    idxCenter = 1+deltaSlice*(regNum-1);
    idxTop    = ceil(idxCenter);
    idxBottom = floor(idxCenter);
    a = idxTop - idxCenter;
    displayData = a*propertyArray(:,:,idxTop) + ...
      (1-a)*propertyArray(:,:,idxBottom);
    displayData = displayData.*dataOptions.mask(:,:,idxTop)...
      .*dataOptions.mask(:,:,idxBottom);
    displayData(dataOptions.mask(:,:,idxTop)...
      .*dataOptions.mask(:,:,idxBottom)==0) = propMin;
    imagesc(slicePlot(regNum),displayData)
    if contourTag
      contourData = ((dataOptions.contourmask(:,:,idxTop) + ...
        dataOptions.contourmask(:,:,idxBottom)) >0)*1.0;
      contour(slicePlot(regNum),contourData,1,...
        'color',plotOptions.contourColor,...
        'linewidth',plotOptions.contourLineWidth)
    end
    axis tight
  end
  
  %this controls the construction of the colorbar. the user can always
  % mess with this later since it is the first child of f. 
  % i.e. f.Children(1) is the colorbar
  cb = colorbar(...
    'Position',...
    [ (2*border+nCol*plotWidth)/figWidth border/figHeight ...
    15/figWidth plotHeight*nRow/figHeight ],...
    'fontsize',plotOptions.colorbarFontSize);
  myExponent = round(...
    log10(mean(abs(propertyArray(logical(dataOptions.mask(:)))))));
  if abs(myExponent)>1
    cb.Ruler.Exponent = myExponent;
  end
  cb.TickLabelInterpreter = 'latex';
  cb.Label.String = plotOptions.colorbarTitle;
  cb.Label.Interpreter = 'latex'; 
  cb.Label.FontSize = plotOptions.colorbarFontSize+2;
  cb.Label.Rotation = -90;
  cb.Label.Position = cb.Label.Position + ...
    [plotOptions.colorbarTitleOffset 0 0];
  
  if ~strcmp(plotOptions.title,'none')
    title(ax,plotOptions.title,...
      'interpreter','latex',...
      'fontsize',plotOptions.titleFontSize)
  end
    
end

function mustBeSameSize(arg,b)
    if any(size(arg)~=size(b))
        error('MASK must be the same size as property array.')
    end
end
