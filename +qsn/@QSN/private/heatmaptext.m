function  heatmaptext(data,varargin)
  %HEATMAPTEXT creates a heatmap with the values shown as text.
  %
  %   HEATMAPTEXT(SEQ) counts the number of occurrences of each codon in the
  %   sequence and displays a formatted table of the result.
  %
  %   HEATMAPTEXT(...,'PRECISION',P) displays the data values with a maximum N
  %   digits of precision.  The default number of digits of precision is 2.
  %
  %   HEATMAPTEXT(SEQ, ... ,'COLORBAR',false) hides the colorbar.
  %
  %   HEATMAPTEXT(SWQ, ... ,'CLIM',[min max]) sets range
  %
  %   HEATMAPTEXT(SWQ, ... ,'CMAP',cmap) sets colormap
  %
  %
  %   Examples:
  %
  %       ih = gallery('invhess',10);
  %       heatmaptext(ih);
  %
  %       figure;
  %       m = gallery('moler',30);
  %       heatmaptext(m);
  %
  %       figure;
  %       heatmaptext(rand(20,10),'fontcolor','r','precision',4);
  %       colormap(bone);
  %
  %   See also COLORMAP, IMAGESC, TEXT.

  %   Copyright 2007 The MathWorks, Inc.
  %   $Revision: 1.1 $  $Date: 2007/08/08 21:56:22 $

  precision = 2;
  cBar = true;
  clim = NaN;
  cmap = colormap();
  if nargin > 1
    if rem(nargin,2)== 0
      error('heatmaptext:IncorrectNumberOfArguments',...
            'Incorrect number of arguments to %s.',mfilename);
    end
    okargs = {'precision', 'colorbar', 'clim', 'cmap'};
    for j=1:2:nargin-2
      pname = varargin{j};
      pval = varargin{j+1};
      k = strmatch(lower(pname), okargs);
      if isempty(k)
        error('heatmaptext:UnknownParameterName',...
              'Unknown parameter name: %s.',pname);
      elseif length(k)>1
        error('heatmaptext:AmbiguousParameterName',...
              'Ambiguous parameter name: %s.',pname);
      else
        switch(k)
          case 1  % frame
            precision = pval;
          case 2
            cBar = pval;
          case 3
            clim = pval;
%          vsdr 4
          case 4
            cmssp = pval;
        end
      end
    end
  end

  if isnan (clim)
    im = imagesc(data);
  else
    im = imagesc(data,clim);
  end
  imAxes =get(im,'parent');
  hFig = get(imAxes,'parent');
  fs = getBestFontSize(imAxes);
  showText = true;
  if fs == 0
    showText = false;
    fs = 9;
  end
  axis off;
  if cBar
    colorbar;
  end
  [rows, cols] = size(data);

  midValue = mean(get(gca,'CLim'));
  ci = (data < midValue) + 1;
  [mx my] = size(cmap);
  cmap = [cmap(1,:); cmap(mx,:)];

  textHandles = zeros(size(data))';
  for i = 1:rows
    for j = 1:cols
      textHandles(j,i) = text(j,i,num2str(data(i,j),precision),...
                              'color',cmap(ci(i,j),:),...
                              'horizontalAlignment','center','verticalAlignment','middle',...
                              'fontsize',fs,'clipping','on','visible','off');
    end
  end
  if showText
    set(textHandles,'visible','on');
  end
  set(imAxes,'UserData',textHandles);

  function fs = getBestFontSize(imAxes)
    % Try to keep font size reasonable for text
    hFig = get(imAxes,'Parent');
    magicNumber = 80;
    nrows = diff(get(imAxes,'YLim'));
    ncols = diff(get(imAxes,'XLim'));
    if ncols < magicNumber && nrows < magicNumber
      ratio = max(get(hFig,'Position').*[0 0 0 1])/max(nrows,ncols);
    elseif ncols < magicNumber
      ratio = max(get(hFig,'Position').*[0 0 0 1])/ncols;
    elseif nrows < magicNumber
      ratio = max(get(hFig,'Position').*[0 0 0 1])/nrows;
    else
      ratio = 1;
    end
    fs = min(9,ceil(ratio/4));    % the gold formula
    if fs < 4
      fs = 0;
    end
  end

end
