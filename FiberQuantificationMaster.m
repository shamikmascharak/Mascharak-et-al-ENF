%% CJUN Rosa Picro Fiber Segmentation: Skeletonization

clear
clc
close all

% Get list of all TIF files in the directory
imagefiles = dir('*.tif');
nfiles = length(imagefiles);

counter = 0;

% Optional scaling factors (to convert pixels to microns)
scale_1 = 1;
scale_2 = 1;

for i = 1:2:nfiles
    
    redfilename = imagefiles(i).name;
    greenfilename = imagefiles(i+1).name;
    
    stored_red_names{i} = redfilename;
    stored_green_names{i} = greenfilename;
    
    stored_red_im{i} = imread(redfilename);
    stored_green_im{i} = imread(greenfilename);
    
    picro_red = 255-imread(redfilename);
    picro_green = 255-imread(greenfilename);
    
    % Remove noise with adaptive filtering if needed
%      picro_red = wiener2(picro_red,[3 3]);
%     picro_green = wiener2(picro_green,[3 3]);
    
    % Binarize, erode, and dilate
%     picro_red_bw = im2bw(imadjust(picro_red,[0.75 1])); % For overexposed
    picro_red_bw = im2bw(imadjust(picro_red));
    
    %picro_green_bw = im2bw(imadjust(picro_green,[0.15 1])); % For overexposed
    picro_green_bw = im2bw(imadjust(picro_green));
    
    % Define diamond structuring element
    seD = strel('diamond',2);    
    picro_red_bw = imerode(picro_red_bw,seD);
%     picro_red_bw = imerode(picro_red_bw,seD); % If extra erosion needed
%     picro_red_bw = imerode(picro_red_bw,seD); % If extra erosion needed
    
%     seD = strel('diamond',1); % Erode green images if noisy
%     picro_green_bw = imerode(picro_green_bw,seD);
    
    % Close bw objects
    picro_red_bw = imclose(picro_red_bw,seD);
    picro_red_bw = imclose(picro_red_bw,seD);
    
    picro_green_bw = imclose(picro_green_bw,seD);
    
    % Skeletonize, separate at branchpoints, and delete single pixel lines
    picro_red_skel = bwmorph(picro_red_bw,'skel',Inf);
    picro_red_branch = bwmorph(picro_red_skel,'branchpoints');
    picro_red_skel(picro_red_branch == 1) = 0;
    picro_red_skel = bwmorph(picro_red_skel,'clean');
    
    stored_red_skel{i} = picro_red_skel;
    
    [L_red, num_red] = bwlabel(picro_red_skel);
    [~, num_red_branch] = bwlabel(picro_red_branch);
    
    props_red = regionprops(picro_red_bw,'Area','EulerNumber','Extent','ConvexArea',...
        'FilledArea','Solidity','Perimeter','Eccentricity','MajorAxisLength',...
        'EquivDiameter','MinorAxisLength');
    props_red_skel = regionprops(picro_red_skel,'Orientation');
    kappa_red = circ_kappa([props_red_skel.Orientation]*pi/180);
    
    stored_red(i,:) = [mean(mean(picro_red)) num_red/scale_2 mean([props_red.MajorAxisLength])*scale_1 ...
        mean([props_red.MinorAxisLength])*scale_1 mean([props_red.MajorAxisLength])/mean([props_red.Perimeter]) ...
        kappa_red num_red_branch/scale_2 mean([props_red.EulerNumber]) mean([props_red.Extent]) ...
        mean([props_red.Perimeter])*scale_1 mean([props_red.Solidity]) mean([props_red.Eccentricity]) ...
        mean([props_red.EquivDiameter])*scale_1];
    
    picro_green_skel = bwmorph(picro_green_bw,'skel',Inf);
    picro_green_branch = bwmorph(picro_green_skel,'branchpoints');
    picro_green_skel(picro_green_branch == 1) = 0;
    picro_green_skel = bwmorph(picro_green_skel,'clean');
    
    stored_green_skel{i} = picro_green_skel;
    
    [L_green, num_green] = bwlabel(picro_green_skel);
    [~, num_green_branch] = bwlabel(picro_green_branch);
    
    props_green = regionprops(picro_green_bw,'Area','EulerNumber','Extent','ConvexArea',...
        'FilledArea','Solidity','Perimeter','Eccentricity','MajorAxisLength',...
        'EquivDiameter','MinorAxisLength');
    props_green_skel = regionprops(picro_green_skel,'Orientation');
    kappa_green = circ_kappa([props_green_skel.Orientation]*pi/180);
    
    stored_green(i,:) = [mean(mean(picro_green)) num_green/scale_2 mean([props_green.MajorAxisLength])*scale_1 ...
        mean([props_green.MinorAxisLength])*scale_1 mean([props_green.MajorAxisLength])/mean([props_green.Perimeter]) ...
        kappa_green num_green_branch/scale_2 mean([props_green.EulerNumber]) mean([props_green.Extent]) ...
        mean([props_green.Perimeter])*scale_1 mean([props_green.Solidity]) mean([props_green.Eccentricity]) ...
        mean([props_green.EquivDiameter])*scale_1];
    
    % Brightness, Number of Fibers/1000 microns^2, Length (micron), Width
    % (micron), Persistence, Angle Randomness, Number of Branchpoints/1000
    % micron^2, Area, Euler Number, Extent, Perimeter, Convex Area, Filled
    % Area, Solidity, Eccentricity, Equivalent Diameter
    
    counter = counter + 1;
    progress = 2*100*counter/(nfiles)
    
end

stored_red(all(stored_red == 0,2),:) = []; % Fiber characteristics for red images
stored_red_names = stored_red_names(~cellfun(@isempty,stored_red_names))';
stored_red_skel = stored_red_skel(~cellfun(@isempty,stored_red_skel))';


stored_green(all(stored_green == 0,2),:) = []; % Fiber characteristics for green images
stored_green_names = stored_green_names(~cellfun(@isempty,stored_green_names))';
stored_green_skel = stored_green_skel(~cellfun(@isempty,stored_green_skel))';

stored_redgreen = [stored_red stored_green];

%% Define tSNE class labels

close all

% 2 week images
class_pt = {'PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS',...
    'PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS',...
    'PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS',...
    'PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS',...
    'PBS','PBS','PBS','PBS','PBS','PBS','PBS',...
    'Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded',...
    'Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded',...
    'Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded',...
    'Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded',...
    'Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded',...
    'Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded',...
    'Vert','Vert','Vert','Vert','Vert','Vert','Vert','Vert','Vert','Vert',...
    'Vert','Vert','Vert','Vert','Vert','Vert','Vert','Vert','Vert','Vert',...
    'Vert','Vert','Vert','Vert','Vert','Vert','Vert','Vert','Vert','Vert',...
    'Vert','Vert','Vert','Vert','Vert','Vert','Vert','Vert','Vert','Vert'};

% 1 month images
class_pt = {'PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS',...
    'PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS',...
    'PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS',...
    'PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS',...
    'PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS',...
    'Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded',...
    'Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded',...
    'Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded',...
    'Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded',...
    'Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded',...
    'Vert','Vert','Vert','Vert','Vert','Vert','Vert','Vert','Vert','Vert',...
    'Vert','Vert','Vert','Vert','Vert','Vert','Vert','Vert','Vert','Vert',...
    'Vert','Vert','Vert','Vert','Vert','Vert','Vert','Vert','Vert','Vert',...
    'Vert','Vert','Vert','Vert','Vert','Vert','Vert','Vert','Vert','Vert',...
    'Vert','Vert','Vert','Vert'};

% 3 month images
class_pt = {'PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS',...
    'PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS',...
    'PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS',...
    'Unwounded','Unwounded','Unwounded','Unwounded','Unwounded',...
    'Unwounded','Unwounded','Unwounded','Unwounded','Unwounded',...
    'Unwounded','Unwounded','Unwounded','Unwounded','Unwounded',...
    'Unwounded','Unwounded','Unwounded','Unwounded','Unwounded',...
    'Unwounded','Unwounded','Unwounded','Unwounded','Unwounded',...
    'Unwounded','Unwounded','Unwounded','Unwounded','Unwounded',...
    'Unwounded','Unwounded','Unwounded','Unwounded','Unwounded',...
    'Unwounded','Unwounded','Unwounded','Unwounded','Unwounded',...
    'Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded','Unwounded',...
    'Vert','Vert','Vert','Vert','Vert','Vert','Vert','Vert','Vert','Vert',...
    'Vert','Vert','Vert','Vert','Vert','Vert','Vert','Vert','Vert','Vert',...
    'Vert','Vert','Vert','Vert','Vert','Vert','Vert','Vert','Vert','Vert',...
    'Vert','Vert','Vert','Vert','Vert','Vert',};

% Class labels for MD Vert experiment
class_pt = {'PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS',...
    'PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS',...
    'PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS',...
    'PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS',...
    'PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS',...
    'Unwounded','Unwounded','Unwounded','Unwounded','Unwounded',...
    'Unwounded','Unwounded','Unwounded','Unwounded','Unwounded',...
    'Unwounded','Unwounded','Unwounded','Unwounded','Unwounded',...
    'Unwounded','Unwounded','Unwounded','Unwounded','Unwounded',...
    'Unwounded','Unwounded','Unwounded','Unwounded','Unwounded',...
    'Unwounded','Unwounded','Unwounded','Unwounded','Unwounded',...
    'Unwounded','Unwounded','Unwounded','Unwounded','Unwounded',...
    'Unwounded','Unwounded','Unwounded','Unwounded','Unwounded',...
    'Unwounded','Unwounded','Unwounded','Unwounded','Unwounded',...
    'Unwounded','Unwounded','Unwounded','Unwounded','Unwounded',...
    'Unwounded','Unwounded','Unwounded','Unwounded','Unwounded',...
    'Unwounded','Unwounded','Unwounded','Unwounded','Unwounded',...
    'Unwounded','Unwounded','Unwounded','Unwounded','Unwounded',...
    'Unwounded','Unwounded','Unwounded','Unwounded','Unwounded',...
    'Unwounded','Unwounded','Unwounded','Unwounded','Unwounded',...
    'Unwounded','Unwounded','Unwounded','Unwounded','Unwounded',...
    'Unwounded','Unwounded','Unwounded','Unwounded','Unwounded',...
    'Unwounded','Unwounded','Unwounded','Unwounded','Unwounded',...
    'Unwounded','Unwounded','Unwounded','Unwounded','Unwounded',...
    'Unwounded','Unwounded','Unwounded',...
    '2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month',...
    '2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month',...
    '2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month',...
    '2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month',...
    '2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month',...
    '2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month',...
    '2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month',...
    '2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month',...
    '2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month',...
    '2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month',...
    '2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month',...
    '2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month',...
    '2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month',...
    '4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month',...
    '4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month',...
    '4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month',...
    '4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month',...
    '4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month',...
    '4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month',...
    '4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month',...
    '4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month',...
    '4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month',...
    '4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month',...
    '4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month',...
    '4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month',...
    '4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month',...
    '4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks',...
    '4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks',...
    '4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks',...
    '4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks',...
    '4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks',...
    '4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks',...
    '4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks',...
    '4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks',...
    '4 Doses - 2 weeks','4 Doses - 2 weeks',...
    '2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks',...
    '2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks',...
    '2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks',...
    '2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks',...
    '2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks',...
    '2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks',...
    '2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks',...
    '2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks',...
    '2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks',...
    '2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks',...
    '2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks',...
    '2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks',...
    '2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks',...
    '2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks',...
    '2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks',...
    '2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks',...
    '2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks',...
    '2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks',...
    '2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks',...
    '2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks','2 Doses - 2 weeks',...
    '2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month',...
    '2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month',...
    '2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month',...
    '2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month',...
    '2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month',...
    '2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month',...
    '2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month',...
    '2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month',...
    '2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month',...
    '2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month',...
    '2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month',...
    '2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month',...
    '2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month','2 Doses - 1 month',...
    '2 Doses - 2 weeks',...
    'PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS',...
    'PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS',...
    'PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS',...
    'PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS','PBS',...
    'PBS','PBS',...
    '4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks',...
    '4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks',...
    '4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks',...
    '4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks',...
    '4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks',...
    '4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks',...
    '4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks',...
    '4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks',...
    '4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks',...
    '4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks',...
    '4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks',...
    '4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks',...
    '4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks',...
    '4 Doses - 2 weeks','4 Doses - 2 weeks','4 Doses - 2 weeks',...
    '4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month',...
    '4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month',...
    '4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month',...
    '4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month',...
    '4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month',...
    '4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month',...
    '4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month',...
    '4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month',...
    '4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month',...
    '4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month',...
    '4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month','4 Doses - 1 month'};

Y = tsne(stored_redgreen);
figure
gscatter(Y(:,1),Y(:,2),class_pt')
