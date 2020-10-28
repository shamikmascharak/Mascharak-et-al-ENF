%% Fiber Quantification Algorithm v2
% Created by Shamik Mascharak
% November 28th, 2019 - Created
% December 1st, 2019 - First iteration completed
% October 28th, 2020 - Uploaded to Github

%% Read files and generate binarized images

clear
clc
close all

% Get list of all .tif files in the directory
% imagefiles = dir('*.tif');
imagefiles = dir('*.png');
nfiles = length(imagefiles);

%% Run color deconvolution

counter = 0;

% Optional scaling factors (to convert pixels to microns)
scale = 1;

% Define target image for image normalization
target_im = imread(imagefiles(1544).name); % Unwounded Extra 8

verbose = 0; % 1 to display results, 0 to not

% Define stain vector for color deconvolution
picro_vect = [0.001 0.70710677 0.70710677;...
    0.70710677 0.001 0.70710677;...
    0.57735026 0.57735026 0.57735026];

for i = 1:nfiles
    filename = imagefiles(i).name;
    stored_names{i} = filename;
    stored_im{i} = imread(filename); % Store the original image
    
    % Normalize image using RGB histogram method
    stored_norm_im{i} = Norm(stored_im{i},target_im,'RGBHist',verbose);
%     stored_norm_im{i} = stored_im{i}; % To bypass normalization

    % Perform color deconvolution (Ruifrok method)
    [Dch M] = Deconvolve(stored_norm_im{i},picro_vect,0);
    
    [red green black] = PseudoColourStains(Dch,M);
    
     picro_red = (255-rgb2gray(red))-(255-rgb2gray(black));
     picro_green = (255-rgb2gray(green))-(255-rgb2gray(black));
    
    % Remove noise with adaptive filtering
    picro_red_noise = wiener2(picro_red,[3 3]);
    picro_green_noise = wiener2(picro_green,[3 3]);
    
    stored_picro_red_noise{i} = picro_red_noise;
    stored_picro_green_noise{i} = picro_green_noise;
    
    picro_red_bw = im2bw(imadjust(picro_red_noise));
    picro_green_bw = im2bw(imadjust(picro_green_noise));
    
    % Define diamond structuring element
    seD = strel('diamond',1);
    picro_red_bw_strel = imerode(picro_red_bw,seD);
    picro_green_bw_strel = imerode(picro_green_bw,seD);
    
    stored_picro_red_strel{i} = picro_red_bw_strel; % Binary map
    stored_picro_green_strel{i} = picro_green_bw_strel; % Binary map
    
    % Skeletonize, separate at branchpoints, and delete single pixel lines
    picro_red_bw_skel = bwmorph(picro_red_bw_strel,'skel',Inf);
    picro_red_bw_branch = bwmorph(picro_red_bw_skel,'branchpoints');
    picro_red_bw_skel(picro_red_bw_branch == 1) = 0;
    picro_red_bw_skel = bwmorph(picro_red_bw_skel,'clean');
    
    stored_picro_red_skel{i} = picro_red_bw_skel; % Skeletonized image
    stored_picro_red_branch{i} = picro_red_bw_branch;
    
    picro_green_bw_skel = bwmorph(picro_green_bw_strel,'skel',Inf);
    picro_green_bw_branch = bwmorph(picro_green_bw_skel,'branchpoints');
    picro_green_bw_skel(picro_green_bw_branch == 1) = 0;
    picro_green_bw_skel = bwmorph(picro_green_bw_skel,'clean');
    
    stored_picro_green_skel{i} = picro_green_bw_skel; % Skeletonized image
    stored_picro_green_branch{i} = picro_green_bw_branch;
    
    counter = counter + 1;
    progress = 100*counter/(nfiles)
    
end

stored_names = stored_names';

%% Quantify image features

clc
counter = 0;

for i = 1:nfiles
    
    red = stored_picro_red_noise{i};
    green = stored_picro_green_noise{i};
    red_bw = stored_picro_red_strel{i};
    green_bw = stored_picro_red_strel{i};
    red_skel = stored_picro_red_skel{i};
    green_skel = stored_picro_green_skel{i};
    red_branch = stored_picro_red_branch{i};
    green_branch = stored_picro_green_branch{i};
    
    % Calculate values from grayscale images
    props_red = regionprops(red_bw,red,'all');
    props_green = regionprops(green_bw,green,'all');
    % Fiber angle-related values
    alpha_red = [props_red.Orientation]*pi/180;
    circ_red = [mean(alpha_red) median(alpha_red) std(alpha_red) skewness(alpha_red)...
        kurtosis(alpha_red) circ_kappa(alpha_red)];

    alpha_green = [props_green.Orientation]*pi/180;
    circ_green = [mean(alpha_green) median(alpha_green) std(alpha_green) skewness(alpha_green)...
        kurtosis(alpha_green) circ_kappa(alpha_green)];

    % Calculate Haralick features from grayscale images (4 different offsets)
    offset = [0 1; -1 1;-1 0;-1 -1];
    glcm_red = graycomatrix(red,'Offset',offset,'Symmetric',true);
    glcm_green = graycomatrix(green,'Offset',offset,'Symmetric',true);
    
    grayprops_red = graycoprops(glcm_red);
    grayprops_green = graycoprops(glcm_green);
    
    % Calculate values from binary images
    props_red_bw = regionprops(red_bw,'all');
    props_green_bw = regionprops(green_bw,'all');
    % Fiber angle-related values
    alpha_red_bw = [props_red_bw.Orientation]*pi/180;
    circ_red_bw = [mean(alpha_red_bw) median(alpha_red_bw) std(alpha_red_bw) skewness(alpha_red_bw)...
        kurtosis(alpha_red_bw) circ_kappa(alpha_red_bw)];
    
    alpha_green_bw = [props_green_bw.Orientation]*pi/180;
    circ_green_bw = [mean(alpha_green_bw) median(alpha_green_bw) std(alpha_green_bw) skewness(alpha_green_bw)...
        kurtosis(alpha_green_bw) circ_kappa(alpha_green_bw)];

    % Calculate values from skeletonized images
    props_red_skel = regionprops(red_skel,'all');
    props_green_skel = regionprops(green_skel,'all');
    % Fiber angle-related values
    alpha_red_skel = [props_red_skel.Orientation]*pi/180;
    circ_red_skel = [mean(alpha_red_skel) median(alpha_red_skel) std(alpha_red_skel) skewness(alpha_red_skel)...
        kurtosis(alpha_red_skel) circ_kappa(alpha_red_skel)];

    alpha_green_skel = [props_green_skel.Orientation]*pi/180;
    circ_green_skel = [mean(alpha_green_skel) median(alpha_green_skel) std(alpha_green_skel) skewness(alpha_green_skel)...
        kurtosis(alpha_green_skel) circ_kappa(alpha_green_skel)];
    
    % Branchpoint-related values
    [L_red, num_red] = bwlabel(red_skel);
    [~, num_red_branch] = bwlabel(red_branch);
    [L_green, num_green] = bwlabel(green_skel);
    [~, num_green_branch] = bwlabel(green_branch);
    
    % Grayscale regionprops, grayscale circ, Haralick features, binary regionprops,
    % binary circ, skeletonized regionprops, skeletonized circ, branchpoints,
    
    quantified{i} = [mean([props_red.Area]) std([props_red.Area]) mean([props_green.Area]) std([props_green.Area])...
        mean([props_red.MajorAxisLength]) std([props_red.MajorAxisLength]) mean([props_green.MajorAxisLength]) std([props_green.MajorAxisLength])...
        mean([props_red.MinorAxisLength]) std([props_red.MinorAxisLength]) mean([props_green.MinorAxisLength]) std([props_green.MinorAxisLength])...
        mean([props_red.Eccentricity]) std([props_red.Eccentricity]) mean([props_green.Eccentricity]) std([props_green.Eccentricity])...
        mean([props_red.ConvexArea]) std([props_red.ConvexArea]) mean([props_green.ConvexArea]) std([props_green.ConvexArea])...
        mean([props_red.Circularity]~=Inf) std([props_red.Circularity]~=Inf) mean([props_green.Circularity]~=Inf) std([props_green.Circularity]~=Inf)...
        mean([props_red.FilledArea]) std([props_red.FilledArea]) mean([props_green.FilledArea]) std([props_green.FilledArea])...
        mean([props_red.EulerNumber]) std([props_red.EulerNumber]) mean([props_green.EulerNumber]) std([props_green.EulerNumber])...
        sum([props_red.EulerNumber]) sum([props_green.EulerNumber])...
        mean([props_red.EquivDiameter]) std([props_red.EquivDiameter]) mean([props_green.EquivDiameter]) std([props_green.EquivDiameter])...
        mean([props_red.Solidity]) std([props_red.Solidity]) mean([props_green.Solidity]) std([props_green.Solidity])...
        mean([props_red.Extent]) std([props_red.Extent]) mean([props_green.Extent]) std([props_green.Extent])...
        mean([props_red.Perimeter]) std([props_red.Perimeter]) mean([props_green.Perimeter]) std([props_green.Perimeter])...
        mean([props_red.PerimeterOld]) std([props_red.PerimeterOld]) mean([props_green.PerimeterOld]) std([props_green.PerimeterOld])...
        mean([props_red.MeanIntensity]) std([props_red.MeanIntensity]) mean([props_green.MeanIntensity]) std([props_green.MeanIntensity])...
        mean([props_red.MinIntensity]) std(double([props_red.MinIntensity])) mean([props_green.MinIntensity]) std(double([props_green.MinIntensity]))...
        mean([props_red.MaxIntensity]) std(double([props_red.MaxIntensity])) mean([props_green.MaxIntensity]) std(double([props_green.MaxIntensity]))...
        mean([props_red.MaxFeretDiameter]) std([props_red.MaxFeretDiameter]) mean([props_green.MaxFeretDiameter]) std([props_green.MaxFeretDiameter])...
        mean([props_red.MaxFeretAngle]) std([props_red.MaxFeretAngle]) mean([props_green.MaxFeretAngle]) std([props_green.MaxFeretAngle])...
        mean([props_red.MinFeretDiameter]) std([props_red.MinFeretDiameter]) mean([props_green.MinFeretDiameter]) std([props_green.MinFeretDiameter])...
        mean([props_red.MinFeretAngle]) std([props_red.MinFeretAngle]) mean([props_green.MinFeretAngle]) std([props_green.MinFeretAngle])...
        circ_red circ_green...
        [grayprops_red.Contrast] [grayprops_red.Correlation] [grayprops_red.Energy] [grayprops_red.Homogeneity]...
        [grayprops_green.Contrast] [grayprops_green.Correlation] [grayprops_green.Energy] [grayprops_green.Homogeneity]...
        mean([props_red_bw.Area]) std([props_red_bw.Area]) mean([props_green_bw.Area]) std([props_green_bw.Area])...
        mean([props_red_bw.MajorAxisLength]) std([props_red_bw.MajorAxisLength]) mean([props_green_bw.MajorAxisLength]) std([props_green_bw.MajorAxisLength])...
        mean([props_red_bw.MinorAxisLength]) std([props_red_bw.MinorAxisLength]) mean([props_green_bw.MinorAxisLength]) std([props_green_bw.MinorAxisLength])...
        mean([props_red_bw.Eccentricity]) std([props_red_bw.Eccentricity]) mean([props_green_bw.Eccentricity]) std([props_green_bw.Eccentricity])...
        mean([props_red_bw.ConvexArea]) std([props_red_bw.ConvexArea]) mean([props_green_bw.ConvexArea]) std([props_green_bw.ConvexArea])...
        mean([props_red_bw.Circularity]~=Inf) std([props_red_bw.Circularity]~=Inf) mean([props_green_bw.Circularity]~=Inf) std([props_green_bw.Circularity]~=Inf)...
        mean([props_red_bw.FilledArea]) std([props_red_bw.FilledArea]) mean([props_green_bw.FilledArea]) std([props_green_bw.FilledArea])...
        mean([props_red_bw.EulerNumber]) std([props_red_bw.EulerNumber]) mean([props_green_bw.EulerNumber]) std([props_green_bw.EulerNumber])...
        sum([props_red_bw.EulerNumber]) sum([props_green_bw.EulerNumber])...
        mean([props_red_bw.EquivDiameter]) std([props_red_bw.EquivDiameter]) mean([props_green_bw.EquivDiameter]) std([props_green_bw.EquivDiameter])...
        mean([props_red_bw.Solidity]) std([props_red_bw.Solidity]) mean([props_green_bw.Solidity]) std([props_green_bw.Solidity])...
        mean([props_red_bw.Extent]) std([props_red_bw.Extent]) mean([props_green_bw.Extent]) std([props_green_bw.Extent])...
        mean([props_red_bw.Perimeter]) std([props_red_bw.Perimeter]) mean([props_green_bw.Perimeter]) std([props_green_bw.Perimeter])...
        mean([props_red_bw.PerimeterOld]) std([props_red_bw.PerimeterOld]) mean([props_green_bw.PerimeterOld]) std([props_green_bw.PerimeterOld])...
        mean([props_red_bw.MaxFeretDiameter]) std([props_red_bw.MaxFeretDiameter]) mean([props_green_bw.MaxFeretDiameter]) std([props_green_bw.MaxFeretDiameter])...
        mean([props_red_bw.MaxFeretAngle]) std([props_red_bw.MaxFeretAngle]) mean([props_green_bw.MaxFeretAngle]) std([props_green_bw.MaxFeretAngle])...
        mean([props_red_bw.MinFeretDiameter]) std([props_red_bw.MinFeretDiameter]) mean([props_green_bw.MinFeretDiameter]) std([props_green_bw.MinFeretDiameter])...
        mean([props_red_bw.MinFeretAngle]) std([props_red_bw.MinFeretAngle]) mean([props_green_bw.MinFeretAngle]) std([props_green_bw.MinFeretAngle])...
        circ_red_bw circ_green_bw...
        mean([props_red_skel.Area]) std([props_red_skel.Area]) mean([props_green_skel.Area]) std([props_green_skel.Area])...
        mean([props_red_skel.MajorAxisLength]) std([props_red_skel.MajorAxisLength]) mean([props_green_skel.MajorAxisLength]) std([props_green_skel.MajorAxisLength])...
        mean([props_red_skel.MinorAxisLength]) std([props_red_skel.MinorAxisLength]) mean([props_green_skel.MinorAxisLength]) std([props_green_skel.MinorAxisLength])...
        mean([props_red_skel.Eccentricity]) std([props_red_skel.Eccentricity]) mean([props_green_skel.Eccentricity]) std([props_green_skel.Eccentricity])...
        mean([props_red_skel.ConvexArea]) std([props_red_skel.ConvexArea]) mean([props_green_skel.ConvexArea]) std([props_green_skel.ConvexArea])...
        mean([props_red_skel.Circularity]~=Inf) std([props_red_skel.Circularity]~=Inf) mean([props_green_skel.Circularity]~=Inf) std([props_green_skel.Circularity]~=Inf)...
        mean([props_red_skel.FilledArea]) std([props_red_skel.FilledArea]) mean([props_green_skel.FilledArea]) std([props_green_skel.FilledArea])...
        mean([props_red_skel.EulerNumber]) std([props_red_skel.EulerNumber]) mean([props_green_skel.EulerNumber]) std([props_green_skel.EulerNumber])...
        sum([props_red_skel.EulerNumber]) sum([props_green_skel.EulerNumber])...
        mean([props_red_skel.EquivDiameter]) std([props_red_skel.EquivDiameter]) mean([props_green_skel.EquivDiameter]) std([props_green_skel.EquivDiameter])...
        mean([props_red_skel.Solidity]) std([props_red_skel.Solidity]) mean([props_green_skel.Solidity]) std([props_green_skel.Solidity])...
        mean([props_red_skel.Extent]) std([props_red_skel.Extent]) mean([props_green_skel.Extent]) std([props_green_skel.Extent])...
        mean([props_red_skel.Perimeter]) std([props_red_skel.Perimeter]) mean([props_green_skel.Perimeter]) std([props_green_skel.Perimeter])...
        mean([props_red_skel.PerimeterOld]) std([props_red_skel.PerimeterOld]) mean([props_green_skel.PerimeterOld]) std([props_green_skel.PerimeterOld])...
        mean([props_red_skel.MaxFeretDiameter]) std([props_red_skel.MaxFeretDiameter]) mean([props_green_skel.MaxFeretDiameter]) std([props_green_skel.MaxFeretDiameter])...
        mean([props_red_skel.MaxFeretAngle]) std([props_red_skel.MaxFeretAngle]) mean([props_green_skel.MaxFeretAngle]) std([props_green_skel.MaxFeretAngle])...
        mean([props_red_skel.MinFeretDiameter]) std([props_red_skel.MinFeretDiameter]) mean([props_green_skel.MinFeretDiameter]) std([props_green_skel.MinFeretDiameter])...
        mean([props_red_skel.MinFeretAngle]) std([props_red_skel.MinFeretAngle]) mean([props_green_skel.MinFeretAngle]) std([props_green_skel.MinFeretAngle])...
        circ_red_skel circ_green_skel...
        num_red num_red_branch num_green num_green_branch];

    counter = counter + 1;
    progress = 100*counter/(nfiles)
    
end

quantified = cell2mat(quantified');

quantified_props = quantified(1:end,1:126);
quantified_props_bw = quantified(1:end,127:224);
quantified_props_skel = quantified(1:end,225:end);

%% Class labels - Mascharak, desJardins-Park, Januszyk et al

for i = 1:length(stored_names)
    if contains(stored_names{i},'- P') & contains(stored_names{i},'POD2')
        class_pt{i} = 'PBS POD2';
        color{i} = [1 0 0];
    end
    if contains(stored_names{i},'- P') & contains(stored_names{i},'POD7')
        class_pt{i} = 'PBS POD7';
        color{i} = [0 1 0];
    end
    if contains(stored_names{i},'- P') & contains(stored_names{i},'POD14')
        class_pt{i} = 'PBS POD14';
        color{i} = [0 0 1];
    end
    if contains(stored_names{i},'- P') & contains(stored_names{i},'POD30')
        class_pt{i} = 'PBS POD30';
        color{i} = [1 1 0];
    end
    if contains(stored_names{i},'- U') | contains(stored_names{i},'Unwounded')
        class_pt{i} = 'Unwounded';
        color{i} = [0 1 1];
    end
    if contains(stored_names{i},'- V') & contains(stored_names{i},'POD2')
        class_pt{i} = 'Vert POD2';
        color{i} = [1 0 1];
    end
    if contains(stored_names{i},'- V') & contains(stored_names{i},'POD7')
        class_pt{i} = 'Vert POD7';
        color{i} = [1 1 1];
    end
    if contains(stored_names{i},'- V') & contains(stored_names{i},'POD14')
        class_pt{i} = 'Vert POD14';
        color{i} = [0.5 0.5 1];
    end
    if contains(stored_names{i},'- V') & contains(stored_names{i},'POD30')
        class_pt{i} = 'Vert POD30';
        color{i} = [0.5 0.5 0.5];
    end
end

%% t-SNE Visualization

Y = tsne(quantified,'Algorithm','exact','Distance','euclidian','Exaggeration',4);
figure
gscatter(Y(:,1),Y(:,2),class_pt')
