%% Fiber Quantification Algorithm v3
% Created by Shamik Mascharak
% November 28th, 2019 - Created
% December 1st, 2019 - First iteration completed
% October 28th, 2020 - Uploaded to Github
% September 22nd, 2022 - Code updated

%% Read files and generate binarized images

clear
clc
close all

% Get list of all .tif files in the directory
imagefiles = dir('*.tif');
% imagefiles = dir('*.png');
nfiles = length(imagefiles);

%% Run color deconvolution

counter = 0;

% Optional scaling factors (to convert pixels to microns)
scale = 1;

% Define target image for image normalization
target_im = imread(imagefiles(1544).name); % Input index from imagefiles variable

verbose = 0; % 1 to display results, 0 to not

% Define stain vector for color deconvolution
trichrome_vect = [0.7995107 0.5913521 0.10528667;...
    0.099971585 0.73738605 0.6680326;...
    0.59227383 0.3264422 0.7366459];

% Define warning counters for images that are blank/empty after processing
warn_count_blue = 0;
warn_count_blue_bw = 0;
warn_count_blue_bw_branch = 0;
warn_count_blue_bw_skel = 0;
warn_count_blue_bw_strel = 0;
warn_count_blue_noise = 0;

for i = 1:nfiles
    filename = imagefiles(i).name;
    stored_names{i} = filename;
    stored_im{i} = imread(filename); % Store the original image
    
    % Normalize image using RGB histogram method
   stored_norm_im{i} = Norm(stored_im{i},target_im,'RGBHist',verbose);
%   stored_norm_im{i} = stored_im{i}; % To bypass normalization

    % Perform color deconvolution (Ruifrok method)
    [Dch M] = Deconvolve(stored_norm_im{i},trichrome_vect,0);
    
    [blue red green] = PseudoColourStains(Dch,M);
    
    trichrome_blue = (255-rgb2gray(blue))-(255-rgb2gray(green));
    
    % Remove noise with adaptive filtering
    trichrome_blue_noise = wiener2(trichrome_blue,[3 3]);
    
    stored_trichrome_blue_noise{i} = trichrome_blue_noise;
    
    trichrome_blue_bw = im2bw(imadjust(trichrome_blue_noise));
    
    % Define diamond structuring element
    seD = strel('diamond',1);
    trichrome_blue_bw_strel = imerode(trichrome_blue_bw,seD);
    
    stored_trichrome_blue_strel{i} = trichrome_blue_bw_strel; % Binary map
    
    % Skeletonize, separate at branchpoints, and delete single pixel lines
    trichrome_blue_bw_skel = bwmorph(trichrome_blue_bw_strel,'skel',Inf);
    trichrome_blue_bw_branch = bwmorph(trichrome_blue_bw_skel,'branchpoints');
    trichrome_blue_bw_skel(trichrome_blue_bw_branch == 1) = 0;
    trichrome_blue_bw_skel = bwmorph(trichrome_blue_bw_skel,'clean');
    
    stored_trichrome_blue_skel{i} = trichrome_blue_bw_skel; % Skeletonized image
    stored_trichrome_blue_branch{i} = trichrome_blue_bw_branch;
    
    % Store indices of images that are blank/empty after processing. These
    % indices can be found in the variables containing the 'warnings_' prefix.
    if not(any(any(trichrome_blue)))
        warnings_trichrome_blue(warn_count_blue+1) = i;
        warn_count_blue = warn_count_blue + 1;
    end
    
    if not(any(any(trichrome_blue_bw)))
        warnings_trichrome_blue_bw(warn_count_blue_bw+1) = i;
        warn_count_blue_bw = warn_count_blue_bw + 1;
    end

    if not(any(any(trichrome_blue_bw_branch)))
        warnings_trichrome_blue_bw_branch(warn_count_blue_bw_branch+1) = i;
        warn_count_blue_bw_branch = warn_count_blue_bw_branch + 1;
    end
    
    if not(any(any(trichrome_blue_bw_skel)))
        warnings_trichrome_blue_bw_skel(warn_count_blue_bw_skel+1) = i;
        warn_count_blue_bw_skel = warn_count_blue_bw_skel + 1;
    end   
    
    if not(any(any(trichrome_blue_bw_strel)))
        warnings_trichrome_blue_bw_strel(warn_count_blue_bw_strel+1) = i;
        warn_count_blue_bw_strel = warn_count_blue_bw_strel + 1;
    end    
    
    if not(any(any(trichrome_blue_noise)))
        warnings_trichrome_blue_noise(warn_count_blue_noise+1) = i;
        warn_count_blue_noise = warn_count_blue_noise + 1;
    end
    
    % Display progress
    counter = counter + 1;
    progress = 100*counter/(nfiles)
    
end

stored_names = stored_names';

%% Quantify image features

clc
counter = 0;

for i = 1:nfiles
    
    blue = stored_trichrome_blue_noise{i};
    blue_bw = stored_trichrome_blue_strel{i};
    blue_skel = stored_trichrome_blue_skel{i};
    blue_branch = stored_trichrome_blue_branch{i};
    
    % Calculate values from grayscale images
    props_blue = regionprops(blue_bw,blue,'all');
    % Fiber angle-related values
    alpha_blue = [props_blue.Orientation]*pi/180;
    circ_blue = [mean(alpha_blue) median(alpha_blue) std(alpha_blue) skewness(alpha_blue)...
        kurtosis(alpha_blue) circ_kappa(alpha_blue)];

    % Calculate Haralick features from grayscale images (4 different offsets)
    offset = [0 1; -1 1;-1 0;-1 -1];
    glcm_blue = graycomatrix(blue,'Offset',offset,'Symmetric',true);
    
    grayprops_blue = graycoprops(glcm_blue);
    
    % Calculate values from binary images
    props_blue_bw = regionprops(blue_bw,'all');
    % Fiber angle-related values
    alpha_blue_bw = [props_blue_bw.Orientation]*pi/180;
    circ_blue_bw = [mean(alpha_blue_bw) median(alpha_blue_bw) std(alpha_blue_bw) skewness(alpha_blue_bw)...
        kurtosis(alpha_blue_bw) circ_kappa(alpha_blue_bw)];

    % Calculate values from skeletonized images
    props_blue_skel = regionprops(blue_skel,'all');
    % Fiber angle-related values
    alpha_blue_skel = [props_blue_skel.Orientation]*pi/180;
    circ_blue_skel = [mean(alpha_blue_skel) median(alpha_blue_skel) std(alpha_blue_skel) skewness(alpha_blue_skel)...
        kurtosis(alpha_blue_skel) circ_kappa(alpha_blue_skel)];
    
    % Branchpoint-related values
    [L_blue, num_blue] = bwlabel(blue_skel);
    [~, num_blue_branch] = bwlabel(blue_branch);
    
    % Grayscale regionprops, grayscale circ, Haralick features, binary regionprops,
    % binary circ, skeletonized regionprops, skeletonized circ, branchpoints,
    
    quantified{i} = [mean([props_blue.Area]) std([props_blue.Area])...
        mean([props_blue.MajorAxisLength]) std([props_blue.MajorAxisLength])...
        mean([props_blue.MinorAxisLength]) std([props_blue.MinorAxisLength])...
        mean([props_blue.Eccentricity]) std([props_blue.Eccentricity])...
        mean([props_blue.ConvexArea]) std([props_blue.ConvexArea])...
        mean([props_blue.Circularity]~=Inf) std([props_blue.Circularity]~=Inf)...
        mean([props_blue.FilledArea]) std([props_blue.FilledArea])...
        mean([props_blue.EulerNumber]) std([props_blue.EulerNumber])...
        sum([props_blue.EulerNumber])...
        mean([props_blue.EquivDiameter]) std([props_blue.EquivDiameter])...
        mean([props_blue.Solidity]) std([props_blue.Solidity])...
        mean([props_blue.Extent]) std([props_blue.Extent])...
        mean([props_blue.Perimeter]) std([props_blue.Perimeter])...
        mean([props_blue.PerimeterOld]) std([props_blue.PerimeterOld])...
        mean([props_blue.MeanIntensity]) std([props_blue.MeanIntensity])...
        mean([props_blue.MinIntensity]) std(double([props_blue.MinIntensity]))...
        mean([props_blue.MaxIntensity]) std(double([props_blue.MaxIntensity]))...
        mean([props_blue.MaxFeretDiameter]) std([props_blue.MaxFeretDiameter])...
        mean([props_blue.MaxFeretAngle]) std([props_blue.MaxFeretAngle])...
        mean([props_blue.MinFeretDiameter]) std([props_blue.MinFeretDiameter])...
        mean([props_blue.MinFeretAngle]) std([props_blue.MinFeretAngle])...
        circ_blue...
        [grayprops_blue.Contrast] [grayprops_blue.Correlation] [grayprops_blue.Energy] [grayprops_blue.Homogeneity]...
        mean([props_blue_bw.Area]) std([props_blue_bw.Area])...
        mean([props_blue_bw.MajorAxisLength]) std([props_blue_bw.MajorAxisLength])...
        mean([props_blue_bw.MinorAxisLength]) std([props_blue_bw.MinorAxisLength])...
        mean([props_blue_bw.Eccentricity]) std([props_blue_bw.Eccentricity])...
        mean([props_blue_bw.ConvexArea]) std([props_blue_bw.ConvexArea])...
        mean([props_blue_bw.Circularity]~=Inf) std([props_blue_bw.Circularity]~=Inf)...
        mean([props_blue_bw.FilledArea]) std([props_blue_bw.FilledArea])...
        mean([props_blue_bw.EulerNumber]) std([props_blue_bw.EulerNumber])...
        sum([props_blue_bw.EulerNumber])...
        mean([props_blue_bw.EquivDiameter]) std([props_blue_bw.EquivDiameter])...
        mean([props_blue_bw.Solidity]) std([props_blue_bw.Solidity])...
        mean([props_blue_bw.Extent]) std([props_blue_bw.Extent])...
        mean([props_blue_bw.Perimeter]) std([props_blue_bw.Perimeter])...
        mean([props_blue_bw.PerimeterOld]) std([props_blue_bw.PerimeterOld])...
        mean([props_blue_bw.MaxFeretDiameter]) std([props_blue_bw.MaxFeretDiameter])...
        mean([props_blue_bw.MaxFeretAngle]) std([props_blue_bw.MaxFeretAngle])...
        mean([props_blue_bw.MinFeretDiameter]) std([props_blue_bw.MinFeretDiameter])...
        mean([props_blue_bw.MinFeretAngle]) std([props_blue_bw.MinFeretAngle])...
        circ_blue_bw...
        mean([props_blue_skel.Area]) std([props_blue_skel.Area])...
        mean([props_blue_skel.MajorAxisLength]) std([props_blue_skel.MajorAxisLength])...
        mean([props_blue_skel.MinorAxisLength]) std([props_blue_skel.MinorAxisLength])...
        mean([props_blue_skel.Eccentricity]) std([props_blue_skel.Eccentricity])...
        mean([props_blue_skel.ConvexArea]) std([props_blue_skel.ConvexArea])...
        mean([props_blue_skel.Circularity]~=Inf) std([props_blue_skel.Circularity]~=Inf)...
        mean([props_blue_skel.FilledArea]) std([props_blue_skel.FilledArea])...
        mean([props_blue_skel.EulerNumber]) std([props_blue_skel.EulerNumber])...
        sum([props_blue_skel.EulerNumber])...
        mean([props_blue_skel.EquivDiameter]) std([props_blue_skel.EquivDiameter])...
        mean([props_blue_skel.Solidity]) std([props_blue_skel.Solidity])...
        mean([props_blue_skel.Extent]) std([props_blue_skel.Extent])...
        mean([props_blue_skel.Perimeter]) std([props_blue_skel.Perimeter])...
        mean([props_blue_skel.PerimeterOld]) std([props_blue_skel.PerimeterOld])...
        mean([props_blue_skel.MaxFeretDiameter]) std([props_blue_skel.MaxFeretDiameter])...
        mean([props_blue_skel.MaxFeretAngle]) std([props_blue_skel.MaxFeretAngle])...
        mean([props_blue_skel.MinFeretDiameter]) std([props_blue_skel.MinFeretDiameter])...
        mean([props_blue_skel.MinFeretAngle]) std([props_blue_skel.MinFeretAngle])...
        circ_blue_skel...
        num_blue num_blue_branch];

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