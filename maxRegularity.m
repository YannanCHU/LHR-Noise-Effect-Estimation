function [M_R, M_R_angle] = maxRegularity(im)
% inputs:
% im: a portion of image
% outputs:
% M_R: the measured maximal regularity
% M_R_angle: the corresponding angle (resolution is 1 degree)
[M,N] = size(im);

% normalised autocorrelation function
Pxy = ifftshift(ifft2(fft2(im) .* conj(fft2(im))));

% polar coordinate
radiusResolution = round(min([M;N])/2+1);
angleResolution = 360;
Ppol = ImToPolar(Pxy, 0, 1, radiusResolution, angleResolution).';
Ppol = Ppol / max(Ppol(:));

% autocorrelation based interaction map
Mpol = 1 - Ppol;

% A row of Mpol(i, j) is called a contrast function. This name reflects the
% close relation between the autocorrelation and the mean square grey-level
% difference (contrast).

% periodic structure will give a contrast function with deep and periodic minima.
% Our definition of regularity quantifies this property.

Rdir = zeros(1,angleResolution); 
for angle_i = 1:1:angleResolution
    %% Procedure 1. Finding the extrema of Fi(d)
    Fi_d0 = Mpol(angle_i,:);
    Fi_d = medfilt1(Fi_d0);

    [Fi_d_max,max_idx,Fi_d_min,min_idx] = extrema(Fi_d);
    if isempty(min_idx == 1) == false
        Fi_d_min(min_idx == 1) = [];
        min_idx(min_idx == 1) = [];
    end
    if isempty(min_idx == radiusResolution) == false
        Fi_d_min(min_idx == radiusResolution) = [];
        min_idx(min_idx == radiusResolution) = [];
    end

    if isempty(max_idx == 1) == false
        Fi_d_max(max_idx == 1) = [];
        max_idx(max_idx == 1) = [];
    end
    if isempty(max_idx == radiusResolution) == false
        Fi_d_max(max_idx == radiusResolution) = [];
        max_idx(max_idx == radiusResolution) = [];
    end

    Nmin = length(Fi_d_min);
    Nmax = length(Fi_d_max);

    %% Procedure 2. Computing intensity regularity Rint
    [max_idx_sorted, sortIdx_max] = sort(max_idx,'ascend');
    % [min_idx_sorted, sortIdx_min] = sort(min_idx,'ascend');

    Fi_d_max_sorted = Fi_d_max(sortIdx_max);
    % Fi_d_min_sorted = Fi_d_min(sortIdx_min);

    Fmin_Fmax_pair = zeros(Nmax,2);
    for left = 1:1:Nmax
        [~, right] = find(Fi_d_max_sorted(left+1:end)>Fi_d_max_sorted(left), 1, 'first');
        right = left + right;

        intervalBegin = max_idx_sorted(left);
        intervalEnd = max_idx_sorted(right);

        if isempty(intervalEnd)
            intervalEnd = radiusResolution;
        end

        min_index_within_interval = min_idx(min_idx > intervalBegin & min_idx < intervalEnd);

        if isempty(min_index_within_interval) == false
            Fmin_temp = min(Fi_d(min_index_within_interval));
            Fmax_temp = Fi_d(intervalBegin);
            %         disp("Fmax'(d) and Fmin'(d): "+Fmax+", " + Fmin);
            Fmin_Fmax_pair(left,:) = [Fmin_temp, Fmax_temp];
        end
    end

    [~, largestAmplitudeIdx] = max(Fmin_Fmax_pair(:,2)-Fmin_Fmax_pair(:,1));
    Fmin = Fmin_Fmax_pair(largestAmplitudeIdx,1);
    Fmax = Fmin_Fmax_pair(largestAmplitudeIdx,2);
    % skip the rectification to reduce computational complexity

    % intensity regularity
    Rint = 1 - Fmin / Fmax;

    if Nmin == 1    % only one minima
        % if all maximum amplitudes are equal, find the first maxima's d
        [~, dmax] = max(Fi_d_max_sorted);
        dmax = max_idx_sorted(dmax)-1;          % -1 is used because the d value should start from 0 instead of 1
        [~, d1] = min(Fi_d_min);
        d1 = min_idx(d1)-1;
        Rpos = 1 - abs(d1 - 2*dmax) / d1;
    elseif Nmin == 0 || Nmax == 0
        Rpos = 0;
        Rint = 0;
    else            % more than one minima
        [F_two_min, index_2_min] = mink(Fi_d_min,2);
        d1_and_d2 = min_idx(index_2_min);
        d1 = min(d1_and_d2) - 1;
        d2 = max(d1_and_d2) - 1;
        F1 = Fi_d(d1+1);
        F2 = Fi_d(d2+1);

        % check whether there are minima within the interval [d1, d2]
        index_temp = find(min_idx > (d1+1) & min_idx < (d2+1));
        if isempty(index_temp)      % no minima between d1 and d2 - normal
            gamma = d1 / d2;  % range from 0 to 1
            % Rpos = 1 - abs(1 - 2*gamma);
            Rpos = 2 * gamma * (gamma > 0 && gamma <= 0.5) + ...
                (2 - 2 * gamma) * (gamma > 0.5 && gamma < 1);
        else            % exist minima between d1 and d2 - special
            gamma = d1 / d2;  % range from 0 to 1
            % Rpos = 1 - abs(1 - 2*gamma);
            % Rpos2 = max([1 - abs(1 - 3*gamma), 0]);
            % Rpos = max(Rpos, Rpos2);
            Rpos = 2 * gamma * (gamma > 0 && gamma <= 0.5) + ...
                (2 - 2 * gamma) * (gamma > 0.5 && gamma < 1);
            Rpos2 = (3*gamma) * (gamma > 0 && gamma <= 1/3) + ...
                (2-3*gamma) * (gamma > 1/3 && gamma <= 2/3) + ...
                0 * (gamma > 2/3 && gamma < 1);
            Rpos = max(Rpos, Rpos2);
        end
    end

    % Rint (intensity) and Rpos (position) are invariant to scaling of d
    % For an angle i, the directional regularity is defined as
    try
        Ri = (Rint*Rpos)^2;
    catch
        disp("Error, the angle is " + angle_i);
    end
    Rdir(angle_i) = Ri;
end

[Tk, Tk_idx, ~, ~] = extrema(Rdir);
Rthr = 0.025;        % the R value smaller than 0.25 cannot be sensed
                    % small R is used for the potential weak structure
smallR_idx = find(Tk < Rthr);
Tk_idx(smallR_idx) = [];
Tk(smallR_idx) = [];        

% maximal regularity and other statistics
[M_R, M_R_idx] = max(Tk);      % maximal regularity
M_R_angle = Tk_idx(M_R_idx);   % angle corresponding to M_R
% mu_R = mean(Tk);    % mean
% var_R = var(Tk);                        % variance
% v_R = length(Tk) / angleResolution;     % maxima density

% fprintf("The maximal regularity is %2.3f, the corresponding angle is %2.1f°\n", M_R, M_R_angle-1);
end