function [sampling_pattern,ky,kz] = caspr_generate_sampling_pattern(ydim,zdim,turbo_factor,nshots,revolutions,oio,vd,varargin)
%%Generate cartesian spiral sampling patterns
% Summary of the function:
%
% 1) Determine optimal number of revolutions (close to the original)
% The number of revolutions dictates potential overlap of points between
% spirals. To minimize these effects we approximate the sampling
% uniformity at the periphery of k-space and evaluate the uniformity for
% 5000 different number of revolutions. Then we select the number of
% revolution that provides the best score.
%
% 2) Interpolate an analytical spiral to a Cartesian grid. Within this
% interpolation we check for duplicate phase encodes within a spiral. If we
% encounter duplicates we iteratively add/subtract 1 from the py/pz coordinate
% untill we find a unique phase encode. In addition, we add small random shifts 
% on the radius of the point before interpolation. This prevents circular gaps from forming in the k-space.  
%
% 3) The phase encode train can simply be converted to a in-out caspr with simple indexing operations. 
% 
% inputs:   ydim = matrix size in Y
%           zdim = matrix size in Z
%           turbo_factor = echo train length
%           nshots = number of echo trains
%           revolutions = approximate number of spiral rotations
%           oio = 0 or 1, with out or in-out respectively.
%
% output:   sampling_pattern (dimensions [ydim, zdim, nshots])

%% Check input and assign default values
[args,fn] = check_arguments(nargin);
for n = 1 : numel(fn)
    eval([fn{n},' = ',num2str(args.(fn{n})),';']);
end

%% 1) Determine optimal number of revolutions
load('random_numb.mat')                     % List of random numbers
n_rand_samp  = 5000;                        % Number of randomly set revolutions to evaluate
golden_angle = (3. - sqrt(5.)) * pi;        % Golden ratio 
ang_thr      = max([turbo_factor / 30 2]);  % Number of points at the end of spiral to consider
cost = [];
for n = 1 : n_rand_samp
    revols(n) = (revolutions + (rand_numb(n)-.5));
    pe_dk        = revols(n) * 2 * pi / turbo_factor;
    azu_angles = [];
    for nn = 1 : round(ang_thr)
        azu_angles = mod(cat(2,azu_angles,[(0:min([nshots nshots])-1) * golden_angle - (nn-1) * pe_dk ]),2*pi);
    end
    
    % Look for an approximately uniform distribution of the phase encodes at the periphery
    cost(n)    = sum(abs(abs(diff(sort(azu_angles))) - (2*pi/(min([nshots nshots])*round(ang_thr)))).^4);
end
[~,minpos] = min(cost);

revolution = revols(minpos) * 2 * pi;
disp(['Revolutions changed from ',num2str(revolutions),' / ',num2str(revolution / (2*pi))])

%% Sample the analytical spiral
pyoffset = round(ydim/2+.5); 
pzoffset = round(zdim/2+.5);

shifty = pyoffset / turbo_factor; % Radial distance between smaples on a spiral interleave
shiftz = pzoffset / turbo_factor;

sampling_pattern = zeros(ydim,zdim,nshots);
ky               = 2 * ydim * zdim * zeros(turbo_factor,nshots);
kz               = 2 * ydim * zdim * ones(turbo_factor,nshots);
for shot=0:(nshots-1)
    ringRadiusY = 0.0;
    ringRadiusZ = 0.0;
    for point=0:(turbo_factor-1)
        theta = (point * revolution / turbo_factor) + shot * golden_angle;
        ringRadiusY = ydim / 2 * (point / turbo_factor)^vd;
        ringRadiusZ = zdim / 2 * (point / turbo_factor)^vd;
        
        % Add random shift in the radial direction to prevent empty rings in k-space
        if point > 0
            vd_fac = (((point+1) / turbo_factor)^vd - ((point-1) / turbo_factor)^vd) / 2 / (1/turbo_factor);
        else
            vd_fac = 0;
        end
        py = (ringRadiusY + shifty * (rand_numb(shot+1)-0.5) * vd_fac) * cos(theta);
        pz = (ringRadiusZ + shiftz * (rand_numb(shot+1)-0.5) * vd_fac) * sin(theta);
               
        
        pyrounded = round(py);
        if pyrounded == -pyoffset
            pyrounded = pyrounded + 1;
        end
        pzrounded = round(pz);
        if pzrounded == -pzoffset
            pzrounded = pzrounded + 1;
        end
            
        % Check if the phase encode is already in the shot, if so start
        % checking neighbors
        while ~isempty(find(ky(:,1+shot) + 1j * kz(:,1+shot) == pyrounded + pyoffset + 1j * (pzrounded + pzoffset)))
            if rand() > 0.5
                pyrounded = pyrounded + (-1)^randi(2);
            else
                pzrounded = pzrounded + (-1)^randi(2);
            end
        end
        
        % Check for boundary conditions
        while pyrounded+pyoffset > ydim
            pyrounded = pyrounded - 1;
        end
        while pzrounded+pzoffset > zdim
            pzrounded = pzrounded - 1;
        end
        
        % Assign to matrix
        ky(1+point,1+shot) = pyrounded + pyoffset;
        kz(1+point,1+shot) = pzrounded + pzoffset;
        sampling_pattern(pyrounded+pyoffset,pzrounded+pzoffset,shot+1)=1;  
    end
end

%% Transform to out-in-out
if oio == 1
    if mod(size(ky,1),2) > 0 % uneven
        ky = ky([end:-2:1 2:2:end],:);
        kz = kz([end:-2:1 2:2:end],:);
    else
        ky = ky([end:-2:1 1:2:end],:);
        kz = kz([end:-2:1 1:2:end],:);
    end
end

%ky = ky(:);
%kz = kz(:);

% Display acceleration
m = sum(sampling_pattern,3);
disp(['Acceleration R = ',num2str(numel(m) / nnz(m))])
end


function [res,fn] = check_arguments(nargin)
    if nargin > 6 
        res = {}
        fn = '';
        return
    end
    if nargin < 7
        res.vd = 1;
    end
    if nargin < 6
        res.oio = 0;
    end
    if nargin < 5
        res.revolutions = 1;
    end
    if nargin < 4
        res.nshots = 100;
    end
    if nargin < 3
        res.turbo_factor = 50;
    end
    if nargin < 2
        res.zdim = 100;
    end
    if nargin < 1
        res.ydim = 100;
    end
    
    fn = fieldnames(res);
end