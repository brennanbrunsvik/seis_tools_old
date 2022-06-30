function [zlayt,zlayb,vlay,varargout] = layerise(Z,V,dvmin,ifplot,varargin)


% Code to take a continuous series of velocity (V) with depth (Z) and
% convert it to several layers of constant velocity. 
% New approach (brb2022.06.29) is to traverse down the model and make a new
% layer when either the velocity or depth (or norm of those two) has
% changed suffiently. 
% 
% INPUTS:
%   Z     - [N x 1] vector of depths (km)
%   V     - [N x 1] vector of velocities (km/s)
%   dvmin - brb2022.06.29 removed temporarily. 
%           Tunable variable specifying layer thickness. If this value is
%           smaller, the continuous profile will be better fit, but with
%           more layers. If the value is larger, the fit will be poorer,
%           but there will be fewer layers. Recommended value: 0.1 (km/s)
%   ifplot - option to plot (true) or not (false)
%   varargin - any number of [N x 1] vectors for other parameters (e.g. Vp, 
%           rho) to be discretised onto the same mesh. The discontinuity 
%           depths will be well-preserved, but because the algorithm isn't
%           actually adapting to these parameters, the jumps across the
%           discontinuities will be less precisely recovered. 
% 
% OUTPUTS:
%   zlayt - [Nlay x 1] vector of depths of layer tops (km)
%   zlayb - [Nlay x 1] vector of depths of layer bases (km)
%   vlay  - [Nlay x 1] vector of velocity values for each layer (km/s)
%   varargout - any number of [Nlay x 1] vector of other parameter values
%           in each layer

%     Z. Eilon 08/2016
%     Modified heavily brb2022.06.29. 

break_v    = 0.005 ; % Tend toward break into a layer every this much velocity change. 
break_z    = 30   ; % Tend toward making all layers smaller than this. 
break_grad = 0.025; % If gradient in v between two provided depths is greater than this, divide a layer here. Mostly this value is irrelevant, we just identify grad==inf, i.e. discontinuities. 

if nargin<4 || isempty(ifplot)
    ifplot=false;
end

dz = diff(Z);
dv = diff(V);
dvdz = dv./dz; % Gradient. 
is_disc = find(isinf(dvdz)); % Top index of discontinuities. dz = 0. 

% Loop to find the depths where we will break up layers. 
% Put in a layer where there is strong gradients. 
% Or if we haven't put in a layer for a long distance. 
% Or a combination of the two. 
% Manually assign first and last Z as layer breaks. 
break_inds = [1]; % Indecies where we will break into layers. 
for iz = 2:(length(Z)-1);

    if dvdz(iz) > break_grad; % Make sure discontinuities (e.g. grad == inf) have a layer break. 
        break_inds = [break_inds; iz]; 
        continue; 
    end 
    
    vchange = abs(V(iz) - V(break_inds(end))); 
    zchange = Z(iz)     - Z(break_inds(end) ); 
    
    tendchange = sqrt( (vchange / break_v)^2 ...
                      +(zchange / break_z)^2); % Will be 1 if vchange or zchange meets our predefined values of when to break layers, but also will be one if a combination of v and z change is high enough. 
    
    if tendchange > 1; 
        break_inds = [break_inds; iz]; 
    end
end
break_inds = unique([break_inds; length(Z)]); % Manually assign last Z as a layer break. 

zlayt = Z(break_inds(1:end-1)); % Top depth of each layer. 
zlayb = Z(break_inds(2:end  )); % Bottom depth of each layer. 

dup_vals = (zlayt == zlayb); % This can catch duplicate depths in one array - Remove accidentally added infinitesimally small layers. 
zlayt = zlayt(~dup_vals); 
zlayb = zlayb(~dup_vals); 

% Look for some obvious signs of mistakes.  
layer_gaps = zlayt(2:end) ~= zlayb(1:end-1); 
if any(layer_gaps); 
    error('The bottom of some layer is not the same as the top of the next layer'); 
elseif Z(end) ~= zlayb(end);
    error('Bottom of new layers not equal to bottom of old layers'); 
elseif Z(1) ~= zlayt(1); 
    error('Top of new layers not equal to top of old layers'); 
end

%%% Get ready to take average of velocity and whatever within each layer. 
nlay               = length(zlayt); % Number of layers in new model
vlay               = nan(size(zlayt)); % Output velocity

% Keep track of whether above or below layer boundary. For discontinuities. 
Z_shift            = Z; 
Z_shift(is_disc  ) = Z_shift(is_disc  ) - 0.00001; 
Z_shift(is_disc+1) = Z_shift(is_disc+1) + 0.00001; 

% Prepare varagin/out for processing.
process_varargin = ~isempty(varargin);  
if process_varargin; % Prepare to process varargin, if it was provided. If needed for testing: varargin = [{V+1} {V+2}]; 
    nargs = length(varargin); 
    varargin_arr = nan(size(Z,1), nargs); % Make into an array for vectorized processing    
    for iv = 1:nargs; 
        varargin_arr(:,iv) = varargin{iv}; 
    end
    varargout_arr = nan(nlay, nargs); 
else
    varargout = {}; 
end

% Important loop - get average of v and varargin within each layer. If needed for testing: disp([Z(in_lay)' ; V(in_lay)']); % See what depths and velocities you are taking. 
for iz = 1:length(zlayt); 
    in_lay = and( ( Z_shift >= zlayt(iz) ) , ( Z_shift <= zlayb(iz) ) ); % Within this layer or not? Note: This realies on Z_shift having discontinuity depths modified. 
    vlay(iz) = mean(V(in_lay)); % Mean of velocity from top to bottom of this layer. 
    if process_varargin; 
        varargout_arr(iz, :) = mean(varargin_arr(in_lay,:),1); % Mean of all varargin varaibles, from top to bottom of this layer. Take mean along first dimension.  
    end
end

if any(isnan(vlay)); 
    error('Did not get velocity for some layer. '); 
end

% Convert varargout_arr to an array of cells. 
if process_varargin; 
    if any(isnan(varargout_arr)); 
        error('A varargout value is nan. Might not have been averaged or provided correctly?'); 
    end
    
    varargout = {}; 
    for iarg = 1:nargs; 
        varargout(:,end+1) = {varargout_arr(:,iarg)}; 
    end
end
         
ifplot = false; if ifplot; warning('Setting ifplot = true'); end; 
if ifplot
    figure(11); clf; set(gcf,'pos', [-581 247 514 796], 'color', 'white'); 
    h = tiledlayout(1, 1+length(varargin),'TileSpacing','compact'); 
    nexttile; hold on; box on; 

    plot(V,Z,'-ko')
    zlayp = reshape([zlayt';zlayb'],2*nlay,1);
    vlayp = reshape([vlay';vlay'],2*nlay,1);
    plot(vlayp,zlayp,'-ro')
    set(gca,'ydir','reverse',...
        'ylim',[0, max(Z)],'xlim',[0.9*min(V) 1.1*max(V)])
    
    % other vars
    for iv = 1:length(varargin)
        nexttile; hold on; box on; 
        plot(varargin{iv},Z,'-ko')
        ivlayp = reshape([varargout{iv}';varargout{iv}'],2*nlay,1);
        plot(ivlayp,zlayp,'-ro')
        set(gca,'ydir','reverse','ylim',[0, max(Z)],...
            'xlim',[0.9*min(varargin{iv}) 1.1*max(varargin{iv})])    
    end
    text(min(V), 250, sprintf('nlay=%1.0f\nbreak z=%1.1f\nbreak v=%1.3f',...
        nlay,break_z,break_v), 'fontsize', 12); 
end

end % End function

% % % % % [zlayt,zlayb,vlay,varargout] = layerise(Z,V,dvmin,ifplot,varargin)
% % % % % 
% % % % % Code to take a continuous series of velocity (V) with depth (Z) and
% % % % % convert it to several layers of constant velocity. The approach is to
% % % % % find discontinuities and constant velocity regions and then fill in any
% % % % % gaps around these with stair-step gradients with constant dV in each step
% % % % % 
% % % % % INPUTS:
% % % % %   Z     - [N x 1] vector of depths (km)
% % % % %   V     - [N x 1] vector of velocities (km/s)
% % % % %   dvmin - tunable variable specifying layer thickness. If this value is
% % % % %           smaller, the continuous profile will be better fit, but with
% % % % %           more layers. If the value is larger, the fit will be poorer,
% % % % %           but there will be fewer layers. Recommended value: 0.1 (km/s)
% % % % %   ifplot - option to plot (true) or not (false)
% % % % %   varargin - any number of [N x 1] vectors for other parameters (e.g. Vp, 
% % % % %           rho) to be discretised onto the same mesh. The discontinuity 
% % % % %           depths will be well-preserved, but because the algorithm isn't
% % % % %           actually adapting to these parameters, the jumps across the
% % % % %           discontinuities will be less precisely recovered. 
% % % % % 
% % % % % OUTPUTS:
% % % % %   zlayt - [Nlay x 1] vector of depths of layer tops (km)
% % % % %   zlayb - [Nlay x 1] vector of depths of layer bases (km)
% % % % %   vlay  - [Nlay x 1] vector of velocity values for each layer (km/s)
% % % % %   varargout - any number of [Nlay x 1] vector of other parameter values
% % % % %           in each layer
% % % % 
% % % % %     Z. Eilon 08/2016
% % % % 
% % % % if nargin<4 || isempty(ifplot)
% % % %     ifplot=false;
% % % % end
% % % % 
% % % % % you can fiddle with these if you want to alter how it identifies
% % % % % discontinuities or constant V regions, but these values should work fine
% % % % % for you
% % % % % discgrad = 0.15; % gradients of more than 0.4 m/s per m are deemed discontinuities
% % % % discgrad = inf; % gradients of more than 0.1 m/s per m are deemed discontinuities
% % % % constgrad = -inf; % gradients of less than 2 m/s per 1000m are deemed constant V
% % % % % discgrad = .1; % gradients of more than 0.1 m/s per m are deemed discontinuities
% % % % % constgrad = .0001; % gradients of less than 2 m/s per 1000m are deemed constant V
% % % % %%%
% % % % % discgrad = 0.1; % gradients of more than 0.1 m/s per m are deemed discontinuities
% % % % % constgrad = 0.001; % gradients of less than 2 m/s per 1000m are deemed constant V
% % % % % brb2022.05.23 if constgrad is large (~0.0001) then propmat often breaks,
% % % % % returns nonsense traces, nans, or similar. Set to extremely small value
% % % % % to make the constgrad thing not take effect. 
% % % % 
% % % % 
% % % % Zdo = Z(:);Vdo = V(:);
% % % % 
% % % % dz = diff(Zdo);
% % % % dv = diff(Vdo);
% % % % grad = dv./dz;
% % % % 
% % % % %% ---------------  DISCONTINUITIES  ---------------  %%
% % % % dind = find(abs(grad) >= discgrad); % find all discontinuities (inc gradients)
% % % % if isempty(dind),dind = find(abs(grad) >= 0.075);end % find all discontinuities (inc gradients)
% % % % if isempty(dind),dind = find(abs(grad) == max(abs(grad)));end % find all discontinuities (inc gradients)
% % % % dind0 = find(isinf(abs(grad))); % find absolute disontinuities (only flat ones)
% % % % dind1 = setdiff(dind,dind0);
% % % % 
% % % % %% clump discontinuities in strings. Fails if strong gradients near Moho
% % % % % % % if ~isempty(dind1)
% % % % % % %     a = diff(dind1);
% % % % % % %     b = find([a;inf]>1);
% % % % % % %     c = diff([0;b]);% length of sequences with discs
% % % % % % %     di1 = cumsum(c); % end points of sequences with discs
% % % % % % %     di0 = di1-c+1; % start points of sequences with discs
% % % % % % %     % re-sort indices back to full set of discontinuities)
% % % % % % %     [~,d1tod_,~] = intersect(dind,dind1);
% % % % % % %     di0 = d1tod_(di0);
% % % % % % %     di1 = d1tod_(di1);
% % % % % % %     % put abs discs back in
% % % % % % %     [~,d0tod_,~] = intersect(dind,dind0);
% % % % % % %     di0 = sort([d0tod_;di0]);
% % % % % % %     di1 = sort([d0tod_;di1]);
% % % % % % % else 
% % % % % % %     di0 = 1;di1 = 1;
% % % % % % % end
% % % % 
% % % % di0 = [1:length(dind)]; 
% % % % di1 = [1:length(dind)]; 
% % % % 
% % % % Ndisc = length(di0); % number of discontinuities
% % % % 
% % % % Zdisc = zeros(Ndisc,1);
% % % % Wdisc = zeros(Ndisc,1);
% % % % Vdto = zeros(Ndisc,1);
% % % % Vdbo = zeros(Ndisc,1);
% % % % kill = [];
% % % % for id = 1:Ndisc
% % % %     Zdisc(id) = mean(Zdo(dind(di0(id)):1+dind(di1(id))));
% % % %     Wdisc(id) = Zdo(dind(di1(id))+1) - Zdo(dind(di0(id)));
% % % %     Vdto(id) = Vdo(dind(di0(id)));
% % % %     Vdbo(id) = Vdo(dind(di1(id))+1);
% % % %     kill = [kill;find(Zdo>=Zdisc(id)-Wdisc(id)/2 & Zdo<=Zdisc(id)+Wdisc(id)/2)];
% % % % end
% % % % 
% % % % %% make discontinuous V and Z
% % % % % kill within zone of discontinuity
% % % % Zdo(kill) = []; 
% % % % Vdo(kill) = [];
% % % % % add in discontinuities explicitly and sort into place
% % % % Zdo = [Zdo;Zdisc;Zdisc];
% % % % Vdo = [Vdo;Vdto;Vdbo];
% % % % [Zdo,isort] = sort(Zdo);
% % % % Vdo = Vdo(isort);
% % % % % add on surface vel if needed
% % % % if min(Zdo)~=0, Zdo = [0;Zdo]; Vdo = [Vdto(1);Vdo];end 
% % % % 
% % % % %% ---------------  CONSTANT LAYERS  ---------------  %%
% % % % 
% % % % % extablish easy constant layers
% % % % dz = diff(Zdo);
% % % % dv = diff(Vdo);
% % % % grad = dv./dz;
% % % % cind = find(abs(grad) <= constgrad); % find constant layers
% % % % if ~isempty(cind)
% % % %     warning('brb2022.05.23 if you put any constant velocity layers here using constgrad, it breaks propmat and gives nonsense traces. Need to fix. ')
% % % % %!%! Removing the constant layer thing. 
% % % %     a = diff(cind);
% % % %     b = find([a;inf]>1);
% % % %     c = diff([0;b]);% length of constant sequences
% % % %     i1 = cumsum(c); % end points of constant sequences
% % % %     i0 = i1-c+1; % start points of constant sequences
% % % %     Nconst = length(i0); % number of constant layers
% % % % %!%! Removing the constant layer thing. 
% % % % %     Nconst=0; % todo later reimpliment the constgrad thing. Removed for now brb2022.05.23. 
% % % % else
% % % %     Nconst=0;
% % % % end
% % % % zlayt = zeros(Nconst,1);
% % % % zlayb = zeros(Nconst,1);
% % % % vlay = zeros(Nconst,1);
% % % % for id = 1:Nconst
% % % %     zlayt(id) = Zdo(cind(i0(id)));
% % % %     zlayb(id) = Zdo(cind(i1(id))+1);
% % % %     vlay(id) = mean(Vdo(cind(i0(id)):1+cind(i1(id))));
% % % % end
% % % %     
% % % % 
% % % % %% insert discontinuities as zero-thickness layers
% % % % [zlayt,isort] = sort([Zdisc;zlayt]);
% % % % zlayb = sort([Zdisc;zlayb]);
% % % % vtemp = [mean([Vdto,Vdbo],2);vlay];
% % % % vlay = vtemp(isort);
% % % % 
% % % % %% any seds/surface?
% % % % if zlayt(1)~=0 && min(zlayt)<7; % only account for sed layer if there is a discontinuity within 7 km of the surface
% % % %     zsb = zlayt(1);
% % % %     zlayt = [0;zsb/2;zlayt];
% % % %     zlayb = [zsb/2;zsb;zlayb];
% % % %     vlay = [Vdo(Z==0);Vdo(find(Zdo==zsb,1));vlay];
% % % % end
% % % % 
% % % % %% find GAPS not yet filled by layers, fill with steps
% % % % gp = find([zlayt;max(Z)] - [0;zlayb] > 0);
% % % % zlayt_temp = [zlayt;max(Z)];
% % % % zlayb_temp = [0;zlayb];
% % % % for igap = 1:length(gp) % loop through regions between discontinuities
% % % %     gpt = zlayb_temp(gp(igap));
% % % %     gpb = zlayt_temp(gp(igap));
% % % %     Zgp = Zdo(Zdo<=gpb & Zdo>=gpt);
% % % %     Vgp = Vdo(Zdo<=gpb & Zdo>=gpt);
% % % %     % kill wrong side of discontinuities if it's there
% % % %     if sum(Zgp==gpt)>1, Zgp(1) = [];Vgp(1) = []; end
% % % %     if sum(Zgp==gpb)>1, Zgp(end) = [];Vgp(end) = []; end
% % % % 
% % % %     %% split gaps into layers and take average V in each layer
% % % %     abscumdV = [0;cumsum(abs(diff(Vgp)))]; % take abs values in case some are negative
% % % %     cumdV = [0;cumsum(diff(Vgp))];
% % % %     nsplit = ceil(abscumdV(end)/dvmin); % nb true layers will be one more than this as we pin at beginning and end
% % % %     absdvsplit = [0:nsplit]*abscumdV(end)/nsplit; %#ok<*NBRAK>
% % % %     isplit = ones(nsplit+1,1); for ii = 2:nsplit+1, isplit(ii) = mindex(abscumdV,absdvsplit(ii)); end
% % % %     vsplit = Vgp(1) + cumdV(isplit);
% % % %     zsplit = mean([[Zgp(1);Zgp(isplit)],[Zgp(isplit);Zgp(end)]],2);
% % % %     zsplitt = zsplit(1:end-1);
% % % %     zsplitb = zsplit(2:end);
% % % %     
% % % %     % remove zero thickness layers
% % % %     wsplit = zsplitb-zsplitt;
% % % %     zsplitt(wsplit==0)=[];
% % % %     zsplitb(wsplit==0)=[];
% % % %     vsplit(wsplit==0)=[];
% % % % 
% % % %     zlayt = [zlayt;zsplitt];
% % % %     zlayb = [zlayb;zsplitb];
% % % %     vlay = [vlay;vsplit];
% % % % end
% % % % 
% % % % %% re-sort layer top & bottoms
% % % % [zlayt,isort] = sort(zlayt);
% % % % zlayb = sort(zlayb);
% % % % vlay = vlay(isort);
% % % % %% remove zero thickness layers (i.e. discontinuities)
% % % % wlay = zlayb-zlayt;
% % % % zlayt(wlay==0)=[];
% % % % zlayb(wlay==0)=[];
% % % % vlay(wlay==0)=[];
% % % % 
% % % % nlay = length(vlay);
% % % % 
% % % % %% Do to other variables, using linterp
% % % % for iv = 1:length(varargin)
% % % %     ival = varargin{iv}; 
% % % %     ivaldo = ival(:);
% % % %     izdo = Z(:);
% % % %    
% % % %     % do discs
% % % %     for id = 1:Ndisc
% % % %         ivaldto(id,1) = ivaldo(dind(id));
% % % %         ivaldbo(id,1) = ivaldo(dind(id)+1);
% % % %     end
% % % %     ivaldo(kill) = [];
% % % %     izdo(kill) = [];
% % % %     % add in discontinuities explicitly and sort into place
% % % %     izdo = [izdo;Zdisc;Zdisc];
% % % %     ivaldo = [ivaldo;ivaldto;ivaldbo];
% % % %     [izdo,isort] = sort(izdo);
% % % %     ivaldo = ivaldo(isort);
% % % %     % add on surface vals if needed
% % % %     if izdo(1)~=0, izdo = [0;izdo]; ivaldo = [ivaldto(1);ivaldo]; end
% % % %     
% % % % 
% % % %     oval = nan(size(vlay));
% % % %     zzz = [min(Z):0.1:max(Z)]; zzz=zzz(:);
% % % %     try
% % % %         val = linterp(izdo,ivaldo,zzz);
% % % %     catch
% % % %         error('Something wrong with linterp')
% % % %     end
% % % %     % first do surface
% % % % %     fprintf('check surface layers being layerised ok..\n')
% % % %     
% % % %     % special dispensation for zero-layer sediments at the top
% % % %     if find(Z==0,1,'last')>1
% % % %         izdo(1:find(Z==0,1,'last')-1) = [];
% % % %         ivaldo(1:find(Z==0,1,'last')-1) = [];
% % % %         ival(1:find(Z==0,1,'last')-1) = [];
% % % %     end
% % % %     
% % % %     oval(1) = ival(1); % was oval(1:2)=ival(1:2);
% % % %     % now average within each layer
% % % %     for ilay = 2:length(vlay) % was ilay = 3:length(...
% % % %         if any(zzz>zlayt(ilay) & zzz<zlayb(ilay)) % check at least one interped node in lay!
% % % %             oval(ilay) = mean(val(zzz>zlayt(ilay) & zzz<zlayb(ilay)));
% % % %         else % layer is thinner than 0.1 km
% % % %             oval(ilay) = linterp(izdo,ivaldo,mean([zlayt(ilay),zlayb(ilay)]));
% % % %         end
% % % %     end
% % % %  
% % % %     varargout{iv} = oval; %#ok<*AGROW>
% % % % end
% % % % 
% % % % %% plots
% % % % if ifplot
% % % %     figure(11); clf
% % % %     subplot(1,1+length(varargin),1),hold on
% % % %     % main var
% % % %     plot(V,Z,'-ko')
% % % % %     plot(Vdo,Zdo,'-ob')
% % % %     zlayp = reshape([zlayt';zlayb'],2*nlay,1);
% % % %     vlayp = reshape([vlay';vlay'],2*nlay,1);
% % % %     plot(vlayp,zlayp,'-ro')
% % % %     set(gca,'ydir','reverse','ylim',[0, max(Z)],'xlim',[0.9*min(V) 1.1*max(V)])
% % % %     % other vars
% % % %     for iv = 1:length(varargin)
% % % %         subplot(1,1+length(varargin),1+iv), hold on
% % % %         plot(varargin{iv},Z,'-ko')
% % % %         ivlayp = reshape([varargout{iv}';varargout{iv}'],2*nlay,1);
% % % %         plot(ivlayp,zlayp,'-ro')
% % % %         set(gca,'ydir','reverse','ylim',[0, max(Z)],'xlim',[0.9*min(varargin{iv}) 1.1*max(varargin{iv})])    
% % % %     end
% % % % end
% % % % 
% % % % end
% % % % 
% % % % 
% % % % function [ YI ] = linterp(X,Y,XI)
% % % % % YI = LINTERP(X,Y,XI) interpolates to find YI, the values of the
% % % % %     underlying function Y at the points in the array XI. X must be a
% % % % %     vector of length N.
% % % % % this function differs from the simple interp1 matlab function in that it
% % % % % can accept a vector x with multiple values of y (e.g. at the top of one
% % % % % layer and the bottom of another)
% % % % %
% % % % % Z. Eilon   May 2015
% % % % 
% % % % % use find_ilay to find where things fit - will not work properly for when
% % % % % a member of X matches a member of XI
% % % % [ilay] = find_ilay(XI,X);
% % % % 
% % % % % linearly interpolate
% % % % YI = (XI - X(ilay)).*(Y(ilay+1)-Y(ilay))./(X(ilay+1)-X(ilay)) + Y(ilay);
% % % % % now sort out coincident elements
% % % % olap = intersect(X,XI);
% % % % for i = 1:length(olap)
% % % % YI(XI==olap(i)) = mean(Y(X==olap(i)));
% % % % end
% % % % 
% % % % end
% % % % 
% % % % function [ minX_ind ] = mindex( X,a )
% % % % % [ minX_ind ] = mindex( X,a )
% % % % %   simple function to return the index of the minimum point in vector X
% % % % %   
% % % % %   if a second argument is given, the function outputs the index of the
% % % % %   point in X closest to the water level, a
% % % % %   basically just outputs the second output of the "min" function,
% % % % %   without giving you the magnitude of the minimum value.
% % % % % 
% % % % %   N.B. can be used as a zero finder if a==0
% % % % %
% % % % %   Intended for use when calling the value in one vector corresponding to
% % % % %   the minimum value in X - i.e more efficient than the clunkier:
% % % % %       Y(find(X==min(X)) or, more often, Y(find((X-a)==min(X-a)))
% % % % %   Instead, can now use
% % % % %       Y(mindex(X)) or Y(mindex(X,a))
% % % % % 
% % % % % Z. Eilon  
% % % % 
% % % % if nargin<2
% % % % [~,minX_ind] = min(X);
% % % % else
% % % % [~,minX_ind] = min(abs(X-a));
% % % % end
% % % % end
% % % % 
% % % % function [ilay] = find_ilay(r,Rb)
% % % % % [ilay] = find_ilay(r,Rb)
% % % % %
% % % % % Function to find the indices that each element of vector r would slot
% % % % % into vector Rb - originally conceived as a solution to the problem of
% % % % % having a series points at different radii and wanting to know which
% % % % % layers each of them were in, where the boundaries of the layers
% % % % % (including the top and bottom) are given by Rb. 
% % % % % 
% % % % % For example, if 3 layers were given by boundaries: [0;100;200;300] the
% % % % % point 50 would be in layer 1, and 217 would be in layer 3.
% % % % % thus, find_ilay([50;217],[0;100;200;300]) = [1;3]
% % % % %
% % % % % If a point in r is on a boundary, it is put into the upper layer, unless
% % % % % it is right at the max, then it is included in outermost layer
% % % % % 
% % % % % Z. Eilon  
% % % % 
% % % % 
% % % % if any(r < min(Rb)) || any(r > max(Rb))
% % % %     error('r must be within extremes of Rb');
% % % % end
% % % % 
% % % % r = r(:);
% % % % Rb = Rb(:);
% % % % 
% % % % N = length(r);
% % % % Nlay = length(Rb);
% % % % 
% % % % [~,ilay] = max((ones(N,1)*Rb' - r*ones(1,Nlay))>0,[],2);
% % % % ilay = ilay-1;
% % % % ilay(ilay==0) = Nlay-1; % if any are on outer edge, say in outermost layer
% % % % 
% % % % end

    