% ===========================================
%     Script to visualize all properties
% ===========================================
%clear all;
clc
% ===========================================
% paths
% ===========================================
test = 'stm_t001';
% ===========================================
stm_frame = 0;
% ===========================================
scripts   = ('/Volumes/SeismoLab/Subduction_sediments/scripts');
colormaps = ('/Volumes/SeismoLab/Subduction_sediments/scripts/colormaps/');
model     = (['/Volumes/SeismoLab/Subduction_sediments/',test]);
% ===========================================
cd(model);
clfall = 2;
if clfall == 0
    % keep figure
elseif clfall == 1
    clf;
elseif clfall == 2
    fh=findall(0,'type','figure');
    for i=1:length(fh)
        clf(fh(i));
        close(fh(i));
    end
end
% ===========================
% --- PLOTTING SPECIFICATIONS
% ===========================
expname     = 'stm';
setup       = 1;
nstart      = 70; 
nend        = nstart;
skip        = 1;
% ===========================
% --- PLOTTING SPECIFICATIONS
% ===========================
sp_var(1)   = 20;
sp_var(2)   = 0;
sp_var(3)   = 0;
% ========================================================================================================================
% 1-vx 2-vy 3-rho 4-temp 5-P 6-Sxx 7-Sxy 8-Nu 9-Exx 10-Exy 11-G 12-Exy_ne 13-G 14-S_ii 15-E_ii 16-E_ne_ii 17-S_ex 20-comp
% ===========================
% --- How to plot it ---
printmod    =  0;
zoom        =  1;
Clr_MAP     =  1;
print_max   =  0; 
% ===========================
% ===========================
% >>> contours
% ===========================
% Temperature contour
temp      = 0;
% composition contour
comp_cont = 0;
% velocity field
vel_field = 0;
% Isotherm location
eq_lines  = [100 200 300 400 500 600 ]; % (from Klingelhoefer, 2010)
%===========================
% Viscosity contours
visc      = 0;
%===========================
% Visco-plastic contours
VP_eii    = 0;
%===========================
% Stress axes
stress_axes = 0;
%===========================
% Tectonic pressure
tectonic_pressure     = 0;
relative_overpressure = 0;
%===========================
axis_wrt_trench =  0;    % 0 = as simulated, 1 = trench at (0,0)
rm_air          =  0;    % Set variable to zero in air, when zoom
% ==== large-scale settings ====
if setup == 1
    % trench coords from comp:
    x_trench = 1100;
    %x_trench = 2200;
    y_trench = 11;
    % Plotting limitations
    if zoom > 0

        
        x_beg = x_trench-200;
        x_end = x_trench+600;
        y_beg = y_trench-10;
        y_end = y_trench+200;
        
        
    end
    % Sticky Air Thickness from init.t3c
    SA_Th   = 20;
end
%==============================
if sp_var(2)==0 || sp_var(3)==0
    plot3var = 0;
else
    plot3var = 1;
end
set(0,'DefaultFigureColormap',feval('jet'))
% Dock figures
if printmod==2 || printmod>3
    set(0,'DefaultFigureWindowStyle','docked');
else
    % Keep as large as possible for paper figure resolution
    set(0,'DefaultFigureWindowStyle','normal');
end
%=================================
% Picture settings
if (exist('./Figures', 'dir') == 0)
    mkdir ./Figures
end
save_path = './Figures/';

if printmod == 3
    format = '-deps2c';
    res = '-r600';
else
    format = '-dpng'; %'-dtiff';
    %format = '-dtiff';  %'-dtiff'
    res = '-r400';
end
% 'new': new file numbering 001..010..999,'old': 1..10..999
fntype     = 'new'; %'old';
% Catch errors
if (nstart > nend)
    error('ERROR: "nend" should be superior to "nstart"!');
end
if (printmod==1 && zoom<=0)
    %error('ERROR in plotted figure: for plotting the non-high resolution area properly you need pcolor (not pcolor) and this is only possible in GUI ! ');
end
%====================
fonttouse = 'Arial';
m1  = 1;
mm1 = 1;
mm2 = 1;
%% Loop over all desired files
for nt = nstart:skip:nend
    % Build filenames for loading
    if (nt<10 && strcmp(fntype, 'new') == 1)
        filename = ([expname '00' num2str(nt) '.gzip.h5']);
    elseif (nt<100 && strcmp(fntype, 'new') == 1)
        filename = ([expname '0' num2str(nt) '.gzip.h5']);
    else
        filename = ([expname num2str(nt) '.gzip.h5']);
    end
    % If the file is missing open the next one
    limit = 0;
    while (exist(filename, 'file') == 0 && limit <4)
        
        disp(['File ', num2str(nt), ' is missing, jump to file ', num2str(nt+1)])
        nt = nt+1;
        if (nt<10 && strcmp(fntype, 'new') == 1)
            filename = ([expname '00' num2str(nt) '.gzip.h5']);
        elseif (nt<100 && strcmp(fntype, 'new') == 1)
            filename = ([expname '0' num2str(nt) '.gzip.h5']);
        else
            filename = ([expname num2str(nt) '.gzip.h5']);
        end
        %filename = (['data/' filename]);
        limit = limit + 1;
    end
    
    if(exist(filename, 'file') == 0)
        disp(['No files were found, check your data files. Exit...']);
        break;
    end
    
    % Read the structure of the hdf5 file
    fileinfo = hdf5info(filename);
    % disp(['Process file number ', num2str(nt)]);
    
    if nt == nstart
        x       = hdf5read(filename,'/ModelGroup/gx'); %x = x; %x = x(2:end)*1e-3;
        y       = hdf5read(filename,'/ModelGroup/gy'); y = y-SA_Th; %y = (y(2:end)*1e-3)-SA_Th;
        param   = hdf5read(filename,'/ModelGroup/Model');
        
        % Retrieves model parameters
        timesum  = param(1);
        xsize = param(2);
        ysize = param(3);
        xnumx = param(4);
        ynumy = param(5);
        % Build staggered nodes coordinates
        % To ensure that plot contents on right location of basic, plot them on basic node coordinates; not on staggered, central nodes coordinates
        % MATLAB already plots the coordinates you give in, in the center of its visualization cell
        x_stag  = 0.5*(x(1:end-1) + x(2:end));
        y_stag  = 0.5*(y(1:end-1) + y(2:end));
        
        % large-scale length to km
        if setup == 1
            x_stag = x_stag/1e3;
            y_stag = y_stag/1e3;
            xsize = xsize/1e3;
            ysize = ysize/1e3;
        end
        % calculate begin and end nodes for zoom visualization with pcolor
        % since pcolor can not handle irregular grids, but pcolor can not export nice enough pictures to AI etc
        if zoom > 0
            for ix = 1:1:xnumx
                if x_stag(ix) >= x_beg
                    nx_beg = ix;
                    break;
                end
            end
            for ix = nx_beg:1:xnumx
                if (ix==xnumx)
                    nx_end = ix-1;
                    break;
                elseif x_stag(ix) >= x_end
                    nx_end = ix;
                    break;
                end
            end
            for iy = 1:1:ynumy
                if (y_stag(iy)>=y_beg)
                    ny_beg = iy;
                    break;
                end
            end
            for iy = ny_beg:1:ynumy
                if (iy==ynumy)
                    ny_end = iy-1;
                    break;
                elseif (y_stag(iy)>=y_end)
                    ny_end = iy;
                    break;
                end
            end
        else
            nx_beg = 1;
            nx_end = xnumx-1;
            ny_beg = 1;
            ny_end = ynumy-1;
        end
        
        % Correct for plotting with respect to trench
        if axis_wrt_trench == 1
            x_stag = x_stag-x_trench;
            y_stag = y_stag-y_trench-SA_Th;
        end
        
        % Make marker grid for composition
        if sp_var(1) == 20 || sp_var(2) == 20 || sp_var(3) == 20 || (rm_air==1 && setup==1)
            % Get minimum gridsize in km, since this is resolution visualization grid
            res_high = (min(y_stag(2:end)-y_stag(1:end-1)));
            
            % Get nr elements for visualization grids
            mxvislim = floor(xsize/(res_high) + 1);
            myvislim = floor(ysize/(res_high) + 1);
            
            % Composition matrix: size of dense regular markergrid
            % Whole for filling matrix correctly; zoom while plotting
            composition = ones(myvislim,mxvislim);
            
            % Calculate begin and end nodes for visualization of composition with pcolor
            % since pcolor can not handle irregular grids, but pcolor can not export nice enough pictures to AI etc
            if zoom > 0
                mnx_beg = floor(x_beg/res_high);
                mnx_end = floor(x_end/res_high);
                mny_beg = floor(y_beg/res_high);
                mny_end = floor(y_end/res_high);
            else
                mnx_beg = 1;
                mnx_end = floor(xsize/res_high);
                mny_beg = 1;
                mny_end = floor(ysize/res_high);
            end
            
            % Make marker grid
            mx = mnx_beg:1:mnx_end;
            my = mny_beg:1:mny_end;
            
            if axis_wrt_trench == 1
                mx = mx*res_high-x_trench;
                my = my*res_high-y_trench-SA_Th;
            elseif axis_wrt_trench == 0
                mx = mx*res_high;
                my = my*res_high;
            end
            
        end
        
    else
        param   = hdf5read(filename,'/ModelGroup/Model');
        timesum  = param(1);
    end
    
    if temp==1
        T = hdf5read(filename,'/NodeGroup/tk');
        T = reshape(T, ynumy, xnumx);         % on basic nodes
    end
    if visc ==1
        nu      = hdf5read(filename,'/NodeGroup/nu');
        nu    = reshape(   nu, ynumy, xnumx);  %on basic nodes
    end
    
    if VP_eii == 1
        sxx_tmp = hdf5read(filename,'/NodeGroup/sxx');
        nd_tmp  = hdf5read(filename,'/NodeGroup/nd');
        e_nexx  = sxx_tmp(1:end,1:end)./(2*nd_tmp(1:end,1:end));
        
        sxy_tmp = hdf5read(filename,'/NodeGroup/sxy');
        nu_tmp  = hdf5read(filename,'/NodeGroup/nu');
        e_nexy  = sxy_tmp(1:end,1:end)./(2*nu_tmp(1:end,1:end));
        
        e_nexx  = reshape(e_nexx, ynumy, xnumx);
        e_nexy  = reshape(e_nexy, ynumy, xnumx);
        
        e_nexx  = e_nexx(1:end-1,1:end-1);     % on staggered nodes
        e_nexy  = e_nexy(1:end,1:end);         % on basic nodes
        
        % Interpolates shear tensor components to staggered nodes, so visualize center of computational cell
        e_nexy  = 0.25*(e_nexy(1:end-1, 1:end-1) + e_nexy(2:end, 2:end) +  e_nexy(2:end, 1:end-1) + e_nexy(1:end-1, 2:end));
        
        % Calculate second invariant
        vp_Eii  = sqrt(e_nexx.^2 + e_nexy.^2);
    end
    
    
    if plot3var == 1
        last_iv = 3;%3;
    else
        last_iv = 1;
    end
    
    for iv = 1:1:last_iv
        vis_var = sp_var(iv);
        
        if vis_var == 1
            % make in cm/year
            var     = hdf5read(filename,'/NodeGroup/vx');
            if setup==1
                var = (var*1e2)*60*60*24*365.25;
            end
            var = reshape(var, ynumy, xnumx);
        elseif vis_var == 2
            var     = hdf5read(filename,'/NodeGroup/vy');
            % switch around, so - is subsidence
            %var = -var;
            if setup==1
                var = (var*1e2)*60*60*24*365.25;
            end
            var = reshape(var, ynumy, xnumx);
        elseif vis_var == 3
            var     = hdf5read(filename,'/NodeGroup/ro');
            var = reshape(var, ynumy, xnumx);
        elseif vis_var == 4
            var     = hdf5read(filename,'/NodeGroup/tk');
            var = reshape(var, ynumy, xnumx);
            var = var-273;
        elseif vis_var == 5
            var     = hdf5read(filename,'/NodeGroup/pr');
            var = reshape(var, ynumy, xnumx);
        elseif vis_var == 6
            var     = hdf5read(filename,'/NodeGroup/sxx');
            var     = reshape(var, ynumy, xnumx);
            var     = var./1e6;
        elseif vis_var == 7
            var     = hdf5read(filename,'/NodeGroup/sxy');
            var = reshape(var, ynumy, xnumx);
        elseif vis_var == 8
            var     = hdf5read(filename,'/NodeGroup/nu');
            var = reshape(var, ynumy, xnumx);
        elseif vis_var == 9
            var     = hdf5read(filename,'/NodeGroup/exx');
            var = reshape(var, ynumy, xnumx);
        elseif vis_var == 10
            var     = hdf5read(filename,'/NodeGroup/exy');
            var = reshape(var, ynumy, xnumx);
        elseif vis_var == 11
            var     = hdf5read(filename,'/NodeGroup/gg');
            var = reshape(var, ynumy, xnumx);
            %exx_ne
        elseif vis_var == 12
            sxx_tmp = hdf5read(filename,'/NodeGroup/sxx');
            nd_tmp  = hdf5read(filename,'/NodeGroup/nd');
            
            sxx_tmp = reshape(sxx_tmp, ynumy, xnumx);
            nd_tmp  = reshape(nd_tmp, ynumy, xnumx);
            
            var     = sxx_tmp(1:end,1:end)./(2*nd_tmp(1:end,1:end));
            %exy_ne
        elseif vis_var == 13
            sxy_tmp = hdf5read(filename,'/NodeGroup/sxy');
            nu_tmp  = hdf5read(filename,'/NodeGroup/nu');
            
            sxy_tmp = reshape(sxy_tmp, ynumy, xnumx);
            nu_tmp  = reshape(nu_tmp, ynumy, xnumx);
            
            var     = sxy_tmp(1:end,1:end)./(2*nu_tmp(1:end,1:end));
            % sii
        elseif vis_var == 14
            sxx_tmp = hdf5read(filename,'/NodeGroup/sxx');
            sxy_tmp = hdf5read(filename,'/NodeGroup/sxy');
            
            sxx_tmp = reshape(sxx_tmp, ynumy, xnumx);
            sxy_tmp = reshape(sxy_tmp, ynumy, xnumx);
            
            sxx_tmp = sxx_tmp(1:end-1,1:end-1);     % on staggered nodes
            sxy_tmp = sxy_tmp(1:end,1:end);         % on basic nodes
            
            % Interpolates shear tensor components to staggered nodes, so visualize center of computational cell
            sxy_tmp   = 0.25*(sxy_tmp(1:end-1, 1:end-1) + sxy_tmp(2:end, 2:end) +  sxy_tmp(2:end, 1:end-1) + sxy_tmp(1:end-1, 2:end));
            
            % Calculate second invariant
            var   = sqrt(sxx_tmp.^2 + sxy_tmp.^2);
            
        elseif vis_var == 15
            exx_tmp = hdf5read(filename,'/NodeGroup/exx');
            exy_tmp = hdf5read(filename,'/NodeGroup/exy');
            
            exx_tmp = reshape(exx_tmp, ynumy, xnumx);
            exy_tmp = reshape(exy_tmp, ynumy, xnumx);
            
            exx_tmp = exx_tmp(1:end-1,1:end-1);     % on staggered nodes
            exy_tmp = exy_tmp(1:end,1:end);         % on basic nodes
            
            % Interpolates shear tensor components to staggered nodes, so visualize center of computational cell
            exy_tmp   = 0.25*(exy_tmp(1:end-1, 1:end-1) + exy_tmp(2:end, 2:end) +  exy_tmp(2:end, 1:end-1) + exy_tmp(1:end-1, 2:end));
            
            % Calculate second invariant
            var   = sqrt(exx_tmp.^2 + exy_tmp.^2);
            %var   = var*8;
            %var(var>1e-15)=1e-10;
            %eii_ne
        elseif vis_var == 16
            sxx_tmp = hdf5read(filename,'/NodeGroup/sxx');
            nd_tmp  = hdf5read(filename,'/NodeGroup/nd');
            e_nexx  = sxx_tmp(1:end,1:end)./(2*nd_tmp(1:end,1:end));
            
            sxy_tmp = hdf5read(filename,'/NodeGroup/sxy');
            nu_tmp  = hdf5read(filename,'/NodeGroup/nu');
            e_nexy  = sxy_tmp(1:end,1:end)./(2*nu_tmp(1:end,1:end));
            
            e_nexx  = reshape(e_nexx, ynumy, xnumx);
            e_nexy  = reshape(e_nexy, ynumy, xnumx);
            
            e_nexx  = e_nexx(1:end-1,1:end-1);     % on staggered nodes
            e_nexy  = e_nexy(1:end,1:end);         % on basic nodes
            
            % Interpolates shear tensor components to staggered nodes, so visualize center of computational cell
            e_nexy  = 0.25*(e_nexy(1:end-1, 1:end-1) + e_nexy(2:end, 2:end) +  e_nexy(2:end, 1:end-1) + e_nexy(1:end-1, 2:end));
            
            % Calculate second invariant
            var   = sqrt(e_nexx.^2 + e_nexy.^2);
        elseif vis_var == 17
            sbrit_tmp = hdf5read(filename,'/NodeGroup/sbrit');
            sbrit_tmp = reshape(sbrit_tmp, ynumy, xnumx);
            sbrit_tmp = sbrit_tmp(1:end-1,1:end-1);
            sxx_tmp = hdf5read(filename,'/NodeGroup/sxx');
            sxy_tmp = hdf5read(filename,'/NodeGroup/sxy');
            sxx_tmp = reshape(sxx_tmp, ynumy, xnumx);
            sxy_tmp = reshape(sxy_tmp, ynumy, xnumx);
            sxx_tmp = sxx_tmp(1:end-1,1:end-1);     % on staggered nodes
            sxy_tmp = sxy_tmp(1:end,1:end);         % on basic nodes
            
            % Interpolates shear tensor components to staggered nodes, so visualize center of computational cell
            sxy_tmp   = 0.25*(sxy_tmp(1:end-1, 1:end-1) + sxy_tmp(2:end, 2:end) +  sxy_tmp(2:end, 1:end-1) + sxy_tmp(1:end-1, 2:end));
            % Calculate second invariant
            sii   = sqrt(sxx_tmp.^2 + sxy_tmp.^2);
            var   = (sbrit_tmp(1:end,1:end)-sii(1:end,1:end))./1e+6;
            
        elseif vis_var == 18
            Markrh = hdf5read(filename,'/NodeGroup/markrh');
            markrh = reshape(Markrh, ynumy, xnumx);
            var   = markrh;
            
        end
        if vel_field
            % make in cm/year
            vx     = hdf5read(filename,'/NodeGroup/vx');
            vx = (vx*1e2)*60*60*24*365.25;
            vx = reshape(vx, ynumy, xnumx);
            vx = vx(1:end-1,2:end);     % on ghost nodes
            % =========================================
            vy     = hdf5read(filename,'/NodeGroup/vy');
            vy = (vy*1e2)*60*60*24*365.25;
            vy = reshape(vy, ynumy, xnumx);
            vy = vy(2:end,1:end-1);     % on ghost node
        end
        %comp
        if vis_var == 20 || (rm_air==1 && setup==1 && iv==1)
            mtype = hdf5read(filename,'/VisMarkerGroup/Mtype');
            cd(colormaps);
            if Clr_MAP == 1
                disp('load colormap')
                load colorMAP_V2.mat;
            else
                load colorMAPgray;
            end
            cd(model);
            % Extract composition from efficient storage
            num       = 1;
            ind       = 1;
            while num<length(mtype)
                value = mtype(num);
                
                if value==-2
                    % Compressed: the next ?? indices are given the color material
                    num_colors  =   mtype(num+1);
                    material    =   mtype(num+2);
                    ind_vec     =   ind:ind+num_colors-1;
                    
                    ind         =   ind_vec(end)+1;
                    num         =   num+3;
                elseif value==-1
                    ind_vec     =   ind;
                    material    =   0;
                    ind         =   ind_vec(end)+1;
                    num         =   num+1;
                else
                    ind_vec     =   ind;
                    material    =   value;
                    
                    ind         =   ind_vec(end)+1;
                    num         =   num+1;
                end
                
                composition(ind_vec) =   material;
            end
        end
        
        % Remove velocities in air for paper plotting to avoid distraction
        if rm_air==1
            if setup==1
                if zoom>0
                    for i=0:1:length(mx)
                        for j=0:1:length(my)
                            % adjust composition and var arrays such that cover same area
                            if composition(mny_beg+j,mnx_beg+i)==0
                                var(ny_beg+j,nx_beg+i) = NaN;
                            end
                        end
                    end
                else
                    display('NOTE: Can not remove air for no zoom (need all within HR grid): will never do for papers anyway!');
                end
            end
        end
        
        
        % Fill corresponding variable for subplotting 3 variables
        if plot3var == 1
            if iv==1 && sp_var(iv)~=20
                var1 = var;
            elseif iv==2 && sp_var(iv)~=20
                var2 = var;
            elseif iv==3 && sp_var(iv)~=20
                var3=var;
            end
        end
    end
    
    %====================================
    %     Plot figure ===================
    %====================================
    if printmod == 1
        figure('Visible', 'Off');
    elseif printmod == 3
        figure('units','normalized','outerposition',[0 0 1 1]);
    else
        figure(1);
        clf;
    end
    
    for iv = 1:1:last_iv
        % Load appropriate variables for this iv
        vis_var = sp_var(iv);
        
        if plot3var == 1
            subplot(3,1,iv);
            
            if iv == 1 && sp_var(iv)~=20
                var = var1;
            elseif iv == 2 && sp_var(iv)~=20
                var = var2;
            elseif iv == 3 && sp_var(iv)~=20
                var = var3;
            end
        end
        
        % Log variables
        if vis_var==8 || vis_var==14 || vis_var==15 || vis_var==16
            %colormap('default');
            cd(colormaps)
            %colormap redblue;
            cd(model)
            if printmod==2 || printmod>3
                pcolor(x_stag(nx_beg:nx_end), y_stag(ny_beg:ny_end), log10(var(ny_beg:ny_end,nx_beg:nx_end)));
            elseif printmod==0 || printmod==1 || printmod==3
                pcolor(x_stag(nx_beg:nx_end), y_stag(ny_beg:ny_end), log10(var(ny_beg:ny_end,nx_beg:nx_end)));
                %pcolor(x_stag(nx_beg:nx_end), y_stag(ny_beg:ny_end),var(ny_beg:ny_end,nx_beg:nx_end));
                %pcolor(x_stag(nx_beg:nx_end), y_stag(ny_beg:ny_end)-20, log10(var(ny_beg:ny_end,nx_beg:nx_end)));
            end
            %freezeColors;
            cd(scripts)
            %load('viridis','viridis');
            cd(model)
            colormap; %(viridis);
            cb=colorbar('SouthOutside');
        elseif vis_var==17
            %colormap('default');
            cd(colormaps)
            colormap redblue;
            cd(model)
            if printmod==2 || printmod>3
                pcolor(x_stag(nx_beg:nx_end), y_stag(ny_beg:ny_end), (var(ny_beg:ny_end,nx_beg:nx_end)));
            elseif printmod==0 || printmod==1 || printmod==3
                pcolor(x_stag(nx_beg:nx_end), y_stag(ny_beg:ny_end), (var(ny_beg:ny_end,nx_beg:nx_end)));
            end
            %cbrb=colorbar;
            % Composition
        elseif vis_var==18
            cd(colormaps)
            colormap redblue;
            cd(model)
            if printmod==2 || printmod>3
                pcolor(x_stag(nx_beg:nx_end), y_stag(ny_beg:ny_end), (var(ny_beg:ny_end,nx_beg:nx_end)));
            elseif printmod==0 || printmod==1 || printmod==3
                pcolor(x_stag(nx_beg:nx_end), y_stag(ny_beg:ny_end), (var(ny_beg:ny_end,nx_beg:nx_end)));
            end
            cbrb=colorbar;
            
        elseif vis_var == 20
            if plot3var == 0
                if Clr_MAP == 1
                    colormap(colorMAP);
                else
                    colormap(colorMAPgray);
                end
            end
            if axis_wrt_trench == 1
                imagesc(mx, my, composition(mny_beg:1:mny_end,mnx_beg:1:mnx_end));
            elseif axis_wrt_trench == 0
                imagesc(mx, my, composition(mny_beg:1:mny_end,mnx_beg:1:mnx_end));
            end
            %freezeColors;
            % All other variables
        elseif vel_field == 1
            if vis_var == 1
                hold on
                cd(colormaps)
                %colormap redblue;
                cd(model)
                if printmod==2 || printmod>3
                    pcolor(x_stag(nx_beg:nx_end), y_stag(ny_beg:ny_end), (var(ny_beg:ny_end,nx_beg:nx_end)));
                elseif printmod==0 || printmod==1 || printmod==3
                    pcolor(x_stag(nx_beg:nx_end), y_stag(ny_beg:ny_end), (var(ny_beg:ny_end,nx_beg:nx_end)));
                end
                %cb=colorbar('SouthOutside');
                %cb=colorbar;
                hold on
                stp  = 11;
                stpy = 21;
                quiver(x_stag(nx_beg+stp:stp:nx_end-6),y_stag(ny_beg+stpy:stp:ny_end-6),vx(ny_beg+stpy:stp:ny_end-6,...
                    nx_beg+stp:stp:nx_end-6)*3,vy(ny_beg+stpy:stp:ny_end-6,nx_beg+stp:stp:nx_end-6)*3,'AutoScale','off','color','y');
                axis([x_stag(nx_beg) x_stag(nx_end) y_stag(ny_beg) y_stag(ny_end)]);
            end
        else
            cd(colormaps)
            colormap hot;%redblue;
            cd(model)
            %colormap('autumn');
            %colormap('default');
            if printmod > 1
                pcolor(x_stag(nx_beg:nx_end), y_stag(ny_beg:ny_end), var(ny_beg:ny_end,nx_beg:nx_end));
            else
                pcolor(x_stag(nx_beg:nx_end), y_stag(ny_beg:ny_end), var(ny_beg:ny_end,nx_beg:nx_end));
            end
            %freezeColors;
            %cb=colorbar;
        end
        
        if comp_cont %composition contour
            hold on
            eq_comp = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18];
            eq_com2 = [0 1];
            [C3, h3] = contour(mx, my, composition(mny_beg:1:mny_end,mnx_beg:1:mnx_end),eq_comp,'color','k','Linewidth',0.03);
            hold off
        end
        
        % Isotherms
        if temp == 1
            hold on;
            [C2,h2] = contour(x_stag(nx_beg:nx_end), y_stag(ny_beg:ny_end), T(ny_beg:ny_end,nx_beg:nx_end)-273.15,eq_lines,'r','Linewidth',1);
            clabel(C2,h2,'color','r','fontsize',10,'rotation',0,'LabelSpacing',500);
            axis ij image;
            hold off;
        end
        % Viscosity
        if visc == 1
            hold on;
            [C3,h3] = contour(x_stag(nx_beg:nx_end),y_stag(ny_beg:ny_end),log10(nu(ny_beg:ny_end,nx_beg:nx_end)),[24 24],'w','Linewidth',0.3);
            axis ij image;
            hold off
        end
        % Visco-plastic
        if VP_eii == 1
            nY_end = 200;
            hold on;
            [C4,h4] = contour(x_stag(nx_beg:nx_end),y_stag(ny_beg:nY_end),log10(vp_Eii(ny_beg:nY_end,nx_beg:nx_end)),[-16 -16],'r','Linewidth',0.3);
            axis ij image;
            hold off
        end
        
        
        %===============================================
        % === Set visual properties
        %===============================================
        if setup==1
            
            if (vis_var < 20)
                cb=colorbar('SouthOutside');
                set(cb, 'Position', [.684 .750 0.22 .0180],'LineWidth',1.3);
                set(cb, 'FontSize', 10, 'FontName', 'Avenir','fontweight','bold');
            end
            
            
            if vis_var == 1
                %ylabel(cb,'Vx (cm/yr)','FontSize', 10, 'FontName', 'Arial','fontweight','bold');
                ylabel(cb,'Horizontal velocity (cm/yr)');
                %vx_lim = 1.5e-9;
                %caxis([-1 +1]);
            elseif vis_var == 2 && vel_field == 1
                ylabel(cbrb,'Upward velocity [cm/yr]')
                caxis ([-5 5]);
            elseif vis_var == 2
                ylabel(cb,'Vertical velocity [cm/yr]')
                caxis([-1 +1]);
            elseif vis_var == 3
                ylabel(cb,'Density [kg/m^3]');
                caxis ([2500 5500]);
            elseif vis_var == 4
                ylabel(cb,'Temperature (^{\circ}C)');
                %caxis ([400 700]);
                caxis ([0 1400]);
            elseif vis_var == 5
                ylabel(cb,'Pressure [Pa]');
                %caxis ([1.3e9 1.4e9]);
                if tectonic_pressure
                    caxis ([-10.9113e+08 +1.9113e+08]);
                end
                if relative_overpressure
                    caxis ([-50 +50]);
                end
            elseif vis_var == 6
                ylabel(cb,'Normal stress (MPa)');
                caxis ([-200 200]);
            elseif vis_var == 7
                ylabel(cb,'Dev. shear stress [Pa]');
                caxis ([-30e6 30e6]);
            elseif vis_var == 8
                %caxis ([17 25]);
                caxis ([20 25]);
                ylabel(cb,'log_{10}(\eta_{eff}) (Pa s)','FontSize', 10, 'FontName', 'Arial','fontweight','bold');
            elseif vis_var == 9
                ylabel(cb,'\epsilon_{xx} [s^{-1}]');
                caxis ([-2e-14 2e-14]);
            elseif vis_var == 10
                ylabel(cb,'\epsilon_{xy} [s^{-1}]');
                caxis ([-2e-14 2e-14]);
            elseif vis_var == 11
                ylabel(cb,'shear modulus [Pa]');
            elseif vis_var == 12
                ylabel(cb,'\epsilon_{vp, xx} [s^{-1}]');
                caxis ([-2e-14 2e-14]);
            elseif vis_var == 13
                ylabel(cb,'\epsilon_{vp, xy} [s^{-1}]');
                caxis ([-2e-14 2e-14]);
            elseif vis_var == 14
                ylabel(cb,'log_{10}(stress II) (MPa)','FontSize', 10, 'FontName', 'Arial','fontweight','bold');
                caxis(gca,[6 10]);
            elseif vis_var == 15
                ylabel(cb,'log10(strain rate II) [s^{-1}]');
                caxis(gca,[-15 -14]);
            elseif vis_var == 16
                ylabel(cb,'log_{10}(\epsilon_{vp, ii}) [s^{-1}]');
                caxis(gca,[-16 -14]);
            elseif vis_var == 17
                caxis(gca,[0 10]);
            elseif vis_var == 20
                if plot3var == 1
                    caxis ([1 10]);
                else
                    caxis ([0 36]);
                end
            end
        end
        
        if vis_var==20
            shading flat
            cbh = findobj( 0,'tag','Colorbar');
            delete(cbh)
            shading interp
        end
        
        if zoom > 0 && sp_var(1) == 20 || sp_var(2) == 20 || sp_var(3) == 20
        elseif zoom < 0 && sp_var(1) == 20 || sp_var(2) == 20 || sp_var(3) == 20
            text(3800,1250,[num2str(timesum/1e6,4),' Myr'],'FontName','Arial','fontweight','bold','HorizontalAlignment','right','color','w','FontSize',13)
        elseif zoom > 0 && sp_var(2) == 17
            text(x_end-20,y_end-20,[num2str(timesum/1e6,4),' Myr'],'FontName','Arial','fontweight','bold','HorizontalAlignment','right','color','w','FontSize',13)
        elseif zoom > 0 && sp_var(2) == 17
            text(x_end-20,y_end-20,[num2str(timesum/1e+6,4),' Myr'],'FontName','Arial','fontweight','bold','HorizontalAlignment','right','color','w','FontSize',13)
        end
        
        
        shading interp
        axis ij image
        set(gca,'FontSize', 12, 'FontName', 'Avenir','fontweight','bold');
        xlabel('Distance (km)','FontSize', 12, 'FontName', 'Avenir','fontweight','bold');
        ylabel('Depth (km)','FontSize', 12, 'FontName', 'Avenir','fontweight','bold');
        set(gca,'XMinorTick','on','YMinorTick','on');
        set(gca,'xTick',0:100:4000);
        set(gca,'yTick',0:50:1400);
        
        if zoom == 0
            set(gca,'xTick',0:500:4000);
            set(gca,'yTick',0:150:1400);
            xlabel('Distance (km)','FontSize', 12, 'FontName', 'Avenir','fontweight','bold');
        end
        
        set(gca, 'Layer', 'top')

        
        % Add title
        if (iv == 1)
            if setup == 1
                tlt = title(['Time: ', num2str(timesum/1e6),' Myr']);
                set(tlt,'FontSize', 12, 'FontName', 'Avenir','fontweight','bold');
            end
        end
    end
    
    


    % ---------------------------------------------------------------
    % --- Save figures ----------------------------------------------
    % ---------------------------------------------------------------
    if printmod >= 1
        % Build up filename
        if vis_var == 1
            fig_var = 'vx';
        elseif vis_var == 2
            fig_var = 'vy';
        elseif vis_var == 3
            fig_var = 'rho';
        elseif vis_var == 4
            fig_var = 'tk';
        elseif vis_var == 5
            fig_var = 'pr';
        elseif vis_var == 6
            fig_var = 'sxx';
        elseif vis_var == 7
            fig_var = 'sxy';
        elseif vis_var == 8
            fig_var = 'nu';
        elseif vis_var == 9
            fig_var = 'exx';
        elseif vis_var == 10
            fig_var = 'exy';
        elseif vis_var == 11
            fig_var = 'gg';
        elseif vis_var == 12
            fig_var = 'exx_ne';
        elseif vis_var == 13
            fig_var = 'exy_ne';
        elseif vis_var == 14
            fig_var = 'sii';
        elseif vis_var == 15
            fig_var = 'eii';
        elseif vis_var == 16
            fig_var = 'eii_ne';
        elseif vis_var == 17
            fig_var = 'strength_excess';
        elseif vis_var == 18
            fig_var = 'markRH';
        elseif vis_var == 20
            fig_var = 'comp';
        end
        if zoom > 0 && plot3var == 0
            ext = '_z_';
        elseif zoom > 0 && plot3var == 1
            ext = '_zsp_';
        elseif zoom <= 0 && plot3var == 0
            ext = '_';
        elseif zoom <= 0 && plot3var == 1
            ext = '_sp_';
        end
        
        if nt<10
            filename_fig    =  [save_path, expname,ext,fig_var,'00', num2str(nt)];
        elseif nt<100 &&  nt>=10
            filename_fig    =  [save_path, expname,ext,fig_var,'0', num2str(nt)];
        else
            filename_fig    =  [save_path, expname,ext,fig_var, num2str(nt)];
        end
        % Note adds last variable to filename in case of subplot
        
        % Can uncomment this print comment if testing many figures
        print(format,res,filename_fig);
        % For paper quality figures:
        % Note pcolor and deps2c allow you to nicely edit figures in Illustrator (prior tp CS6)
        
        % Data save so that can very nice plots in a different script
        % if desired (if printmod=3 quality or complexity is not enough)
        if printmod==3 && nt==nstart
            if zoom <= 0
                file_save = ['hq_fig_data_',fig_var,'_',char(expname),num2str(nt),'.mat'];
            else
                file_save = ['hq_fig_data_z_',fig_var,'_',char(expname),num2str(nt),'.mat'];
            end
            if vis_var==20
                if setup == 1
                    %save(file_save,'mx','my','composition','T','x','y','x_stag','y_stag','res_high','y_trench','x_trench','SA_Th','nx_beg','ny_beg','nx_end','ny_end','mnx_beg','mny_beg','mnx_end','mny_end','eq_lines');
                end
            else
                if setup == 1
                elseif temp == 1
                    save(file_save,'var','T','x','y','x_stag','y_stag','y_trench','x_trench','SA_Th','nx_beg','ny_beg','nx_end','ny_end','eq_lines');
                elseif temp == 0
                    save(file_save,'var','x','y','x_stag','y_stag','y_trench','x_trench','SA_Th','nx_beg','ny_beg','nx_end','ny_end','eq_lines');
                else
                    save(file_save,'var','x','y','x_stag','y_stag','y_trench','x_trench','nx_beg','ny_beg','nx_end','ny_end');
                end
            end
        end
        
        if nt~=nend
            %close(gcf);
        end
        %disp(['Image ', num2str(nt), ' printed and saved.'])
    end
    
end





