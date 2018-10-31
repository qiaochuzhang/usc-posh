function []=plotter_grid(switches,nsteps)

savefigs = switches(2);
NS = 9; %for now, number of stages is manually inserted; will be automated in the future.

if savefigs==1
    warning off;
    dirname=['figures','_np',num2str(nsteps)];
    mkdir(dirname);
    warning on;
end


% Figure parameters ======================================================%
def_size = 1; %switch to return to default figure size
LW  = 1; %linewidth for 2D plots
FS	= 14;  %fontsize for axis labels
FSN = 14;  %fontsize for tick mark labels
FST = 16;  %fontsize for figure title
frcx = 0.80; %fraction of screen length along x-direction
frcy = 0.65; %fraction of screen length along y-direction
rpy = 0.1; %range percentage for y-axis limits (shortcut)
fullscreen = get(0,'ScreenSize'); %get dimensions of monitor screen
pxi=(1-frcx)*fullscreen(3)/2; %lower left corner x-coordinate
pyi=(1-frcy)*fullscreen(4)/2; %lower left corner y-coordinate
pxs=frcx*fullscreen(3); %custom base dimension of figure
pys=frcy*fullscreen(4); %custom height dimension of figure
fps=1; %animation frames per second
quality=100; %animation quality percent

azim=90+45; %azimuthal viewing angle
pol=45; %polar viewing angle

load('seeddata.mat');

for k=1:nsteps
    
    %vary Wp
    if def_size==1; fig1=figure(1); else
    fig1=figure('Position',[pxi,pyi,pxs,pys]); end        
    hold on; grid on;
    v1=2; v2=3; v3=1;
    temp=permute(metric_grid,[v1,v2,v3]);
    p1=surf(var_arrays{v1},var_arrays{v2},temp(:,:,k)');
    axis([min(var_arrays{v1}),max(var_arrays{v1}),min(var_arrays{v2}),max(var_arrays{v2}),min(metric_grid(:)),max(metric_grid(:))]);
    view([azim,pol]);
    shading interp;
    % axis image;
    xlabel(['W_{n} (nm)']);
    ylabel(['W_{n,b} (nm)']);
    zlabel(['f_{0} (MHz)']);
    title({['VCO output frequency vs. transistor widths'];['N_{s} = ',...
        num2str(NS),', W_{p} = ',num2str(var_arrays{v3}(k)),' nm']});
    hold off;
    
    %vary Wn
    if def_size==1; fig2=figure(2); else
    fig2=figure('Position',[pxi,pyi,pxs,pys]); end        
    hold on; grid on; 
    v1=1; v2=3; v3=2;
    temp=permute(metric_grid,[v1,v2,v3]);
    p1=surf(var_arrays{v1},var_arrays{v2},temp(:,:,k)');
    axis([min(var_arrays{v1}),max(var_arrays{v1}),min(var_arrays{v2}),max(var_arrays{v2}),min(metric_grid(:)),max(metric_grid(:))]);
    view([azim,pol]);
    shading interp;
    % axis image;
    xlabel(['W_{p} (nm)']);
    ylabel(['W_{n,b} (nm)']);
    zlabel(['f_{0} (MHz)']);
    title({['VCO output frequency vs. transistor widths'];['N_{s} = ',...
        num2str(NS),', W_{n} = ',num2str(var_arrays{v3}(k)),' nm']});
    hold off;
    
    %vary Wnb
    if def_size==1; fig3=figure(3); else
    fig3=figure('Position',[pxi,pyi,pxs,pys]); end        
    hold on; grid on; 
    v1=1; v2=2; v3=3;
    temp=permute(metric_grid,[v1,v2,v3]);
    p1=surf(var_arrays{v1},var_arrays{v2},temp(:,:,k)');
    axis([min(var_arrays{v1}),max(var_arrays{v1}),min(var_arrays{v2}),max(var_arrays{v2}),min(metric_grid(:)),max(metric_grid(:))]);
    view([azim,pol]);
    shading interp;
    % axis image;
    xlabel(['W_{p} (nm)']);
    ylabel(['W_{n} (nm)']);
    zlabel(['f_{0} (MHz)']);
    title({['VCO output frequency vs. transistor widths'];['N_{s} = ',...
        num2str(NS),', W_{n,b} = ',num2str(var_arrays{v3}(k)),' nm']});
    hold off;
    
    if savefigs==1           
        for f=[1:3] %can generalize later...
            fs=num2str(f);
            fname=['fig',fs]; %filename
            if exist(fname)==1
                saveas(eval(fname),[dirname,'/',fname,'_k',num2str(k),'.fig']); %save as *.fig file
                saveas(eval(fname),[dirname,'/',fname,'_k',num2str(k),'.png']); %save as *.png file
            end
        end
        disp('Done saving figures!');
        close all;
    end
    
    clf;
    
end

disp('Generate VTK file (4D figure)'); %can generalize later...
mat2vtk([dirname,'/metric_grid.vtk'],metric_grid,'binary');

end

function mat2vtk(filename,matrix,format)
% Writes a 3D matrix as a *.VTK file as input for Paraview.
% Usage: mat2vtk('example.vtk',matrix,'binary');

% Get the matrix dimensions.
[Nx,Ny,Nz] = size(matrix);

% Open the file.
fid = fopen(filename,'w');
if fid == -1
    error('Cannot open file for writing.');
end

switch format
    case 'ascii'
        fprintf(fid,'# vtk DataFile Version 2.0\n');
        fprintf(fid,'Volume example\n');
        fprintf(fid,'ASCII\n');
        fprintf(fid,'DATASET STRUCTURED_POINTS\n');
        fprintf(fid,'DIMENSIONS %d %d %d\n',Nx,Ny,Nz);
        fprintf(fid,'ASPECT_RATIO %d %d %d\n',1,1,1);
        fprintf(fid,'ORIGIN %d %d %d\n',0,0,0);S
        fprintf(fid,'POINT_DATA %d\n',Nx*Ny*Nz);
        fprintf(fid,'SCALARS Pressure int 1\n');
        fprintf(fid,'LOOKUP_TABLE default\n');
        fwrite(fid, num2str(matrix(:)'));
    case 'binary'
        fprintf(fid,'# vtk DataFile Version 2.0\n');
        fprintf(fid,'Volume example\n');
        fprintf(fid,'BINARY\n');
        fprintf(fid,'DATASET STRUCTURED_POINTS\n');
        fprintf(fid,'DIMENSIONS %d %d %d\n',Nx,Ny,Nz);
        fprintf(fid,'ASPECT_RATIO %d %d %d\n',1,1,1);
        fprintf(fid,'ORIGIN %d %d %d\n',0,0,0);
        fprintf(fid,'POINT_DATA %d\n',Nx*Ny*Nz);
        fprintf(fid,'SCALARS Pressure float 1\n');
        fprintf(fid,'LOOKUP_TABLE default\n');
        fwrite(fid, matrix(:),'float','ieee-be');
    otherwise
        error('wrong input dummy :P');
end

% Close the file.
fclose(fid);

end