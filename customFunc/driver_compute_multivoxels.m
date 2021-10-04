%close all;


%   femesh: struct with fields
%       ncompartment: [1 x 1]
%       nboundary: [1 x 1]
%       points: {1 x ncompartment}[3 x npoint(icmpt)]
%       facets: {ncompartment x nboundary}[3 x nfacet(icmpt, ibdry)]
%       elements: {1 x ncompartment}[4 x nelement(icmpt)]
%       point_map: {1 x ncompartment}[1 x npoint(icmpt)]


for icmpt = 1:ncompartment
    xmin(icmpt) = min(min(min(femesh.points{icmpt}(1,:))));
    ymin(icmpt) = min(min(min(femesh.points{icmpt}(2,:))));
    zmin(icmpt) = min(min(min(femesh.points{icmpt}(3,:))));
    
    xmax(icmpt) = max(max(max(femesh.points{icmpt}(1,:))));
    ymax(icmpt) = max(max(max(femesh.points{icmpt}(2,:))));
    zmax(icmpt) = max(max(max(femesh.points{icmpt}(3,:))));
end

xxmin = min(xmin);  yymin = min(ymin);  zzmin = min(zmin);
xxmax = max(xmax);  yymax = max(ymax);  zzmax = max(zmax);

voxels_nx = 4; voxels_ny = 2; voxels_nz = 2;
voxels_xlim = linspace(xxmin,xxmax,voxels_nx);
voxels_ylim = linspace(yymin,yymax,voxels_ny);
voxels_zlim = linspace(zzmin,zzmax,voxels_nz);

%voxels_xlim = [xxmin,-25, 25, xxmax];
%voxels_ylim = [yymin,-25, 25, yymax];
%voxels_zlim = [zzmin,-25, 25, zzmax];

voxels_xlim = linspace(xxmin,xxmax,3);
voxels_ylim = linspace(yymin,yymax,3);
voxels_zlim = linspace(zzmin,zzmax,3);

voxels_xlim = [xxmin:1:xxmax];
voxels_ylim = [yymin:1:yymax];
voxels_zlim = [zzmin:1:zzmax];

voxels_nx = length(voxels_xlim);
voxels_ny = length(voxels_ylim);
voxels_nz = length(voxels_zlim);


% voxels_ylim = [yymin,25,yymax];
% voxels_ny = length(voxels_ylim);

iv = 0;
ng=size(btpde.magnetization,4);
nsimul = size(btpde.magnetization,2);
Signal=zeros(voxels_nx-1,voxels_ny-1,voxels_nz-1,nsimul,ng,ncompartment);
clear signal_voxel volumes_voxel signal_voxel_normalize;
p=progressbar(prod([voxels_nx-1,voxels_ny-1,voxels_nz-1]));
parfor ii = 1:voxels_nx-1
    for jj = 1:voxels_ny-1
        for kk = 1:voxels_nz-1
            %iv = iv+1;
            percent=p.progress;
            xmin = voxels_xlim(ii);  xmax = voxels_xlim(ii+1);
            ymin = voxels_ylim(jj);  ymax = voxels_ylim(jj+1);
            zmin = voxels_zlim(kk);  zmax = voxels_zlim(kk+1);
            for icmpt = 1:ncompartment
                ll = find(femesh.points{icmpt}(1,:)>=xmin ...
                    & femesh.points{icmpt}(1,:)<=xmax ...
                    & femesh.points{icmpt}(2,:)>=ymin ...
                    & femesh.points{icmpt}(2,:)<=ymax ...
                    & femesh.points{icmpt}(3,:)>=zmin ...
                    & femesh.points{icmpt}(3,:)<=zmax);
                
                M = btpde.massmat_cmpts{icmpt}(ll,ll);
                %volumes_voxel(icmpt,iv) = sum(sum(M));
                %clear y;
                for isimul = 1:nsimul
                    for ig=1:ng
                    y = btpde.magnetization{icmpt,isimul,1,ig}(ll);
                    %signal_voxel(icmpt,isimul,iv,ig) = sum(M*y,1);
                    Signal(ii,jj,kk,isimul,ig,icmpt)=sum(M*y,1);
                    %signal_voxel_normalize(icmpt,isimul,iv) ...
                    %     = signal_voxel(icmpt,isimul,iv)/volumes_voxel(icmpt,iv);
                    end
                end
            end
            
        end
    end
end
percent=p.stop;

Signalorig=Signal;

Signal=abs(sum(Signal,6));

Sigsave=zeros(voxels_nx-1,voxels_ny-1,voxels_nz-1,ng*(nsimul-1)+1);

for isimul=1:nsimul
    if isimul==1
        Sigsave(:,:,:,1)=squeeze(mean(Signal(:,:,:,isimul,:),5));
    else
        Sigsave(:,:,:,1+1+(isimul-2)*ng:(isimul-1)*ng+1)=squeeze(Signal(:,:,:,isimul,:));
    end
end
save_avw(abs(Sigsave),'~/GoogleDrive/Neuron.nii.gz','f',[1 1 1 1]); 


Signal=squeeze(abs(Signalorig(:,:,:,:,:,1)));
Sigsave=zeros(voxels_nx-1,voxels_ny-1,voxels_nz-1,ng*(nsimul-1)+1);

for isimul=1:nsimul
    if isimul==1
        Sigsave(:,:,:,1)=squeeze(mean(Signal(:,:,:,isimul,:),5));
    else
        Sigsave(:,:,:,1+1+(isimul-2)*ng:(isimul-1)*ng+1)=squeeze(Signal(:,:,:,isimul,:));
    end
end
save_avw(abs(Sigsave),'~/GoogleDrive/Neuron_IC.nii.gz','f',[1 1 1 1]); 

% iv=0;
% for ii = 1:voxels_nx-1
%     for jj = 1:voxels_ny-1
%         for kk = 1:voxels_nz-1
%             iv = iv+1;
%             Sigsave(ii,jj,kk,1)=mean(squeeze(signal_voxel(1,1,iv,:)));
%             Sigsave(ii,jj,kk,2:301)=squeeze(signal_voxel(1,2,iv,:));
%             Sigsave(ii,jj,kk,302:601)=squeeze(signal_voxel(1,2,iv,:));
%             Sigsave(ii,jj,kk,602:901)=squeeze(signal_voxel(1,2,iv,:));
%         end
%     end
% end
% save_avw(abs(Sigsave),'~/GoogleDrive/Neuron100.nii.gz','f',[1 1 1 1]);            
            



% Plot the finite element mesh
plot_femesh(femesh, setup.pde.compartments);
for ifig = 1:2
    figure(ifig);
    for ii = 1:voxels_nx
        xval = voxels_xlim(ii);
        pt1 = [xval,yymin,zzmin];
        pt2 = [xval,yymin,zzmax];
        pt3 = [xval,yymax,zzmax];
        pt4 = [xval,yymax,zzmin];
        face = [pt1;pt2;pt3;pt4];
        h=patch('Xdata',face(:,1),'Ydata',face(:,2),'Zdata',face(:,3));
        set(h,'facealpha',0.1);
        view([1,1,1]);
    end
    for jj = 1:voxels_ny
        yval = voxels_ylim(jj);
        pt1 = [xxmin,yval,zzmin];
        pt2 = [xxmin,yval,zzmax];
        pt3 = [xxmax,yval,zzmax];
        pt4 = [xxmax,yval,zzmin];
        face = [pt1;pt2;pt3;pt4];
        h=patch('Xdata',face(:,1),'Ydata',face(:,2),'Zdata',face(:,3));
        set(h,'facealpha',0.1);
        view([1,1,1]);
    end
    for kk = 1:voxels_nz
        zval = voxels_zlim(kk);
        pt1 = [xxmin,yymin,zval];
        pt2 = [xxmin,yymax,zval];
        pt3 = [xxmax,yymax,zval];
        pt4 = [xxmax,yymin,zval];
        face = [pt1;pt2;pt3;pt4];
        h=patch('Xdata',face(:,1),'Ydata',face(:,2),'Zdata',face(:,3));
        set(h,'facealpha',0.1);
        view([1,1,1]);
    end
%     view([0 90]);
end

% format compact;
% 
% disp('signal one voxel');
% disp(real(btpde.signal));
% disp('signal multiple voxels');
% disp(real(signal_voxel));

% signal_normalize = real(btpde.signal./volumes');
% disp('one voxel normalized');
% disp(real(btpde.signal./volumes'));
% disp('multiple voxels normalized');
% disp(real(signal_voxel_normalize));
% disp('volume fractions of geometry in the voxels');
% volfrac_voxel = full(volumes_voxel)./volumes';
% disp(volfrac_voxel);


% for icmpt = 1:ncompartment
%     figure; hold on;
%     iplot = 0;
%     for iv = 1:num_voxels
%         h = plot(real(signal_voxel_normalize(icmpt,:,iv))','o-');
%         set(h,'linewidth',1,'markersize',15);
%         iplot = iplot + 1;
%         legend_vec{iplot} = ['voxel ',num2str(iv), ', voxel vf = ', num2str(volfrac_voxel(icmpt,iv),2)];
%     end
%     h = plot(signal_normalize(icmpt,:),'x-.');
%     set(h,'linewidth',1,'markersize',15);
%     iplot = iplot + 1;
%     legend_vec{iplot} = ['voxel whole geometry'];
%     legend(legend_vec{1:iplot});
%     xlabel('simulation number');
%     ylabel('signal / volume geometry in voxel');
%     title(['Cmpt = ',num2str(icmpt)]);
%     
% end

