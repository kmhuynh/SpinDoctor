function [SIG_cmpts,mean_SIG_cmpts,ADC_cmpts]=simAxonSig(rcell,kcell,bvecs,blist,smalldelta,bigdelta,tempdir)

% set up structural files
setup=setup_1axon_custom(rcell,kcell,blist,smalldelta,bigdelta,bvecs,tempdir)

% Set up the PDE model in the geometrical compartments.
setup.pde = prepare_pde(setup);

% Prepare experiments (gradient sequence, bvalues, qvalues, solvers)
setup = prepare_experiments(setup);

% Create or load finite element mesh
[femesh, surfaces, cells] = create_geometry(setup);

% Get volume and surface area quantities from mesh
[volumes, surface_areas] = get_vol_sa(femesh);

% Compute volume weighted mean of diffusivities over compartments
% Take trace of each diffusion tensor, divide by 3
mean_diffusivity = trace(sum(setup.pde.diffusivity .* shiftdim(volumes, -1), 3)) ...
    / (3 * sum(volumes));

% Get sizes
ncompartment = length(setup.pde.compartments);
namplitude = length(setup.gradient.values);
nsequence = length(setup.gradient.sequences);
ndirection = size(setup.gradient.directions, 2);

% %% Perform HADC experiment
% if isfield(setup, "hadc")
%     % Solve HADC model
%     hadc = solve_hadc(femesh, setup);
% end
% 
% ADC_allcmpts=hadc.adc_allcmpts;

%% Perform MF experiments
if isfield(setup, "mf")
    % Perform Laplace eigendecomposition
    eiglim = length2eig(setup.mf.length_scale, mean_diffusivity);
    lap_eig = compute_laplace_eig(femesh, setup.pde, eiglim, setup.mf.neig_max);

    % Compute length scales of eigenvalues
    lap_eig.length_scales = eig2length(lap_eig.values, mean_diffusivity);

%     % Compute the JN value that relates the eigenmodes to their contribution
%     % to the Matrix Formalism signal for a diffusion-encoding sequence
%     mf_jn = compute_mf_jn(lap_eig.values, setup);

    % Compute MF magnetization
    mf = solve_mf(femesh, setup, lap_eig);

    % Fit ADC from MF signal
    mf_fit = fit_signal(mf.signal, mf.signal_allcmpts, setup.gradient.bvalues);

end
ADC_cmpts=squeeze(mf_fit.adc);
mf.signal=abs(mf.signal);
mf.signal=permute(mf.signal,[4 1 2 3]);
SIG_cmpts=zeros(1+300*numel(find(setup.gradient.values>0)),ncompartment);
mean_SIG_cmpts=zeros(length(setup.gradient.values),ncompartment);
for i=1:ncompartment
    for ib=1:length(setup.gradient.values)
        if setup.gradient.values(ib)==0
           SIG_cmpts(1,i)=mean(mf.signal(:,i,ib),1);
           mean_SIG_cmpts(ib,i)=SIG_cmpts(1,i);
        else            
           SIG_cmpts(1+[1:ndirection]+(find(setup.gradient.values==setup.gradient.values(ib))-2)*length(ndirection),i)=mf.signal(:,i,find(setup.gradient.values==setup.gradient.values(ib)));
           mean_SIG_cmpts(ib,i)=mean(SIG_cmpts(1+[1:ndirection]+(find(setup.gradient.values==setup.gradient.values(ib))-2)*length(ndirection),i));
        end
    end
end
mean_SIG_cmpts=mean_SIG_cmpts./max(mean_SIG_cmpts);
rmdir([setup.name '_dir'],'s');
delete([setup.name '_cells']);

