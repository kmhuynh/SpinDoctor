function [SOL,SIG_cmpts,SIG_allcmpts,difftime,ctime] ...
    = SIG_BTPDE(experiment,mymesh,DIFF_cmpts,kappa_bdys,IC_cmpts,OUTPUT_MAGNETIZATION)

% solve Bloch-Torrey equation
% 
% Input:
%     1. experiment is a structure with 10 elements:
%         ngdir_total 
%         gdir        
%         sdeltavec   
%         bdeltavec   
%         seqvec       
%         npervec     
%         rtol        
%         atol        
%         qvalues     
%         bvalues        
%     2. mymesh is a structure with 10 elements:
%         Nnode
%         Nele
%         Nface
%         Pts_cmpt_reorder
%         Ele_cmpt_reorder
%         Pts_ind
%         Pts_boundary_reorder
%         Fac_boundary_reorder
%         Nboundary
%         Ncmpt
%     3. DIFF_cmpts
%     4. kappa_bdys
%     5. IC_cmpts
%     6. OUTPUT_MAGNETIZATION
%
% Output:
%     1. SOL
%     2. SIG_cmpts
%     3. SIG_allcmpts
%     4. difftime
%     5. ctime

global FEM_M FEM_K FEM_A FEM_Q FEM_G
global QVAL UG 
global BDELTA SDELTA SEQ OGSEPER 

bvalues = experiment.bvalues;
qvalues = experiment.qvalues;
gdir = experiment.gdir;
sdeltavec = experiment.sdeltavec;
bdeltavec = experiment.bdeltavec;
seqvec = experiment.seqvec;
npervec = experiment.npervec;
ODEsolve_rtol = experiment.rtol;
ODEsolve_atol = experiment.atol;
yes = 1;  no = 0;

UG = gdir;
UG = UG/norm(UG);

disp(['Solving Bloch-Torrey PDE, direction: ',num2str(UG(1)),',',num2str(UG(2)),',',num2str(UG(3))]);

if (mymesh.Ncmpt == 1 || abs(max(kappa_bdys)) <= 1e-16)
    DO_COUPLING = no;
else
    DO_COUPLING = yes;
end

nexperi = length(sdeltavec);
%disp(['Simulating ',num2str(nexperi),' experiments']);
nb = size(qvalues,2);

ctime = nan*ones(nexperi,nb);

%disp('Setting up FEM matrices');
%tic
for icmpt = 1:mymesh.Ncmpt
    coordinates = mymesh.Pts_cmpt_reorder{icmpt};
    elements = mymesh.Ele_cmpt_reorder{icmpt};
    facets = [];
    GX = -sqrt(-1)*UG*coordinates;
    FEM_MAT{icmpt}.Q = sparse(length(coordinates),length(coordinates));
    for iboundary = 1:mymesh.Nboundary
        if (kappa_bdys(iboundary) ~= 0)
            neumann = mymesh.Fac_boundary_reorder{icmpt}{iboundary}';
            neumann_nodes = unique(neumann);
            %zvalues = mymesh.Pts_cmpt_reorder{1}(3,neumann_nodes);
            %mvalues = sqrt((20-abs(zvalues)));
            coeffs_flux_matrix = zeros(max(neumann_nodes),1);
            coeffs_flux_matrix(neumann_nodes) = kappa_bdys(iboundary); %*mvalues
            if ~isempty(neumann)
                FEM_MAT{icmpt}.Q = FEM_MAT{icmpt}.Q + flux_matrixP1_3D(neumann,coordinates',coeffs_flux_matrix);
            end
        end
    end
    [FEM_MAT{icmpt}.K,volumes] = stiffness_matrixP1_3D(elements',coordinates',DIFF_cmpts(icmpt));
    FEM_MAT{icmpt}.M = mass_matrixP1_3D(elements',volumes);
    FEM_MAT{icmpt}.A = mass_matrixP1_3D(elements',volumes,GX');    
    % calculate the IC
    IC_Pts{icmpt} = IC_cmpts(icmpt)*ones(size(mymesh.Pts_cmpt_reorder{icmpt},2),1);
end
%toc


if (DO_COUPLING == yes)
    disp('Coupling FEM matrices');
    tic
    [FEMcouple_MAT,FEMcouple_ind0,FEMcouple_indf] ...
        = generate_FEM_coupling(FEM_MAT,mymesh.Ncmpt,mymesh.Nboundary,...
        mymesh.Pts_cmpt_reorder,mymesh.Ele_cmpt_reorder,mymesh.Pts_ind,mymesh.Pts_boundary_reorder,mymesh.Fac_boundary_reorder);
    toc
    % IC for coupled case
    IC_couple = zeros(size(FEMcouple_MAT.M,1),1);
    for icmpt = 1:mymesh.Ncmpt
        IC_couple(FEMcouple_ind0(icmpt):FEMcouple_indf(icmpt),1) = IC_Pts{icmpt};
    end
else
    FEMcouple_MAT = [];
    FEMcouple_ind0 = [];
    FEMcouple_indf = [];    
end

%% solve ODE
for iexperi = 1:nexperi
    SDELTA = sdeltavec(iexperi);
    BDELTA = bdeltavec(iexperi);
    TE = SDELTA+BDELTA;
    SEQ = seqvec(iexperi);% for choosing case PGSE, OGSEcos or OGSEsin
    omega = 2*pi*npervec(iexperi)/SDELTA;
    OGSEPER = 1./omega*2*pi;%% set up number for OGSE   
    
    disp(['Experiment: sdelta ',num2str(SDELTA), ' bdelta ',num2str(BDELTA)]);
    
    TLIST = [0,TE];
    
    for ib = 1:nb     
        
        tic;        
        
        b_start_time = clock;
        
        % global variable setting QVAL for ODE time stepping
        
        QVAL = qvalues(iexperi,ib);               
		disp(['            qvalue ',num2str(QVAL,'%.1e'), ' bvalue ', num2str(bvalues(iexperi,ib),'%.1e')]);        
        
        difftime(iexperi) = seqdifftime;
        
        %% Solving for case of coupling between compartments.      
        if (DO_COUPLING == yes)
            FEM_M = FEMcouple_MAT.M;
            FEM_K = FEMcouple_MAT.K;
            FEM_A = FEMcouple_MAT.A; %*seqprofile(state.time)*QVAL;
            FEM_Q = FEMcouple_MAT.Q;
            FEM_G = sparse(zeros(size(FEM_M,1),1));
            
            options = odeset('Mass',FEM_M,'AbsTol',ODEsolve_atol,'RelTol',ODEsolve_rtol,'Vectorized','on','Stats','off',...
                'Jacobian',@odejac_bt_includeb);            
            disp('***Coupled: start ode solve ode15s'); tic
            sol = ode15s(@odefun_bt_includeb,TLIST,IC_couple,options);
            disp('***Coupled: end ode solve ode15s'); toc
            for icmpt = 1:mymesh.Ncmpt
                YOUT{iexperi}{ib}{icmpt} = sol.y(FEMcouple_ind0(icmpt):FEMcouple_indf(icmpt),:);
                TOUT{iexperi}{ib}{icmpt} = sol.x;
                MT{iexperi}{ib}{icmpt} = sum(FEM_MAT{icmpt}.M*YOUT{iexperi}{ib}{icmpt},1);
                if (OUTPUT_MAGNETIZATION)
                    SOL{iexperi}{ib}{icmpt} = YOUT{iexperi}{ib}{icmpt}(:,end);
                else
                    SOL{iexperi}{ib}{icmpt} = [];
                end
                YOUT{iexperi}{ib}{icmpt} = [];
                TOUT{iexperi}{ib}{icmpt} = [];
            end
        else
            %% Solving for case of no coupling between compartments.            
            for icmpt = 1:mymesh.Ncmpt
                FEM_M = FEM_MAT{icmpt}.M;
                FEM_K = FEM_MAT{icmpt}.K;
                FEM_A = FEM_MAT{icmpt}.A;
                FEM_Q = FEM_MAT{icmpt}.Q;
                FEM_G = sparse(zeros(size(FEM_M,1),1));
                
                options = odeset('Mass',FEM_M,'AbsTol',ODEsolve_atol,'RelTol',ODEsolve_rtol,'Vectorized','on','Stats','off',...
                    'Jacobian',@odejac_bt_includeb);
                %disp('***Uncoupled: start ode solver ode15s'); tic
                ICC = IC_Pts{icmpt};
                if (max(abs(ICC))<=1e-16)
                    sol.y = zeros(size(ICC,1),size(TLIST,2));
                    sol.x = TLIST;
                else
                    sol = ode15s(@odefun_bt_includeb,TLIST,ICC,options);
                end
                %disp('***Uncoupled: end ode solver ode15s'); toc
                YOUT{iexperi}{ib}{icmpt} = sol.y;
                TOUT{iexperi}{ib}{icmpt} = sol.x;
                MT{iexperi}{ib}{icmpt} = sum(FEM_MAT{icmpt}.M*YOUT{iexperi}{ib}{icmpt},1);  
                
                if (OUTPUT_MAGNETIZATION)
                    SOL{iexperi}{ib}{icmpt} = YOUT{iexperi}{ib}{icmpt}(:,end);
                else
                    SOL{iexperi}{ib}{icmpt} = [];
                end
                
                YOUT{iexperi}{ib}{icmpt} = [];
                TOUT{iexperi}{ib}{icmpt} = [];
            end            
        end
        
        ctime(iexperi,ib)=etime(clock, b_start_time);
        
        toc 
    end
end

SIG_cmpts = zeros(mymesh.Ncmpt, nexperi, nb);
SIG_allcmpts = zeros(nexperi, nb);
for iexperi = 1:nexperi 
    for ib = 1:nb
        for icmpt = 1:mymesh.Ncmpt
            SIG_cmpts(icmpt,iexperi,ib) = MT{iexperi}{ib}{icmpt}(end);
            SIG_allcmpts(iexperi,ib) = SIG_allcmpts(iexperi,ib) + SIG_cmpts(icmpt,iexperi,ib);
        end
    end
end