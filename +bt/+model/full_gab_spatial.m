classdef full_gab_spatial < bt.model.template_spatial
	% This is a template function that allows rapid prototyping of new spatial
	% fitting by leveraging model.params.spatial_spectrum with absolutely
	% no regard for performance

	properties
		% No new properties required by this derived class
        params_typical
	end

	methods
		function self = full_gab_spatial() % Constructor
			self.name = 'full_gab_spatial';
			
            self.param_names = {'Gee','Gei','Ges','Gse','Gsr','Gsn','Gre','Grs','Alpha','Beta','t0','EMGa',...
                                'Gee_Xamp','Gei_Xamp','Gsn_Xamp','t0_Xamp','Alpha_Xamp','Beta_Xamp',...
                                'Gee_Yamp','Gei_Yamp','Gsn_Yamp','t0_Yamp','Alpha_Yamp','Beta_Yamp'};                              
			self.param_symbols = {'G_{ee}','G_{ei}','G_{es}','G_{se}','G_{sr}','G_{sn}','G_{re}','G_{rs}','\alpha','\beta','t_0','A_{EMG}',...
                                'Gee_{Xamp}','Gei_{Xamp}','Gsn_{Xamp}','t_{0Xamp}','Alpha_{Xamp}','Beta_{Xamp}',...
                                'Gee_{Yamp}','Gei_{Yamp}','Gsn_{Yamp}','t_{0Yamp}','Alpha_{Yamp}','Beta_{Yamp}'};                             
			self.param_units = {'','','','','','','','','s^{-1}','s^{-1}','ms','',...
                                '','','','ms','s^{-1}','s^{-1}',   '','','','ms','s^{-1}','s^{-1}'};                            
            self.initial_step_size = [0.4    0.4     0.5     0.5    0.4    0.4    0.1     0.3      5     40    0.005  0.05];
		   	self.limits =            [eps    -20     eps     eps    -20    eps    eps     eps     10    100    0.075     0 ;...
                                      20    -eps      20      20   -eps    40     10      20    100    800     0.14     1 ];
            self.initial_step_size = [self.initial_step_size  0.1  0.1  0.2  0.001 0.5  3    0.1  0.1  0.2  0.001 0.5  3]; 
            self.limits =            [self.limits             [-4   -4  -10  -0.03 -40 -300   -4   -4  -10  -0.03 -40 -300;...
                                                                4    4   10   0.03  40  300    4    4   10   0.03  40  300]];                                                                    
                                                           
            self.n_params = length(self.param_names);                                
			self.n_fitted = self.n_params;
			self.skip_fit = zeros(1,self.n_fitted);    
			self.p = model.params();
			self.p.disable_set = true;
            self.params_typical = [];
                        
        end

        function params = fullgab2full_params(self,pars)
            params = [pars(1:2)  pars(3)*pars(4)  pars(3)*pars(5)*pars(7)  pars(5)*pars(8)  pars(9:end)];
        end

		function valid = validate_params(self,pars) % Check if parameters are valid
			% This implements a basic range check. Derived classes should impose their own checks as well
			% Probably simplest to overload this function
			valid = ~(pars(1)/pars(2) > -0.5 || (pars(1) + pars(2)) > 1 || pars(10)/pars(9) > 20 || any(pars > self.limits(2,:) | pars < self.limits(1,:) ) );
            
            % check_cosine_validity?
		end

		function set_params(self,pars)
			% This function is what decides what the spatial variations are actually going to be
			%e.g. p.apply_variation('cosine','t0',pars(10));
            self.p = p_from_fitted_params(self,self.p,pars);
            self.p.apply_variation('g_ee','cosine_2d',pars(13),pars(19));
            self.p.apply_variation('g_ei','cosine_2d',pars(14),pars(20));
            self.p.apply_variation('g_sn','cosine_2d',pars(15),pars(21));
            self.p.apply_variation('t0','cosine_2d',pars(16),pars(22));
            self.p.apply_variation('alpha','cosine_2d',pars(17),pars(23));
            self.p.apply_variation('beta','cosine_2d',pars(18),pars(24));
        end
		
		function [initial_values,prior_pp] = initialize_fit(self,target_f,target_P) % Return the first posterior/prior and initial values
			[f1,P1,idx,db_data,min_chisq] = db_fit.quick_fit(self,target_f,mean(target_P,2)); % Fit the first (or only) spectrum given
			p = model.params(db_data.iswake(idx));
            initial_values = [db_data.gab(idx,:) p.alpha(1) p.beta(1) p.t0 eps zeros(1,12)];
			prior_pp = self.uniform_priors();
		end
		
		function xyz = get_xyz(self,params) % Turn whatever model parameters are specified into XYZ - params is a matrix
            params = self.fullgab2full_params(params);
			xyz(:,1) = params(:,1)./(1-params(:,2));
			xyz(:,2) = (params(:,3)+params(:,4))./(1-params(:,5))./(1-params(:,2));
			xyz(:,3) = -params(:,5).*params(:,6).*params(:,7)./((params(:,6)+params(:,7)).^2);
        end
        
        function set_electrodes(self,EEG)

			set_electrodes@bt.model.template(self,{EEG.chanlocs.labels})
            
            %Convert the EEGLAB polar coordinates to X and Y coordinates
            %for details see electrode_positions function
            theta = 90-[EEG.chanlocs.theta]; % Because in the model framework, y=0 is the front
            radius = [EEG.chanlocs.radius];         
            self.output_x = self.p.Lx/2 * (1 + radius.*cosd(theta))';
            self.output_y = self.p.Ly/2 * (1 - radius.*sind(theta))';
            
        end
        
        function [chisq,P, chisq_eeg,chisq_lfp,P_lfp] = objective(self,pars) % Calculate the objective

            is_fit_phi_lfp = true;

            p = p_from_fitted_params(self,model.params,pars);
            [~, is_realistic] = set_nus_with_realistic_phia(p,self.params_typical.phia,0);
            if ~is_realistic
                P = [];
                chisq = NaN;
                return;
            end
            
            [chisq,P] = objective@bt.model.template_spatial(self,pars);

            if is_fit_phi_lfp && ~isnan(chisq)
                [chisq_lfp,P_lfp] = objective_lfp(self,pars);
                chisq_eeg = chisq;
                chisq = chisq_eeg + chisq_lfp;
            end
        end

        function [chisq_lfp,P_lfp] = objective_lfp(self,pars)
           
            pars = self.fullgab2full_params(pars);

            spec_w = self.target_f*2*pi;
            spec_Mtot = exp(1i*spec_w*pars(8));
		    spec_L = 1./((1-1i*spec_w/pars(6)).*(1-1i*spec_w/pars(7)));
		    Jei_oneminus = 1-spec_L.*pars(2);
		    Jsrs_oneminus = 1-spec_L.*spec_L.*pars(5);
		    re2 = self.p.re.^2;
            spec_gamma_prefactor = (1-1i*spec_w/self.p.gammae).^2;
		    q2re2 = (spec_gamma_prefactor - 1./(Jei_oneminus).*(spec_L.*pars(1) + ((spec_L.*spec_L.*pars(3) + spec_L.*spec_L.*spec_L.*pars(4)).*spec_Mtot)./(Jsrs_oneminus)));
            T_prefactor = (1-spec_L.*pars(1)-spec_L.*pars(2)).*spec_L.*self.p.phin./(Jei_oneminus.*Jsrs_oneminus); %does not include stimulation phiw
            
            Lx = 0.02; % linear dimensions of the STN/ANT
            kmax = 4;
            dk = 2*pi/Lx;
            m_rows = -kmax:kmax;
            n_cols = -kmax:kmax;
            [kxa,kya] = meshgrid(dk*m_rows,dk*n_cols);
            k2 = kxa.^2+kya.^2;
            k2u = unique(k2(:));
            k2u = [k2u histc(k2(:),k2u)];
%             k2u = [0,1];

            P_lfp = zeros(size(spec_w));
            for j = 1:size(k2u,1)
                P_lfp = P_lfp + k2u(j,2).*abs(T_prefactor./(k2u(j,1)*re2+q2re2)).^2;
            end

            P_lfp = P_lfp./utils.mex_trapz(self.target_f(self.weights>0),P_lfp(self.weights>0));
            normalization_target_lfp = utils.mex_trapz(self.target_f(self.weights>0),self.target_P_lfp(self.weights>0));
            P_lfp = normalization_target_lfp*P_lfp;
            sqdiff_lfp = (abs(P_lfp-self.target_P_lfp)./self.target_P_lfp).^2;
            chisq_lfp = sum(sqdiff_lfp(:).*self.weights(:));

        end
		
	end
end



