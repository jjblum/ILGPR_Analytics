classdef LGP < handle
    
    %{
        Represents one gaussian process.
        Uses GPML library
    %}
    
    properties
        X        % data locations matrix
        Z        % data value vector
        W        % length scale matrix
        hyp      % vector of log(hyperparameters)
        u        % center location
        L        % lower triangular (Cholesky) decomp of K
        R        % permutation matrix
        K        % covariance matrix
        invK     % inverse of covariance matrix + sensor noise on diagonal
        N        % number of training points
        NMAX     % maximum number of training points
        NSTART   % number of points needed before computation actually begins
        started  % boolean, false while N < NSTART
        data     % vector of Datum objects
        f        % most recent likelihood result
        ID       % the index of the LGP itself (order of creation)
        alpha    % the prediction vector = inv(K + sigma_n^2*I)*Z
    end
    
    methods
        function obj = LGP(ID,Datum)            
            obj.ID = ID;
            obj.X = Datum.getX();
            obj.Z = Datum.getZ();
            obj.hyp.mean = 0;
            obj.hyp.cov = [log(5) log(5) log(1)]'; % log(length scale 1), log(length scale 2), log(process variability)
            obj.hyp.lik = log(exp(1)); % sigma_n, noise in sensor
            obj.W = diag((1./exp(obj.hyp.cov(1:2))).^2);
            obj.u = Datum.getX();
            obj.R = [];
            obj.K = [];
            obj.L = [];
            obj.N = 1;
            obj.NMAX = 100;
            obj.NSTART = 10;
            obj.started = 0;
            obj.data = cell(1,1);
            obj.data{1} = Datum;
            obj.f = 0;
        end
        
        function newDatum(self,Datum)
            self.N = self.N + 1;
            self.data{self.N} = Datum;
            updateX(self,Datum.getX());
            updateZ(self,Datum.getZ());
            updateCenter(self);
            choleskyUpdate(self);
        end
        
        function [zHat, sHat] = predict(self,x)
            Xx = horzcat(self.X,x);
            a = bsxfun(@minus,Xx,mean(Xx,2));
            oldX = a(:,1:end-1);
            newX = a(:,end);
            K_new = self.covarianceVector(oldX,newX);
            zHat = K_new'*self.alpha;                
            
            k_new = exp(self.hyp.cov(3))^2 + exp(self.hyp.lik) + 1e-5;
%             KinvRowsVector = kron(eye(length(self.alpha)),self.Z')\self.alpha;
%             KinvMatrix = reshape(KinvRowsVector,length(self.alpha),length(self.alpha))'; % note the last transpose
            % this KinvMatrix is NOT unique. So I don't know if we can actually use it. Probably not.
            KinvMatrixByChol = self.L'\inv(self.L);           
            sHat = k_new - K_new'*KinvMatrixByChol*K_new;        
        end        
        
        function optimizeHyperparameters(self)  
            covfunc = @covSEard;
            meanfunc = @meanConst; 
            likfunc = @likGauss;
            self.hyp.mean = mean(self.Z); % guess that mean function is just the mean of the training data
            [self.hyp, f_log, iterations] = minimize(self.hyp, @gp, -10, @infExact, meanfunc, covfunc, likfunc, self.X', self.Z);
            self.W = diag((1./exp(self.hyp.cov(1:2))).^2);
            firstCholesky(self);
        end        
        
        function updateX(self, x)
            self.X = horzcat(self.X,x);            
        end
        
        function updateZ(self, z)
            self.Z = vertcat(self.Z,z);
        end
        
        function getR(self, m)
            % m = index of value being replaced in covariance matrix
        end
        
        function updateCenter(self)
            self.u = mean(self.X,2);
        end
        
        function choleskyUpdate(self)
            
            if self.N < self.NSTART
                return;
            elseif self.N == self.NSTART
                firstCholesky(self);    
                self.started = 1;
%             elseif self.N > self.NMAX
                % the incremental cholesky update with removal (Find R and include it)
            else
                % the incremental cholesky update without removal
                
                % update L
                self.incrementalCholeskyUpdate();                
                
            end
            self.optimizeHyperparameters();
        end                
        
        function K_new = covarianceVector(self,oldX,newX)
            oldX_minus_newX = bsxfun(@minus,oldX,newX);
            K_new = exp(self.hyp.cov(3))^2*exp(-1/2*sum(oldX_minus_newX'*self.W.*oldX_minus_newX',2));            
        end
        
        function incrementalCholeskyUpdate(self)

            a = bsxfun(@minus,self.X,mean(self.X,2));
            oldX = a(:,1:end-1);
            newX = a(:,end);
            K_new = self.covarianceVector(oldX,newX);
            k_new = exp(self.hyp.cov(3))^2 + exp(self.hyp.lik) + 1e-5;
            l = self.forwardSubstitution(self.L,K_new);
            lstar = sqrt(k_new - norm(l,2)^2);
            
            if ~isreal(lstar)
                disp(sprintf('WARNING: LGP %d l* element is imaginary\n',self.ID));
                keyboard                
            end
            
            self.L = vertcat(horzcat(self.L,zeros(self.N-1,1)),horzcat(l',lstar));            
            
            self.K = self.L*self.L';
            
            if ~isreal(self.L)
                disp(sprintf('WARNING: LGP %d L matrix has imaginary elements\n',self.ID));
                keyboard
            end
                
            if (rcond(self.L) < 1e-12 )
                disp(sprintf('WARNING: LGP %d L matrix is singular\n',self.ID));
                keyboard
            end
            if (rcond(self.K) < 1e-12 )
                disp(sprintf('WARNING: LGP %d K matrix is singular\n',self.ID));
                keyboard
            end            
            
            self.alpha = self.backwardSubstitution(self.L',self.forwardSubstitution(self.L,self.Z));
        end
        
        function firstCholesky(self)
            firstK(self);
            self.L = chol(self.K,'lower');
            self.alpha = self.K\self.Z;
        end
        
        function firstK(self)
            self.K = covSEard(self.hyp.cov,self.X');
            self.K = self.K + (exp(self.hyp.lik) + 1e-5)*eye(size(self.K));
        end
        
        function x = forwardSubstitution(self,L,b)
            n = length(b);
            x = zeros(n,1);
            x(1) = b(1)/L(1,1);
            for j = 2:n
                x(j) = (b(j)-L(j,1:j-1)*x(1:j-1))/L(j,j);
            end
        end
        
        function x = backwardSubstitution(self,U,b)
            n = length(b);
            x = zeros(n,1);
            x(n) = b(n)/U(n,n);
            for j = n-1:-1:1
                x(j) = (b(j)-U(j,j+1:n)*x(j+1:n))/U(j,j);
            end
        end        
        
        function f = getF(self)
            f = self.f;
        end
        
        function updateF(self,x)
            %             keyboard
            % ratio = norm(x-self.u)*sqrt(self.W(1,1))
            % rough_f = exp(-1/2*ratio)
            self.f = exp(-1/2*(x-self.u)'*self.W*(x-self.u));
            % approximately exp(-1/2*[ratio of squared distance to squared length scale]). e.g. if ratio is 5, f ~ 0.1
            % So if we want a new LGP when more than 2*length scale away, need to set cutoff to about 0.4
        end
        
        function started = isStarted(self)
            started = self.started;
        end
    end
    
end